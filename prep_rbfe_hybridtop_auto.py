import openfe
import pathlib
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from openff.units import unit
import copy


receptor_path = "../structures/5H9P_prepped.pdb"
ligands_path = "../structures/lactoamides.sdf"

mapper = "kartograf"  # options: "lomap" or "kartograf"
scorer = "lomap"    # options: "lomap" or "kartograf_rmsd"
allow_partial_fused_rings = True
Nrepeats = 3



def get_Nperturbed(mapping) -> int:
    
    Nperturbed_A = len(list(mapping.componentA_unique))
    Nperturbed_B = len(list(mapping.componentB_unique))
    Nperturbed = Nperturbed_A + Nperturbed_B
    
    return Nperturbed



# Loading the ligands
from rdkit import Chem
print("Loading ligands from SDF...")
supp = Chem.SDMolSupplier(ligands_path, removeHs=False)
ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in supp]


# Charging the ligands
print("Assigning partial charges to ligands...")
from openfe.protocols.openmm_utils.omm_settings import OpenFFPartialChargeSettings
from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges

charge_settings = OpenFFPartialChargeSettings(partial_charge_method="am1bcc", off_toolkit_backend="ambertools")

charged_ligands = bulk_assign_partial_charges(
    molecules=ligands,
    overwrite=False,
    method=charge_settings.partial_charge_method,
    toolkit_backend=charge_settings.off_toolkit_backend,
    generate_n_conformers=charge_settings.number_of_conformers,
    nagl_model=charge_settings.nagl_model,
    processors=1
)


# Define mapper
if mapper == "kartograf":
    from kartograf import KartografAtomMapper
    from kartograf.atom_mapping_scorer import MappingRMSDScorer
    mapper = KartografAtomMapper(atom_map_hydrogens=True, allow_partial_fused_rings=allow_partial_fused_rings)
if mapper == "lomap":
    mapper = openfe.LomapAtomMapper(max3d=1.0, element_change=False)

# Define scorer
if scorer == "kartograf_rmsd":
    scorer = MappingRMSDScorer()
if scorer == "lomap":
    scorer = openfe.lomap_scorers.default_lomap_score




# Creating the ligand network
print("Creating ligand network...")
network_planner = openfe.ligand_network_planning.generate_minimal_redundant_network
ligand_network = network_planner(
    ligands=charged_ligands,
    mappers=[mapper],
    scorer=scorer
)

# Save ligand network into picture
from openfe.utils.atommapping_network_plotting import plot_atommapping_network
network_fig = plot_atommapping_network(ligand_network)
network_fig.savefig("ligand_network.png", dpi=300, bbox_inches="tight")


# Writing the AlchemicalNetwork to disk
print("Writing mappings to disk...")
mappings_dir = pathlib.Path("mappings")
mappings_dir.mkdir(exist_ok=True)

# Write the graphml file
with open("ligand_network.graphml", mode='w') as f:
    f.write(ligand_network.to_graphml())






# Defining settings
# Settings for non-compilated edges (default with reduced solvent padding)
settings_normal = RelativeHybridTopologyProtocol.default_settings()
settings_normal.solvation_settings.solvent_padding = 2 * unit.nanometer
settings_normal.protocol_repeats = 1

# Settings for edges where the charge changes (20 ns, alchemical ion, 22 windows)
settings_charge_change = RelativeHybridTopologyProtocol.default_settings()
settings_charge_change.solvation_settings.solvent_padding = 2 * unit.nanometer
settings_charge_change.alchemical_settings.explicit_charge_correction = True
settings_charge_change.simulation_settings.production_length = 20 * unit.nanosecond
settings_charge_change.simulation_settings.n_replicas = 22
settings_charge_change.lambda_settings.lambda_windows = 22
settings_charge_change.protocol_repeats = 1

# Settings for edges with a large perturbed region (10 ns, 22 windows) - not used currently
settings_large = RelativeHybridTopologyProtocol.default_settings()
settings_large.solvation_settings.solvent_padding = 2 * unit.nanometer
settings_large.simulation_settings.production_length = 10 * unit.nanosecond
settings_large.simulation_settings.n_replicas = 22
settings_large.lambda_settings.lambda_windows = 22
settings_large.protocol_repeats = 1









# Defining solvent and protein components
solvent = openfe.SolventComponent()
protein = openfe.ProteinComponent.from_pdb_file(receptor_path)


# Creating the alchemical network
print("Creating alchemical network...")
transformations = []
for mapping in ligand_network.edges:
    print(mapping.componentA.name, "->", mapping.componentB.name)
    print("Score:", mapping.annotations)
    
    Nperturbed = get_Nperturbed(mapping)
    print("Nperturbed", Nperturbed)
    charge_difference = mapping.get_alchemical_charge_difference()
    print("Charge difference:", charge_difference)
    
    # Save mapping figure
    mapping.draw_to_file(f"mappings/mapping_{mapping.componentA.name}_{mapping.componentB.name}.png") 
    

    # Choose settings based on the type of transformation        
    if abs(charge_difference) > 1e-3:
        print("Using charge change settings")
        settings = copy.deepcopy(settings_charge_change)
    elif Nperturbed > 15:
        print("Using large perturbed region settings")
        settings = copy.deepcopy(settings_large)
    else:
        print("Using normal settings")
        settings = copy.deepcopy(settings_normal)                  

         
    # Creating chemical systems
    for leg in ['solvent', 'complex']:
        print("leg", leg)
        
        # use the solvent and protein created above
        sysA_dict = {'ligand': mapping.componentA,
                     'solvent': solvent}
        sysB_dict = {'ligand': mapping.componentB,
                     'solvent': solvent}

        if leg == 'complex':
            # If this is a complex transformation we add in the protein to the chemical states, and decrease solvent padding
            sysA_dict['protein'] = protein
            sysB_dict['protein'] = protein
            settings.solvation_settings.solvent_padding = 1 * unit.nanometer
        else:
            settings.solvation_settings.solvent_padding = 2 * unit.nanometer
        
        print("padding settings:", settings.solvation_settings.solvent_padding)

        # Create protocol using the settings defined above
        protocol = RelativeHybridTopologyProtocol(settings)
        
        # we don't have to name objects, but it can make things (like filenames) more convenient
        sysA = openfe.ChemicalSystem(sysA_dict, name=f"{mapping.componentA.name}_{leg}")
        sysB = openfe.ChemicalSystem(sysB_dict, name=f"{mapping.componentB.name}_{leg}")

        # Create transformation
        transformation = openfe.Transformation(
            stateA=sysA,
            stateB=sysB,
            mapping=mapping,
            protocol=protocol,
            name=f"rbfe_{sysA.name}_{sysB.name}"
        )
        transformations.append(transformation)
        
        print()

network = openfe.AlchemicalNetwork(transformations)


# Writing the AlchemicalNetwork to disk
print("Writing transformations to disk...")
transformation_dir = pathlib.Path("transformations")
transformation_dir.mkdir(exist_ok=True)

# Write out each transformation
for transformation in network.edges:
    transformation.to_json(transformation_dir / f"{transformation.name}.json")





# Writing the slurm submission scripts
for mapping in ligand_network.edges:
    print(mapping.componentA.name, "->", mapping.componentB.name)
    Nperturbed = get_Nperturbed(mapping)
    charge_difference = mapping.get_alchemical_charge_difference()



    for leg in ['solvent', 'complex']:
        print("leg", leg)

        if leg == 'solvent':
            if abs(charge_difference) > 1e-3 or Nperturbed > 15:
                time = "12"
                partition = "b32_128_gpu_long"
            else:
                time = "04"
                partition = "b32_128_gpu"
        
        if leg == 'complex':
            if abs(charge_difference) > 1e-3 or Nperturbed > 15:
                time = "160"
                partition = "b32_128_gpu_extra"
            else:
                time = "04"
                partition = "b32_128_gpu"

        for repeat in range(1, Nrepeats + 1):
            print("repeat", repeat)

            # Variables to insert
            job_name = f"rbfe_{mapping.componentA.name}_{leg}_{mapping.componentB.name}_{leg}"
            json_in = f"transformations/{job_name}.json"
            json_out = f"results/repeat{repeat}/{job_name}.json"

            # Write to file
            filename = f"{job_name}_{repeat}.sh"
            with open(filename, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"#SBATCH -J {job_name}\n")
                f.write(f"#SBATCH --time={time}:00:00 --mem=55G --partition={partition}\n")
                f.write("#SBATCH --gres=gpu:1\n")
                f.write("#SBATCH --nodes=1\n")
                f.write("#SBATCH --mail-type=NONE\n")
                f.write("eval \"$(micromamba shell hook --shell=bash)\"\n")
                f.write("micromamba activate openfe_env\n")
                f.write(f"openfe quickrun {json_in} -o {json_out} -d results/repeat{repeat}/{job_name} > {job_name}_{repeat}.log")
