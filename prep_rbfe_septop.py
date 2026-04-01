import openfe
import pathlib
from openfe.protocols.openmm_septop import SepTopProtocol
from openff.units import unit

receptor_path = "../structures/5H9P_prepped.pdb"
ligands_path = "../structures/lactoamides_keyleg.sdf"

mapper = "lomap"  # options: "lomap" or "kartograf"
scorer = "lomap"    # options: "lomap" or "kartograf_rmsd"
Nrepeats = 1
N_lambda_windows = 22 




# Load the ligands from an SDF file
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
network_planner = openfe.ligand_network_planning.generate_minimal_spanning_network

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
    




# Creating chemical systems
solvent = openfe.SolventComponent()
protein = openfe.ProteinComponent.from_pdb_file(receptor_path)
systemA = openfe.ChemicalSystem({
    'ligand': charged_ligands[0],
    'solvent': solvent,
    'protein': protein
})
systemB = openfe.ChemicalSystem({
    'ligand': charged_ligands[1],
    'solvent': solvent,
    'protein': protein
})




# Creating settings
settings = SepTopProtocol.default_settings()
# Run only a single repeat
settings.protocol_repeats = 1
# Change the min and max distance between protein and ligand atoms for Boresch restraints to avoid periodicity issues
settings.complex_restraint_settings.host_min_distance = 0.5 * unit.nanometer
settings.complex_restraint_settings.host_max_distance = 1.5 * unit.nanometer
# Set the equilibration time to 2 ns (which is also the default)
settings.solvent_simulation_settings.equilibration_length = 2000 * unit.picosecond
settings.complex_simulation_settings.equilibration_length = 2000 * unit.picosecond

print("settings:", settings)

# Creating a protocol
protocol = SepTopProtocol(settings)


# Creating the alchemical network
transformations = []
for edge in ligand_network.edges:
    # use the solvent and protein created above
    sysA_dict = {'ligand': edge.componentA,
                 'protein': protein,
                 'solvent': solvent}
    sysB_dict = {'ligand': edge.componentB,
                 'protein': protein,
                 'solvent': solvent}

    # we don't have to name objects, but it can make things (like filenames) more convenient
    sysA = openfe.ChemicalSystem(sysA_dict, name=f"{edge.componentA.name}")
    sysB = openfe.ChemicalSystem(sysB_dict, name=f"{edge.componentB.name}")

    transformation = openfe.Transformation(
        stateA=sysA,
        stateB=sysB,
        mapping=None,
        protocol=protocol,  # use protocol created above
        name=f"rbfe_{sysA.name}_{sysB.name}"
    )
    transformations.append(transformation)

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

    charge_difference = mapping.get_alchemical_charge_difference()
    print("Charge difference:", charge_difference)


    time = "160"
    partition = "b32_128_gpu_extra"

    for repeat in range(1, Nrepeats + 1):
        print("repeat", repeat)

        # Variables to insert
        job_name = f"rbfe_{mapping.componentA.name}_{mapping.componentB.name}"
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
