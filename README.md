# openfe-scripts

Small helper scripts that extend the standard OpenFE CLI workflow with a few extra utilities for setting up and analysing OpenFE calculations.

Current scripts include:

- `prep_rbfe_hybridtop.py` for preparing hybrid-topology RBFE networks and writing transformation JSON files
- `prep_rbfe_septop.py` for preparing separated-topology (SepTop) RBFE networks and writing transformation JSON files
- `workup.py` for summarising and analysing completed OpenFE calculations
- `plot_network.py` for plotting ligand-network summaries with RBFE results

## Installation

There is no special installation for this repository.

The scripts are intended to be run from a vanilla OpenFE conda/mamba environment with the usual OpenFE dependencies available, for example:

```bash
micromamba activate openfe_env
```

or the equivalent conda environment you already use for OpenFE.

## `prep_rbfe_hybridtop.py`

This script prepares an RBFE hybrid-topology setup from:

- a receptor PDB file
- a ligand SDF file
- a selected network-planning mode

It writes:

- `network_setup/ligand_network.png`
- `network_setup/ligand_network.graphml`
- `network_setup/transformations/*.json`
- mapping PNG files into `mappings/`

### Basic usage

Minimal example using the default settings:

```bash
python3 prep_rbfe_hybridtop.py \
  --rec example_small/structures/5H9P_prepped.pdb \
  --ligs example_small/structures/thiodigalactoside_thiolactose_Cfiller.sdf \
  --mapper kartograf \
  --scorer lomap \
  --network minimal_redundant \
  --windows 11 \
  --windowtime 5 \
  --output-dir .
```

Important defaults:

- `--network minimal_redundant`
- `--mapper kartograf`
- `--scorer lomap`
- `--windows 11`
- `--windowtime 5`

Argument meanings:

- `--rec`: path to the prepared receptor PDB file
- `--ligs`: path to the ligand SDF file
- `--mapper`: atom mapper used to generate atom mappings between ligands
- `--scorer`: scorer used for automatic network planning; ignored for `--network custom`
- `--network`: network-construction mode, one of `minimal_spanning`, `minimal_redundant`, `radial`, `maximal`, or `custom`
- `--windows`: number of lambda windows and OpenMM replicas
- `--windowtime`: production simulation time per window in nanoseconds
- `--output-dir`: directory where `network_setup/` and `mappings/` will be written

Additional optional arguments:

- `--central-ligand`: central ligand name for `--network radial`
- `--custom-network`: YAML file with custom ligand pairs for `--network custom`

### Custom network from ligand names

You can also define the network explicitly with a YAML file and run:

```bash
python3 prep_rbfe_hybridtop.py \
  --rec example_small/structures/5H9P_prepped.pdb \
  --ligs example_small/structures/thiodigalactoside_thiolactose_Cfiller.sdf \
  --network custom \
  --custom-network custom_network.yaml
```

An example file is included in this repository as `custom_network.yaml`.

The YAML file must contain a top-level `edges` key with ligand-name pairs:

```yaml
edges:
  - [thiodigalactoside, thiolactose]
  - [thiolactose, Cfiller]
```

Ligand names in the YAML must match the molecule names in the input SDF.

### Help

To see all available options:

```bash
python3 prep_rbfe_hybridtop.py --help
```

## `prep_rbfe_septop.py`

This script prepares an RBFE separated-topology (SepTop) setup from:

- a receptor PDB file
- a ligand SDF file
- a selected network-planning mode

It writes:

- `network_setup/ligand_network.png`
- `network_setup/ligand_network.graphml`
- `network_setup/transformations/*.json`
- mapping PNG files into `mappings/`

Unlike the hybrid-topology script, each SepTop transformation JSON represents the full RBFE cycle for one ligand pair. The solvent and complex legs are handled internally by the SepTop protocol, so the script writes one transformation per ligand-network edge with `mapping=None`.

### Basic usage

Minimal example using the default network-planning settings:

```bash
python3 prep_rbfe_septop.py \
  --rec example_small/structures/5H9P_prepped.pdb \
  --ligs example_small/structures/thiodigalactoside_thiolactose_Cfiller.sdf \
  --mapper kartograf \
  --scorer lomap \
  --network minimal_redundant \
  --protocol-repeats 1 \
  --host-min-distance 0.5 \
  --host-max-distance 1.5 \
  --output-dir .
```

Useful optional arguments:

- `--windows`: resize the default SepTop lambda schedules and set the replica count
- `--windowtime`: set the production length per lambda window in nanoseconds
- `--equilibration-time`: set the alchemical equilibration length in nanoseconds
- `--central-ligand`: central ligand name for `--network radial`
- `--custom-network`: YAML file with custom ligand pairs for `--network custom`

### Analysis note

SepTop result analysis is best kept in a separate follow-on workflow. The existing `workup.py` is built around hybrid-topology solvent/complex leg JSONs and does not yet understand SepTop result structure.

For available options:

```bash
python3 prep_rbfe_septop.py --help
```
