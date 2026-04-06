# openfe-scripts

A collection of helper scripts that extend the [OpenFE](https://github.com/OpenFreeEnergy/openfe) CLI workflow with extra utilities for setting up, running, and analyzing Relative Free Energy (RFE) calculations.

## Current Scripts

- **`prep_rbfe_hybridtop.py`**: Prepare hybrid-topology RBFE networks and write transformation JSON files.
- **`prep_rbfe_septop.py`**: Prepare separated-topology (SepTop) RBFE networks and write transformation JSON files.
- **`workup_hybridtop.py`**: Summarize and analyze completed OpenFE calculations, generating CSV reports and convergence plots.
- **`plot_network.py`**: Render ligand-network summaries with ΔΔG results and quality metrics.

## Installation

These scripts are designed to run within a standard OpenFE environment. No special installation is required beyond ensuring the necessary dependencies (OpenFE, RDKit, OpenFF) are available.

```bash
micromamba activate openfe_env
```

---

## Network Preparation

### Hybrid Topology (`prep_rbfe_hybridtop.py`)

Prepares a standard RBFE hybrid-topology setup from a receptor PDB and ligand SDF.

#### Usage Example
```bash
python3 prep_rbfe_hybridtop.py \
  --rec protein.pdb \
  --ligs ligands.sdf \
  --mapper kartograf \
  --network minimal_spanning \
  --windows 11 \
  --windowtime 5.0
```

#### Key Arguments
- `--mapper`: `kartograf` (default) or `lomap`.
- `--scorer`: `lomap` (default) or `kartograf_rmsd`.
- `--network`: `minimal_spanning` (default), `minimal_redundant`, `radial`, `maximal`, or `custom`.
- `--windows`: Number of lambda windows/replicas (default: 11).
- `--windowtime`: Production time per window in ns (default: 5.0).
- `--custom-network`: YAML file for `--network custom` mode.

### Separated Topology (`prep_rbfe_septop.py`)

Prepares a SepTop-based RBFE network (available since OpenFE 1.7.0). SepTop handles solvent and complex legs internally, so each transformation JSON represents a full RBFE cycle with `mapping=None`.

#### Usage Example
```bash
python3 prep_rbfe_septop.py \
  --rec protein.pdb \
  --ligs ligands.sdf \
  --network minimal_spanning \
  --protocol-repeats 1 \
  --host-min-distance 0.5 \
  --host-max-distance 1.5
```

#### Key Arguments
- All standard network args from the hybrid script.
- `--windows`: Number of lambda windows. If omitted, uses protocol defaults.
- `--windowtime`: Production time per window. If omitted, uses protocol defaults.
- `--equilibration-time`: Alchemical equilibration length (default: 2.0 ns).
- `--protocol-repeats`: Number of protocol repeats (default: 1).
- `--host-min-distance` / `--host-max-distance`: Restraint distance parameters (default: 0.5/1.5 nm).

---

## Analysis and Plotting

### `workup_hybridtop.py`

Processes results from the `results/` directory (organized in `repeatN` subfolders) and generates a comprehensive analysis.

```bash
python3 workup_hybridtop.py results/
```

It writes CSV summaries to `workup/` including:
- `summary_energies.csv`: Aggregated DDG results.
- `summary_convergence.csv`: Quality metrics (MBAR/HREX overlap, hysteresis).
- `summary_convergence_heatmap.png`: A visual summary of simulation quality.

### `plot_network.py`

Renders a publication-quality ligand network plot showing calculated ΔΔG values and color-coded edge quality metrics.

```bash
python3 plot_network.py workup/ --results-dir results/
```

- **Output**: `summary_ligand_network.png`
- **Visualization**: Shows aligned ligand structures as nodes and DDG values with error bars on edges. Edge colors (Green/Orange/Red) reflect the worst-performing overlap metric.

---

## Help

For a full list of options for any script:
```bash
python3 <script_name>.py --help
```
