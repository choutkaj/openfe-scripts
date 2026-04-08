# openfe-scripts

A collection of helper scripts that extend the [OpenFE](https://github.com/OpenFreeEnergy/openfe) CLI workflow with extra utilities for setting up, running, and analyzing Relative Free Energy (RFE) calculations.

## Current Scripts

- **`prep_rbfe_hybridtop.py`**: Prepare hybrid-topology RBFE networks and write transformation JSON files.
- **`prep_rbfe_septop.py`**: Prepare separated-topology (SepTop) RBFE networks and write transformation JSON files.
- **`workup_hybridtop.py`**: Summarize and analyze completed OpenFE calculations, generating CSV reports and convergence plots.
- **`workup_septop.py`**: Analyze SepTop result directories following the official OpenFE SepTop analysis tutorial and write TSV summaries.
- **`plot_network.py`**: Render ligand-network summaries with ΔΔG results and quality metrics.

## Installation

These scripts should run within a standard OpenFE environment. Just activate your OpenFE environment and run the scripts. For example:

```bash
micromamba activate openfe
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

Prepares a SepTop-based RBFE network. Experimental.

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
- `summary_energies.csv`: Aggregated ΔΔG results.
- `summary_convergence.csv`: Quality metrics (MBAR/HREX overlap, hysteresis).
- `summary_convergence_heatmap.png`: A visual summary of simulation quality.

Here is an example of the output heatmap with quality and convergence metrics. Still a lot of red, yayks! Gotta sample more.

![Heatmap with metrics](https://github.com/choutkaj/openfe-scripts/blob/main/summary_convergence_heatmap.png)


### `workup_septop.py`

Processes one or more SepTop result directories exactly along the lines of the
official OpenFE SepTop analysis tutorial and writes:
- `ddg.tsv`: Edgewise relative binding free energies.
- `dg.tsv`: MLE-derived absolute binding free energies.
- `ddg_raw.tsv`: Per-leg raw SepTop free energies.

If your OpenFE outputs are organized as `results/repeat1`, `results/repeat2`, ...
you can point the script at the parent `results/` directory and it will analyze
all repeat subfolders automatically.

```bash
python3 workup_septop.py results/ --output-dir workup_septop/
```

### `plot_network.py`

Renders a ligand network plot showing calculated ΔΔG values and color-coded edge quality metrics. The network plotting is still a bit rough around the edges.

```bash
python3 plot_network.py workup/ --results-dir results/
```


## Help

For a full list of options for any script:
```bash
python3 <script_name>.py --help
```
