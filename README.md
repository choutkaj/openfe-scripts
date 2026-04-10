# openfe-scripts

A collection of helper scripts that extend the [OpenFE](https://github.com/OpenFreeEnergy/openfe) CLI workflow with extra utilities for setting up, running, and analyzing Relative Free Energy (RFE) calculations.

## Current Scripts

- **`prep-rbfe-hybridtop.py`**: Prepare hybrid-topology RBFE networks and write transformation JSON files.
- **`prep-rbfe-septop.py`**: Prepare separated-topology (SepTop) RBFE networks and write transformation JSON files.
- **`workup-hybridtop.py`**: Summarize and analyze completed OpenFE calculations, generating CSV reports and convergence plots.
- **`workup-septop.py`**: Analyze SepTop result directories following the official OpenFE SepTop analysis tutorial and write TSV summaries.
- **`plot-network.py`**: Render ligand-network summaries with ΔΔG results and quality metrics.

## Installation

These scripts should run within a standard OpenFE environment. Just activate your OpenFE environment and run the scripts. For example:

```bash
micromamba activate openfe
```

---

## Network Preparation

### Hybrid Topology (`prep-rbfe-hybridtop.py`)

Prepares a standard RBFE hybrid-topology setup from a receptor PDB and ligand SDF.

#### Usage Example
```bash
python3 prep-rbfe-hybridtop.py \
  --rec protein.pdb \
  --ligs ligands.sdf \
  --partial-charge-method am1bcc \
  --mapper kartograf \
  --network minimal_spanning \
  --small-molecule-forcefield openff-2.2.1 \
  --windows 11 \
  --windowtime 5.0
```

#### Arguments and Defaults
- Required: `--rec` receptor PDB path, `--ligs` ligand SDF path.
- `--mapper`: `kartograf` by default. Choices: `kartograf`, `lomap`.
- `--scorer`: `lomap` by default. Choices: `lomap`, `kartograf_rmsd`.
- `--network`: `minimal_spanning` by default. Choices: `minimal_spanning`, `minimal_redundant`, `radial`, `maximal`, `custom`.
- `--custom-network`: no default; only valid and required when `--network custom` is selected.
- `--central-ligand`: no default; if omitted for `radial`, the first ligand in the SDF is used.
- `--windows`: `11` by default.
- `--windowtime`: `5.0` ns by default.
- `--small-molecule-forcefield`: `openff-2.2.1` by default. Confirmed working values: `openff-2.2.1`, `gaff-2.2.20`.
- `--partial-charge-method`: `am1bcc` by default. Choices: `am1bcc`, `nagl`. If `nagl` is selected, the OpenFE/OpenFF default NAGL model is used.
- `--output-dir`: current directory (`.`) by default.

### Separated Topology (`prep-rbfe-septop.py`)

Prepares a SepTop-based RBFE network. Experimental.

#### Usage Example
```bash
python3 prep-rbfe-septop.py \
  --rec protein.pdb \
  --ligs ligands.sdf \
  --partial-charge-method am1bcc \
  --network minimal_spanning \
  --small-molecule-forcefield openff-2.2.1 \
  --protocol-repeats 1 \
  --host-min-distance 0.5 \
  --host-max-distance 1.5
```

#### Arguments and Defaults
- Required: `--rec` receptor PDB path, `--ligs` ligand SDF path.
- `--mapper`: `kartograf` by default. Choices: `kartograf`, `lomap`.
- `--scorer`: `lomap` by default. Choices: `lomap`, `kartograf_rmsd`.
- `--network`: `minimal_spanning` by default. Choices: `minimal_spanning`, `minimal_redundant`, `radial`, `maximal`, `custom`.
- `--custom-network`: no default; only valid and required when `--network custom` is selected.
- `--central-ligand`: no default; if omitted for `radial`, the first ligand in the SDF is used.
- `--windows`: no default; if omitted, the SepTop protocol default lambda schedule is used.
- `--windowtime`: no default; if omitted, the SepTop protocol default production length is used.
- `--small-molecule-forcefield`: `openff-2.2.1` by default. Confirmed working values so far: `openff-2.2.1`, `gaff-2.2.20`.
- `--partial-charge-method`: `am1bcc` by default. Choices: `am1bcc`, `nagl`. If `nagl` is selected, the OpenFE/OpenFF default NAGL model is used.
- `--equilibration-time`: `2.0` ns by default.
- `--protocol-repeats`: `1` by default.
- `--host-min-distance`: `0.5` nm by default.
- `--host-max-distance`: `1.5` nm by default.
- `--output-dir`: current directory (`.`) by default.

---

## Analysis and Plotting

### `workup-hybridtop.py`

Processes results from the `results/` directory (organized in `repeatN` subfolders) and generates a comprehensive analysis.

```bash
python3 workup-hybridtop.py results/
```

It writes CSV summaries to `workup/` including:
- `summary_energies.csv`: Aggregated ΔΔG results.
- `summary_convergence.csv`: Quality metrics (MBAR/HREX overlap, hysteresis).
- `summary_convergence_heatmap.png`: A visual summary of simulation quality.

Here is an example of the output heatmap with quality and convergence metrics. Still a lot of red, yayks! Gotta sample more.

![Heatmap with metrics](https://github.com/choutkaj/openfe-scripts/blob/main/summary_convergence_heatmap.png)

#### Arguments and Defaults
- `results_dir`: optional positional argument; defaults to `results`.
- `--output-dir`: no explicit default; if omitted, output is written to `<results_dir>/../workup`.
- `--no-export`: off by default.
- `--no-plots`: off by default.


### `workup-septop.py`

Processes one or more SepTop result directories exactly along the lines of the
official OpenFE SepTop analysis tutorial and writes:
- `ddg.tsv`: Edgewise relative binding free energies.
- `dg.tsv`: MLE-derived absolute binding free energies.
- `ddg_raw.tsv`: Per-leg raw SepTop free energies.

If your OpenFE outputs are organized as `results/repeat1`, `results/repeat2`, ...
you can point the script at the parent `results/` directory and it will analyze
all repeat subfolders automatically.

```bash
python3 workup-septop.py results/ --output-dir workup_septop/
```

#### Arguments and Defaults
- `results_dirs`: one or more required positional directories.
- `--output-dir`: current directory (`.`) by default.

### `plot-network.py`

Renders a ligand network plot showing calculated ΔΔG values and color-coded edge quality metrics. The network plotting is still a bit rough around the edges.

```bash
python3 plot-network.py workup/ --results-dir results/
```

#### Arguments and Defaults
- `workup_dir`: optional positional argument; defaults to `example/workup`.
- `--results-dir`: no explicit default; if omitted, the script uses `<workup_dir>/../results` when that directory exists, otherwise ligand depictions are skipped.
- `--output`: no explicit default; if omitted, output is written to `<workup_dir>/summary_ligand_network.png`.


## Help

For a full list of options for any script:
```bash
python3 <script_name>.py --help
```
