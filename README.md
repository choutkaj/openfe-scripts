# openfe-scripts

A collection of helper scripts that extend the [OpenFE](https://github.com/OpenFreeEnergy/openfe) CLI workflow with extra utilities for setting up, running, and analyzing OpenFE calculations.

## Current Scripts

- **`prep-rbfe-hybridtop.py`**: Prepare hybrid-topology RBFE networks and write transformation JSON files.
- **`prep-rbfe-septop.py`**: Prepare separated-topology (SepTop) RBFE networks and write transformation JSON files.
- **`workup-hybridtop.py`**: Summarize and analyze completed OpenFE calculations, generating CSV reports and convergence plots.
- **`workup-septop.py`**: Analyze SepTop result directories following the official OpenFE SepTop analysis tutorial and write TSV summaries.
- **`plot-network.py`**: Render ligand-network summaries with ΔΔG results and quality metrics.

## Installation

These scripts are meant to be used from within a standard OpenFE environment.
The recommended setup is to clone this repository and install it in
editable mode inside your OpenFE environment. 

```bash
micromamba activate openfe    # modify if you use something else than micromamba
git clone https://github.com/choutkaj/openfe-scripts.git
cd openfe-scripts
python -m pip install --no-build-isolation -e .
```

The `--no-build-isolation` flag helps in a vanilla OpenFE environment by
preventing `pip` from trying to download build-time packaging tools.

After installation, the commands can be run directly:

```bash
prep-rbfe-septop.py --help
prep-rbfe-hybridtop.py --help
workup-septop.py --help
workup-hybridtop.py --help
plot-network.py --help
```

When you update the repository later with `git pull`, the installed commands
will automatically use the updated code because the install is editable.

If you prefer not to install the repository, you can still run the files
directly from the checkout with `python3 /path/to/openfe-scripts/<script>.py`.

The examples below assume the editable install, so the commands are shown
without `python3`.

---

## Network Preparation

### Hybrid Topology (`prep-rbfe-hybridtop.py`)

Prepares a standard RBFE hybrid-topology setup from a receptor PDB and SDF with ligands.

#### Minimal Example
```bash
prep-rbfe-hybridtop.py \
  --rec protein.pdb \
  --ligs ligands.sdf \
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

Prepares a SepTop RBFE network.

#### Minimal Example
```bash
prep-rbfe-septop.py \
  --rec protein.pdb \
  --ligs ligands.sdf \
  --solvent-windows 27 \
  --complex-windows 19 \
  --solvent-windowtime 10.0 \
  --complex-windowtime 10.0 \
  --smart-lambda
```

#### Arguments and Defaults
- Required: `--rec` receptor PDB path, `--ligs` ligand SDF path.
- `--mapper`: `kartograf` by default. Choices: `kartograf`, `lomap`.
- `--scorer`: `lomap` by default. Choices: `lomap`, `kartograf_rmsd`.
- `--network`: `minimal_spanning` by default. Choices: `minimal_spanning`, `minimal_redundant`, `radial`, `maximal`, `custom`.
- `--custom-network`: no default; only valid and required when `--network custom` is selected.
- `--central-ligand`: no default; if omitted for `radial`, the first ligand in the SDF is used.
- `--solvent-windows`: no explicit script default; if omitted, the solvent leg keeps the current OpenFE SepTop default lambda schedule with `27` windows.
- `--complex-windows`: no explicit script default; if omitted, the complex leg keeps the current OpenFE SepTop default lambda schedule with `19` windows.
- `--smart-lambda`: off by default; when used together with `--solvent-windows` and/or `--complex-windows`, extra windows are concentrated in regions where the combined lambda schedule changes most instead of stretching the defaults by window index.
- `--solvent-windowtime`: no explicit script default; if omitted, the solvent leg keeps the OpenFE SepTop default multistate production length of `10.0` ns.
- `--complex-windowtime`: no explicit script default; if omitted, the complex leg keeps the OpenFE SepTop default multistate production length of `10.0` ns.
- Note: the separate OpenFE endstate-equilibration settings are currently `equilibration_length_nvt = 0.1` ns, `equilibration_length = 0.1` ns, and `production_length = 2.0` ns for both solvent and complex. These are different settings from `--solvent-windowtime` / `--complex-windowtime`.
- `--small-molecule-forcefield`: no explicit script default; if omitted, the current OpenFE default is `openff-2.2.1`. Confirmed working values so far: `openff-2.2.1`, `gaff-2.2.20`.
- `--partial-charge-method`: `am1bcc` by default. Choices: `am1bcc`, `nagl`. If `nagl` is selected, the OpenFE/OpenFF default NAGL model is used.
- `--equilibration-time`: no explicit script default; if omitted, the current OpenFE multistate equilibration length is `1.0` ns for both solvent and complex.
- `--protocol-repeats`: `1` by default in this script. OpenFE SepTop itself currently defaults to `3`, but the wrapper intentionally overrides that to `1`.
- `--host-min-distance`: no explicit script default; if omitted, the current OpenFE SepTop default is `0.5` nm.
- `--host-max-distance`: no explicit script default; if omitted, the current OpenFE SepTop default is `1.5` nm.
- `--output-dir`: current directory (`.`) by default.

---

## Analysis and Plotting

### `workup-hybridtop.py`

Processes results from the `results/` directory (organized in `repeatN` subfolders) and generates a comprehensive analysis.

```bash
workup-hybridtop.py results/
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
workup-septop.py results/ --output-dir workup_septop/
```

#### Arguments and Defaults
- `results_dirs`: one or more required positional directories.
- `--output-dir`: current directory (`.`) by default.

### `plot-network.py`

Renders a ligand network plot showing calculated ΔΔG values and color-coded edge quality metrics. The network plotting is still a bit rough around the edges.

```bash
plot-network.py workup/ --results-dir results/
```

#### Arguments and Defaults
- `workup_dir`: optional positional argument; defaults to `example/workup`.
- `--results-dir`: no explicit default; if omitted, the script uses `<workup_dir>/../results` when that directory exists, otherwise ligand depictions are skipped.
- `--output`: no explicit default; if omitted, output is written to `<workup_dir>/summary_ligand_network.png`.


## Help

For a full list of options for any script:
```bash
<script_name>.py --help
```
