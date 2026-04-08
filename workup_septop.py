#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Sequence


def _load_runtime_dependencies() -> None:
    global np
    global pd
    global unit
    global format_estimate_uncertainty
    global _collect_result_jsons
    global load_json
    global Measurement
    global FEMap

    if "np" in globals():
        return

    import numpy as np_module
    import pandas as pd_module
    from cinnabar import FEMap as FEMap_class
    from cinnabar import Measurement as Measurement_class
    from openfecli.commands.gather import _collect_result_jsons as collect_result_jsons
    from openfecli.commands.gather import format_estimate_uncertainty as format_uncertainty
    from openfecli.commands.gather import load_json as load_json_func
    from openff.units import unit as unit_module

    np = np_module
    pd = pd_module
    unit = unit_module
    format_estimate_uncertainty = format_uncertainty
    _collect_result_jsons = collect_result_jsons
    load_json = load_json_func
    Measurement = Measurement_class
    FEMap = FEMap_class


def existing_directory(path_value: str) -> Path:
    path = Path(path_value)
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"Directory does not exist: {path}")
    return path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze OpenFE SepTop result directories following the official "
            "SepTop analysis tutorial."
        ),
    )
    parser.add_argument(
        "results_dirs",
        nargs="+",
        type=existing_directory,
        help=(
            "One or more OpenFE SepTop result directories. You may pass repeat "
            "folders directly, for example results/repeat1 results/repeat2, or "
            "a parent directory such as results/ containing repeat subfolders."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory where ddg.tsv, dg.tsv, and ddg_raw.tsv will be written.",
    )
    return parser.parse_args(argv)


def _expand_results_directories(paths: Sequence[Path]) -> list[Path]:
    """
    Expand any parent results directory that contains repeat subfolders.

    This keeps the tutorial behaviour for explicit repeat directories while also
    supporting the common OpenFE layout of `results/repeatN`.
    """
    expanded: list[Path] = []
    seen: set[Path] = set()

    for path in paths:
        repeat_dirs = sorted(
            child for child in path.iterdir() if child.is_dir() and child.name.startswith("repeat")
        )

        if repeat_dirs:
            candidates = repeat_dirs
        else:
            candidates = [path]

        for candidate in candidates:
            resolved_candidate = candidate.resolve()
            if resolved_candidate in seen:
                continue
            seen.add(resolved_candidate)
            expanded.append(candidate)

    return expanded


def _load_valid_result_json(
    fpath: os.PathLike | str,
) -> tuple[tuple[str, str] | None, dict | None]:
    """Load the data from a results JSON into a dict."""
    result = load_json(fpath)
    try:
        names = _get_names(result)
    except (ValueError, IndexError):
        print(f"{fpath}: Missing ligand names. Skipping.")
        return None, None
    if result["estimate"] is None:
        errormsg = f"{fpath}: No 'estimate' found, assuming to be a failed simulation."
        raise ValueError(errormsg)
    return names, result


def _get_legs_from_result_jsons(
    result_fns: list[Path],
) -> dict[tuple[str, str], dict[str, list]]:
    """
    Iterate over a list of result JSONs and populate a dict of dicts with all data
    needed for results processing.
    """
    from collections import defaultdict

    ddgs = defaultdict(lambda: defaultdict(list))

    for result_fn in result_fns:
        names, result = _load_valid_result_json(result_fn)
        if names is None:
            continue

        ddgs[names]["overall"].append([result["estimate"], result["uncertainty"]])
        proto_key = [
            k for k in result["unit_results"].keys() if k.startswith("ProtocolUnitResult")
        ]
        for p in proto_key:
            outputs = result["unit_results"][p]["outputs"]
            if "unit_estimate" in outputs:
                simtype = outputs["simtype"]
                dg = outputs["unit_estimate"]
                dg_error = outputs["unit_estimate_error"]
                ddgs[names][simtype].append([dg, dg_error])
            elif "standard_state_correction_A" in outputs:
                corr_A = outputs["standard_state_correction_A"]
                corr_B = outputs["standard_state_correction_B"]
                ddgs[names]["standard_state_correction_A"].append(
                    [corr_A, 0 * unit.kilocalorie_per_mole]
                )
                ddgs[names]["standard_state_correction_B"].append(
                    [corr_B, 0 * unit.kilocalorie_per_mole]
                )

    return ddgs


def _get_names(result: dict) -> tuple[str, str]:
    """Get the ligand names from a unit's results data."""
    try:
        nm = list(result["unit_results"].values())[0]["name"]
    except KeyError as exc:
        raise ValueError("Failed to guess names") from exc

    toks = nm.split(",")
    toks = toks[1].split()
    return toks[1], toks[3]


def _error_std(result: dict[str, list]) -> float:
    """Calculate the error of the estimate as the std of the repeats."""
    return np.std([value[0].m for value in result["overall"]])


def _error_mbar(result: dict[str, list]) -> float:
    """
    Calculate the error of the estimate using the reported MBAR errors.
    This also takes into account repeated edges by using the average MBAR error.
    """
    complex_errors = [value[1].m for value in result["complex"]]
    solvent_errors = [value[1].m for value in result["solvent"]]
    return np.sqrt(np.mean(complex_errors) ** 2 + np.mean(solvent_errors) ** 2)


def extract_results_dict(
    results_files: list[os.PathLike | str],
) -> dict[tuple[str, str], dict[str, list]]:
    """
    Get a dictionary of SepTop results from a list of directories.
    """
    result_fns = _collect_result_jsons(results_files)
    sim_results = _get_legs_from_result_jsons(result_fns)
    return sim_results


def generate_ddg(results_dict: dict[tuple[str, str], dict[str, list]]) -> "pd.DataFrame":
    """Compute DDG values for the given results."""
    data = []
    repeats = {len(value["overall"]) for value in results_dict.values()}
    error_func = _error_mbar if 1 in repeats else _error_std

    for ligpair, results in sorted(results_dict.items()):
        ddg = np.mean([value[0].m for value in results["overall"]])
        error = error_func(results)
        mean_value, uncertainty = format_estimate_uncertainty(ddg, error, unc_prec=2)
        data.append((ligpair[0], ligpair[1], mean_value, uncertainty))

    return pd.DataFrame(
        data,
        columns=[
            "ligand_i",
            "ligand_j",
            "DDG(i->j) (kcal/mol)",
            "uncertainty (kcal/mol)",
        ],
    )


def generate_dg_mle(results_dict: dict[tuple[str, str], dict[str, list]]) -> "pd.DataFrame":
    """Compute MLE-derived DG values for the given results."""
    ddg_dataframe = generate_ddg(results_dict)
    fe_results = []

    for _, row in ddg_dataframe.iterrows():
        lig_a, lig_b, ddg_bind, bind_unc = row.tolist()
        fe_results.append(
            Measurement(
                labelA=lig_a,
                labelB=lig_b,
                DG=ddg_bind * unit.kilocalorie_per_mole,
                uncertainty=bind_unc * unit.kilocalorie_per_mole,
                computational=True,
            )
        )

    femap = FEMap()
    for entry in fe_results:
        femap.add_measurement(entry)

    femap.generate_absolute_values()
    df = femap.get_absolute_dataframe()
    df = df.iloc[:, :3]
    df.rename({"label": "ligand"}, axis="columns", inplace=True)
    return df


def generate_dg_raw(results_dict: dict[tuple[str, str], dict[str, list]]) -> "pd.DataFrame":
    """Get the raw DG values for every leg in the SepTop transformation cycles."""
    data = []

    for ligpair, results in sorted(results_dict.items()):
        for simtype, repeats in sorted(results.items()):
            if simtype == "overall":
                continue

            for repeat in repeats:
                mean_value, uncertainty = format_estimate_uncertainty(
                    repeat[0].m,
                    repeat[1].m,
                    unc_prec=2,
                )
                data.append((simtype, ligpair[0], ligpair[1], mean_value, uncertainty))

    return pd.DataFrame(
        data,
        columns=[
            "leg",
            "ligand_i",
            "ligand_j",
            "DG(i->j) (kcal/mol)",
            "uncertainty (kcal/mol)",
        ],
    )


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    _load_runtime_dependencies()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    input_dirs = _expand_results_directories([Path(path) for path in args.results_dirs])
    results_dict = extract_results_dict(input_dirs)

    df_ddg = generate_ddg(results_dict)
    df_dg = generate_dg_mle(results_dict)
    df_raw = generate_dg_raw(results_dict)

    ddg_path = output_dir / "ddg.tsv"
    dg_path = output_dir / "dg.tsv"
    raw_path = output_dir / "ddg_raw.tsv"

    df_ddg.to_csv(ddg_path, sep="\t", lineterminator="\n", index=False)
    df_dg.to_csv(dg_path, sep="\t", lineterminator="\n", index=False)
    df_raw.to_csv(raw_path, sep="\t", lineterminator="\n", index=False)

    print(f"Wrote {ddg_path}")
    print(f"Wrote {dg_path}")
    print(f"Wrote {raw_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
