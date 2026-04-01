from __future__ import annotations

import argparse
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import numpy as np
import pandas as pd
import yaml
from plot_network import plot_summary_ligand_network


ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"
TIME_UNIT_TO_NS = {
    "femtosecond": 1e-6,
    "picosecond": 1e-3,
    "nanosecond": 1.0,
}
TARGET_HYSTERESIS_FRACTIONS = (0.5, 0.7, 0.9)

KBT_TO_KCAL_MOL_AT_298K = 0.592
QUALITY_RULES = {
    "general": None,
    "repeats": {
        "direction": "higher_is_better",
        "good": 3,
        "ok": 2,
        "description": "bad < 2, ok = 2, good ≥ 3",
        "good_label": "≥ 3",
        "ok_label": "= 2",
        "bad_label": "< 2",
    },
    "lambda_windows": {
        "direction": "higher_is_better",
        "good": 11,
        "ok": 11,
        "description": "bad < 11, good ≥ 11",
        "good_label": "≥ 11",
        "ok_label": "",
        "bad_label": "< 11",
    },
    "simulation_time_per_window": {
        "direction": "higher_is_better",
        "good": 11,
        "ok": 11,
        "description": "bad < 11, good ≥ 11",
        "good_label": "≥ 11",
        "ok_label": "",
        "bad_label": "< 11",
    },
    "perturbed_heavy_atoms": {
        "direction": "lower_is_better",
        "good": 8,
        "ok": 14,
        "description": "good ≤ 8, ok 8 < x ≤ 14, bad > 14",
        "good_label": "≤ 8",
        "ok_label": "(8, 14]",
        "bad_label": "> 14",
    },
    "mbar_overlap": {
        "direction": "higher_is_better",
        "good": 0.10,
        "ok": 0.03,
        "description": "bad < 0.03, ok 0.03 ≤ x < 0.10, good ≥ 0.10",
        "good_label": "≥ 0.10",
        "ok_label": "[0.03, 0.10)",
        "bad_label": "< 0.03",
    },
    "hrex_overlap": {
        "direction": "higher_is_better",
        "good": 0.20,
        "ok": 0.10,
        "description": "bad < 0.10, ok 0.10 ≤ x < 0.20, good ≥ 0.20",
        "good_label": "≥ 0.20",
        "ok_label": "[0.10, 0.20)",
        "bad_label": "< 0.10",
    },
    "hysteresis": {
        "direction": "lower_is_better",
        "good": 1.0 * KBT_TO_KCAL_MOL_AT_298K,
        "ok": 2.0 * KBT_TO_KCAL_MOL_AT_298K,
        "description": "good ≤ 1 kBT, ok 1 < x ≤ 2 kBT, bad > 2 kBT",
        "good_label": "≤ 1 kBT",
        "ok_label": "(1, 2] kBT",
        "bad_label": "> 2 kBT",
    },
}
QUALITY_REFERENCE_RULES = [
    ("Repeats", "repeats"),
    ("Lambda windows", "lambda_windows"),
    ("Sim. time / window", "simulation_time_per_window"),
    ("Perturbed h. atoms", "perturbed_heavy_atoms"),
    ("MBAR overlap", "mbar_overlap"),
    ("HREX overlap", "hrex_overlap"),
    ("Hysteresis", "hysteresis"),
]

QUALITY_METRIC_SPECS = {
    "n_repeats": {
        "label": "Repeats",
        "format": "{:.0f}",
        "rule": "repeats",
    },
    "perturbed_heavy_atoms_total": {
        "label": "Perturbed heavy atoms",
        "format": "{:.0f}",
        "rule": "perturbed_heavy_atoms",
    },
    "n_lambda_windows": {
        "label": "λ windows",
        "format": "{:.0f}",
        "rule": "lambda_windows",
    },
    "simulation_time_ns_per_window": {
        "label": "Sim. time / window (ns)",
        "format": "{:.1f}",
        "rule": "simulation_time_per_window",
    },
    "complex_mbar_offdiag_min": {
        "label": "Complex MBAR min",
        "format": "{:.2f}",
        "rule": "mbar_overlap",
    },
    "solvent_mbar_offdiag_min": {
        "label": "Solvent MBAR min",
        "format": "{:.2f}",
        "rule": "mbar_overlap",
    },
    "complex_replica_exchange_offdiag_min": {
        "label": "Complex HREX min",
        "format": "{:.2f}",
        "rule": "hrex_overlap",
    },
    "solvent_replica_exchange_offdiag_min": {
        "label": "Solvent HREX min",
        "format": "{:.2f}",
        "rule": "hrex_overlap",
    },
    "complex_hysteresis_50pct_kcal_mol": {
        "label": "Complex hysteresis @ 50%",
        "format": "{:.1f}",
        "rule": "hysteresis",
    },
    "solvent_hysteresis_50pct_kcal_mol": {
        "label": "Solvent hysteresis @ 50%",
        "format": "{:.1f}",
        "rule": "hysteresis",
    },
    "complex_hysteresis_70pct_kcal_mol": {
        "label": "Complex hysteresis @ 70%",
        "format": "{:.1f}",
        "rule": "hysteresis",
    },
    "solvent_hysteresis_70pct_kcal_mol": {
        "label": "Solvent hysteresis @ 70%",
        "format": "{:.1f}",
        "rule": "hysteresis",
    },
    "complex_hysteresis_90pct_kcal_mol": {
        "label": "Complex hysteresis @ 90%",
        "format": "{:.1f}",
        "rule": "hysteresis",
    },
    "solvent_hysteresis_90pct_kcal_mol": {
        "label": "Solvent hysteresis @ 90%",
        "format": "{:.1f}",
        "rule": "hysteresis",
    },
    "min_leg_mbar_offdiag_min": {
        "label": "Worst MBAR min",
        "format": "{:.2f}",
        "rule": "mbar_overlap",
    },
    "min_leg_replica_exchange_offdiag_min": {
        "label": "Worst HREX min",
        "format": "{:.2f}",
        "rule": "hrex_overlap",
    },
    "max_leg_hysteresis_50pct_kcal_mol": {
        "label": "Worst hysteresis @ 50%",
        "format": "{:.1f}",
        "rule": "hysteresis",
    },
}
QUALITY_COLORS = {
    "good": "#b7e4c7",
    "ok": "#ffd6a5",
    "bad": "#ffadad",
    "missing": "#ffffff",
}
QUALITY_COLUMN_GROUPS = [
    (
        "General",
        [
            "n_repeats",
            "n_lambda_windows",
            "simulation_time_ns_per_window",
            "perturbed_heavy_atoms_total",
        ],
    ),
    (
        "MBAR overlap",
        [
            "complex_mbar_offdiag_min",
            "solvent_mbar_offdiag_min",
            "min_leg_mbar_offdiag_min",
        ],
    ),
    (
        "HREX overlap",
        [
            "complex_replica_exchange_offdiag_min",
            "solvent_replica_exchange_offdiag_min",
            "min_leg_replica_exchange_offdiag_min",
        ],
    ),
    (
        "Hysteresis [kcal/mol]",
        [
            "complex_hysteresis_50pct_kcal_mol",
            "complex_hysteresis_70pct_kcal_mol",
            "complex_hysteresis_90pct_kcal_mol",
            "solvent_hysteresis_50pct_kcal_mol",
            "solvent_hysteresis_70pct_kcal_mol",
            "solvent_hysteresis_90pct_kcal_mol",
            "max_leg_hysteresis_50pct_kcal_mol",
        ],
    ),
]
EMPHASIZED_QUALITY_COLUMNS = {
    "min_leg_mbar_offdiag_min",
    "min_leg_replica_exchange_offdiag_min",
    "max_leg_hysteresis_50pct_kcal_mol",
}
QUALITY_STATUS_ORDER = {
    "good": 0,
    "ok": 1,
    "bad": 2,
}
LEG_RESULTS_EXPORT_COLUMNS = (
    "ligand_i",
    "ligand_j",
    "leg",
    "repeat",
    "repeat_index",
    "result_json",
    "dg_kcal_mol",
    "dg_error_kcal_mol",
    "n_lambda_windows",
    "simulation_time_ns_per_window",
    "equilibration_time_ns_per_window",
    "total_time_ns_per_window",
    "planned_total_sampling_ns",
    "mbar_overlap_scalar",
    "mbar_second_eigenvalue",
    "replica_exchange_second_eigenvalue",
    "replica_exchange_spectral_gap",
    "hysteresis_profile_fractions",
    "hysteresis_profile_kcal_mol",
    "throughput_ns_per_day",
    "number_of_uncorrelated_samples",
    "statistical_inefficiency",
    "hysteresis_50pct_kcal_mol",
    "hysteresis_50pct_fraction_actual",
    "hysteresis_70pct_kcal_mol",
    "hysteresis_70pct_fraction_actual",
    "hysteresis_90pct_kcal_mol",
    "hysteresis_90pct_fraction_actual",
    "mbar_offdiag_values",
    "mbar_offdiag_min",
    "mbar_offdiag_max",
    "mbar_offdiag_mean",
    "mbar_offdiag_upper_values",
    "mbar_offdiag_lower_values",
    "replica_exchange_offdiag_values",
    "replica_exchange_offdiag_min",
    "replica_exchange_offdiag_max",
    "replica_exchange_offdiag_mean",
    "replica_exchange_offdiag_upper_values",
    "replica_exchange_offdiag_lower_values",
    "ligand_i_heavy_atoms",
    "ligand_j_heavy_atoms",
    "mapped_heavy_atoms_core",
    "perturbed_heavy_atoms_ligand_i",
    "perturbed_heavy_atoms_ligand_j",
    "perturbed_heavy_atoms_total",
    "protein_RMSD_mean",
    "protein_RMSD_max",
    "ligand_RMSD_mean",
    "ligand_RMSD_max",
    "ligand_COM_drift_mean",
    "ligand_COM_drift_max",
    "protein_2D_RMSD_mean",
    "protein_2D_RMSD_max",
)
REPEAT_RESULTS_SHARED_EXPORT_COLUMNS = (
    "ligand_i",
    "ligand_j",
    "repeat",
    "repeat_index",
    "n_legs_found",
    "legs_found",
    "ligand_i_heavy_atoms",
    "ligand_j_heavy_atoms",
    "mapped_heavy_atoms_core",
    "perturbed_heavy_atoms_ligand_i",
    "perturbed_heavy_atoms_ligand_j",
    "perturbed_heavy_atoms_total",
    "ddg_kcal_mol",
    "ddg_error_kcal_mol",
    "min_leg_mbar_offdiag_min",
    "min_leg_replica_exchange_offdiag_min",
    "max_leg_hysteresis_50pct_kcal_mol",
)
REPEAT_RESULTS_LEG_EXPORT_FIELDS = (
    "result_json",
    "dg_kcal_mol",
    "dg_error_kcal_mol",
    "n_lambda_windows",
    "simulation_time_ns_per_window",
    "equilibration_time_ns_per_window",
    "total_time_ns_per_window",
    "planned_total_sampling_ns",
    "mbar_overlap_scalar",
    "mbar_second_eigenvalue",
    "replica_exchange_second_eigenvalue",
    "replica_exchange_spectral_gap",
    "hysteresis_profile_fractions",
    "hysteresis_profile_kcal_mol",
    "throughput_ns_per_day",
    "number_of_uncorrelated_samples",
    "statistical_inefficiency",
    "hysteresis_50pct_kcal_mol",
    "hysteresis_50pct_fraction_actual",
    "hysteresis_70pct_kcal_mol",
    "hysteresis_70pct_fraction_actual",
    "hysteresis_90pct_kcal_mol",
    "hysteresis_90pct_fraction_actual",
    "mbar_offdiag_values",
    "mbar_offdiag_min",
    "mbar_offdiag_max",
    "mbar_offdiag_mean",
    "mbar_offdiag_upper_values",
    "mbar_offdiag_lower_values",
    "replica_exchange_offdiag_values",
    "replica_exchange_offdiag_min",
    "replica_exchange_offdiag_max",
    "replica_exchange_offdiag_mean",
    "replica_exchange_offdiag_upper_values",
    "replica_exchange_offdiag_lower_values",
    "ligand_i_heavy_atoms",
    "ligand_j_heavy_atoms",
    "mapped_heavy_atoms_core",
    "perturbed_heavy_atoms_ligand_i",
    "perturbed_heavy_atoms_ligand_j",
    "perturbed_heavy_atoms_total",
    "protein_RMSD_mean",
    "protein_RMSD_max",
    "ligand_RMSD_mean",
    "ligand_RMSD_max",
    "ligand_COM_drift_mean",
    "ligand_COM_drift_max",
    "protein_2D_RMSD_mean",
    "protein_2D_RMSD_max",
)
LEG_SUMMARY_EXPORT_COLUMNS = (
    "ligand_i",
    "ligand_j",
    "leg",
    "n_repeats",
    "dg_mean_kcal_mol",
    "dg_std_kcal_mol",
    "dg_sem_kcal_mol",
    "dg_error_mean_kcal_mol",
    "mbar_overlap_scalar_mean",
    "mbar_offdiag_min_mean",
    "mbar_offdiag_max_mean",
    "replica_exchange_offdiag_min_mean",
    "replica_exchange_offdiag_max_mean",
    "hysteresis_50pct_mean_kcal_mol",
    "hysteresis_70pct_mean_kcal_mol",
    "hysteresis_90pct_mean_kcal_mol",
    "throughput_ns_per_day_mean",
)
SUMMARY_ENERGIES_EXPORT_COLUMNS = (
    "ligand_i",
    "ligand_j",
    "ddg_mean_kcal_mol",
    "ddg_std_kcal_mol",
    "ddg_sem_kcal_mol",
    "ddg_error_mean_kcal_mol",
    "complex_dg_mean_kcal_mol",
    "solvent_dg_mean_kcal_mol",
)
SUMMARY_CONVERGENCE_EXPORT_COLUMNS = (
    "ligand_i",
    "ligand_j",
    "n_repeats",
    "n_lambda_windows",
    "simulation_time_ns_per_window",
    "perturbed_heavy_atoms_total",
    "complex_mbar_offdiag_min",
    "solvent_mbar_offdiag_min",
    "complex_replica_exchange_offdiag_min",
    "solvent_replica_exchange_offdiag_min",
    "complex_hysteresis_50pct_kcal_mol",
    "complex_hysteresis_70pct_kcal_mol",
    "complex_hysteresis_90pct_kcal_mol",
    "solvent_hysteresis_50pct_kcal_mol",
    "solvent_hysteresis_70pct_kcal_mol",
    "solvent_hysteresis_90pct_kcal_mol",
    "min_leg_mbar_offdiag_min",
    "min_leg_replica_exchange_offdiag_min",
    "max_leg_hysteresis_50pct_kcal_mol",
)
REPEAT_RESULTS_EXPORT_COLUMNS = (
    *REPEAT_RESULTS_SHARED_EXPORT_COLUMNS,
    *tuple(
        f"{leg_name}_{field}"
        for leg_name in ("complex", "solvent")
        for field in REPEAT_RESULTS_LEG_EXPORT_FIELDS
    ),
)
CSV_EXPORT_COLUMNS = {
    "leg_results.csv": LEG_RESULTS_EXPORT_COLUMNS,
    "repeat_results.csv": REPEAT_RESULTS_EXPORT_COLUMNS,
    "leg_summary.csv": LEG_SUMMARY_EXPORT_COLUMNS,
    "summary_energies.csv": SUMMARY_ENERGIES_EXPORT_COLUMNS,
    "summary_convergence.csv": SUMMARY_CONVERGENCE_EXPORT_COLUMNS,
}


@dataclass
class WorkupResults:
    leg_results: pd.DataFrame
    repeat_results: pd.DataFrame
    leg_summary: pd.DataFrame
    transformation_summary: pd.DataFrame
    summary_energies: pd.DataFrame
    summary_convergence: pd.DataFrame
    results_root: Path | None = None

    def export_csv(self, output_dir: str | Path) -> dict[str, Path]:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        exports = {
            "leg_results.csv": self.leg_results,
            "repeat_results.csv": self.repeat_results,
            "leg_summary.csv": self.leg_summary,
            "summary_energies.csv": self.summary_energies,
            "summary_convergence.csv": self.summary_convergence,
        }

        written_paths: dict[str, Path] = {}
        for filename, frame in exports.items():
            destination = output_path / filename
            export_columns = [column for column in CSV_EXPORT_COLUMNS[filename] if column in frame.columns]
            export_frame = frame.loc[:, export_columns].copy()
            export_frame.to_csv(destination, index=False)
            written_paths[filename] = destination

        return written_paths

    def export_plots(self, output_dir: str | Path) -> dict[str, Path]:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        written_paths: dict[str, Path] = {}
        convergence_plot = output_path / "summary_convergence_heatmap.png"
        plot_summary_convergence_heatmap(self.summary_convergence, convergence_plot)
        written_paths["summary_convergence_heatmap.png"] = convergence_plot

        network_plot = output_path / "summary_ligand_network.png"
        plot_summary_ligand_network(
            self.transformation_summary,
            network_plot,
            results_root=self.results_root,
        )
        if network_plot.exists():
            written_paths["summary_ligand_network.png"] = network_plot

        return written_paths

    def pivot_metric(
        self,
        metric: str = "ddg_mean_kcal_mol",
        frame: str = "transformation_summary",
        index: str = "ligand_i",
        columns: str = "ligand_j",
    ) -> pd.DataFrame:
        source = getattr(self, frame)
        return source.pivot(index=index, columns=columns, values=metric)


def load_results(results_dir: str | Path) -> WorkupResults:
    results_path = Path(results_dir)
    result_files = sorted(results_path.glob("repeat*/rbfe_*.json"))

    if not result_files:
        raise FileNotFoundError(f"No result JSON files found under {results_path}")

    leg_records = [parse_leg_result(path, results_path) for path in result_files]
    leg_results = _records_to_frame(
        leg_records,
        sort_columns=["repeat_index", "ligand_i", "ligand_j", "leg"],
    )
    repeat_results = build_repeat_results(leg_results)
    leg_summary = build_leg_summary(leg_results)
    transformation_summary = build_transformation_summary(repeat_results)
    summary_energies = build_summary_energies(transformation_summary)
    summary_convergence = build_summary_convergence(transformation_summary)

    return WorkupResults(
        leg_results=leg_results,
        repeat_results=repeat_results,
        leg_summary=leg_summary,
        transformation_summary=transformation_summary,
        summary_energies=summary_energies,
        summary_convergence=summary_convergence,
        results_root=results_path.resolve(),
    )


def parse_leg_result(result_path: Path, results_root: Path) -> dict[str, Any]:
    payload = json.loads(result_path.read_text())
    analysis_result = _find_analysis_result(payload)
    setup_result = analysis_result["inputs"]["setup_results"]
    simulation_result = analysis_result["inputs"]["simulation_results"]
    protocol_settings = analysis_result["inputs"]["protocol"]["settings"]
    simulation_settings = protocol_settings["simulation_settings"]
    lambda_settings = protocol_settings["lambda_settings"]

    state_a_name = setup_result["inputs"]["stateA"]["name"]
    state_b_name = setup_result["inputs"]["stateB"]["name"]
    ligand_i, leg_i = _split_state_name(state_a_name)
    ligand_j, leg_j = _split_state_name(state_b_name)

    if leg_i != leg_j:
        raise ValueError(
            f"Mismatched legs in {result_path.name}: {state_a_name} vs {state_b_name}"
        )

    transformation_id = f"{ligand_i}__{ligand_j}"
    repeat_name, repeat_index = _parse_repeat_folder(result_path.parent.name)

    unit_mbar = analysis_result["outputs"]["unit_mbar_overlap"]
    mbar_matrix = decode_ndarray(unit_mbar["matrix"])
    mbar_eigenvalues = decode_ndarray(unit_mbar["eigenvalues"])

    repex_stats = analysis_result["outputs"]["replica_exchange_statistics"]
    repex_matrix = decode_ndarray(repex_stats["matrix"])
    repex_eigenvalues = decode_ndarray(repex_stats["eigenvalues"])

    convergence = analysis_result["outputs"]["forward_and_reverse_energies"]
    fractions = decode_ndarray(convergence["fractions"])
    forward_dgs = decode_ndarray(convergence["forward_DGs"]["magnitude"])
    reverse_dgs = decode_ndarray(convergence["reverse_DGs"]["magnitude"])
    hysteresis_profile = np.abs(forward_dgs - reverse_dgs)
    hysteresis_metrics = extract_hysteresis_metrics(fractions, hysteresis_profile)

    simulation_time_per_window_ns = quantity_to_ns(simulation_settings["production_length"])
    equilibration_time_per_window_ns = quantity_to_ns(
        simulation_settings["equilibration_length"]
    )
    total_time_per_window_ns = (
        simulation_time_per_window_ns + equilibration_time_per_window_ns
    )
    lambda_windows = int(lambda_settings["lambda_windows"])

    real_time_entry = load_last_real_time_analysis(result_path.with_suffix(""))
    analysis_start = _parse_datetime(analysis_result["start_time"])
    analysis_end = _parse_datetime(analysis_result["end_time"])

    structural_metrics = load_structural_metrics(results_root, analysis_result)
    perturbed_heavy_atom_metrics = load_perturbed_heavy_atom_metrics(
        result_path,
        setup_result,
        results_root,
    )
    mbar_offdiag_metrics = extract_adjacent_diagonal_metrics("mbar", mbar_matrix)
    repex_offdiag_metrics = extract_adjacent_diagonal_metrics(
        "replica_exchange",
        repex_matrix,
    )

    return {
        "transformation_id": transformation_id,
        "transformation": f"{ligand_i} -> {ligand_j}",
        "ligand_i": ligand_i,
        "ligand_j": ligand_j,
        "leg": leg_i,
        "repeat": repeat_name,
        "repeat_index": repeat_index,
        "repeat_id": analysis_result["outputs"]["repeat_id"],
        "generation": analysis_result["outputs"]["generation"],
        "result_json": str(result_path),
        "transformation_dir": str(result_path.with_suffix("")),
        "analysis_start": analysis_start.isoformat(),
        "analysis_end": analysis_end.isoformat(),
        "analysis_walltime_hours": _elapsed_hours(analysis_start, analysis_end),
        "dg_kcal_mol": quantity_to_float(analysis_result["outputs"]["unit_estimate"]),
        "dg_error_kcal_mol": quantity_to_float(
            analysis_result["outputs"]["unit_estimate_error"]
        ),
        "n_lambda_windows": lambda_windows,
        "simulation_time_ns_per_window": simulation_time_per_window_ns,
        "equilibration_time_ns_per_window": equilibration_time_per_window_ns,
        "total_time_ns_per_window": total_time_per_window_ns,
        "planned_total_sampling_ns": lambda_windows * total_time_per_window_ns,
        "mbar_overlap_scalar": float(unit_mbar["scalar"]),
        "mbar_second_eigenvalue": _second_largest_eigenvalue(mbar_eigenvalues),
        "replica_exchange_second_eigenvalue": _second_largest_eigenvalue(
            repex_eigenvalues
        ),
        "replica_exchange_spectral_gap": float(
            1.0 - _second_largest_eigenvalue(repex_eigenvalues)
        ),
        "hysteresis_profile_fractions": _serialize_vector(fractions),
        "hysteresis_profile_kcal_mol": _serialize_vector(hysteresis_profile),
        "reported_production_iterations": float(
            analysis_result["outputs"]["production_iterations"]
        ),
        "reported_equilibration_iterations": float(
            analysis_result["outputs"]["equilibration_iterations"]
        ),
        "completed_iterations": _maybe_nested_get(real_time_entry, "iteration"),
        "percent_complete": _maybe_nested_get(real_time_entry, "percent_complete"),
        "throughput_ns_per_day": _maybe_nested_get(
            real_time_entry, "timing_data", "ns_per_day"
        ),
        "average_seconds_per_iteration": _maybe_nested_get(
            real_time_entry, "timing_data", "average_seconds_per_iteration"
        ),
        "number_of_uncorrelated_samples": _maybe_nested_get(
            real_time_entry, "mbar_analysis", "number_of_uncorrelated_samples"
        ),
        "n_equilibrium_iterations_detected": _maybe_nested_get(
            real_time_entry, "mbar_analysis", "n_equilibrium_iterations"
        ),
        "statistical_inefficiency": _maybe_nested_get(
            real_time_entry, "mbar_analysis", "statistical_inefficiency"
        ),
        "real_time_standard_error_kT": _maybe_nested_get(
            real_time_entry, "mbar_analysis", "standard_error_in_kT"
        ),
        **hysteresis_metrics,
        **mbar_offdiag_metrics,
        **repex_offdiag_metrics,
        **perturbed_heavy_atom_metrics,
        **structural_metrics,
        "simulation_checkpoint": str(simulation_result["outputs"]["checkpoint"]["path"]),
    }


def build_repeat_results(leg_results: pd.DataFrame) -> pd.DataFrame:
    if leg_results.empty:
        return pd.DataFrame()

    base_columns = {
        "transformation_id",
        "transformation",
        "ligand_i",
        "ligand_j",
        "repeat",
        "repeat_index",
        "leg",
    }
    repeat_records: list[dict[str, Any]] = []

    grouped = leg_results.groupby(
        ["transformation_id", "transformation", "ligand_i", "ligand_j", "repeat", "repeat_index"],
        sort=True,
        dropna=False,
    )

    for group_key, group_frame in grouped:
        transformation_id, transformation, ligand_i, ligand_j, repeat, repeat_index = (
            group_key
        )
        row: dict[str, Any] = {
            "transformation_id": transformation_id,
            "transformation": transformation,
            "ligand_i": ligand_i,
            "ligand_j": ligand_j,
            "repeat": repeat,
            "repeat_index": repeat_index,
            "n_legs_found": int(len(group_frame)),
            "legs_found": ",".join(sorted(group_frame["leg"].tolist())),
        }
        for shared_column in (
            "ligand_i_heavy_atoms",
            "ligand_j_heavy_atoms",
            "mapped_heavy_atoms_core",
            "perturbed_heavy_atoms_ligand_i",
            "perturbed_heavy_atoms_ligand_j",
            "perturbed_heavy_atoms_total",
        ):
            shared_values = group_frame[shared_column].dropna().unique().tolist()
            if len(shared_values) == 1:
                row[shared_column] = shared_values[0]

        by_leg = {
            record["leg"]: record for record in group_frame.to_dict(orient="records")
        }
        complex_record = by_leg.get("complex")
        solvent_record = by_leg.get("solvent")

        if complex_record and solvent_record:
            row["ddg_kcal_mol"] = (
                complex_record["dg_kcal_mol"] - solvent_record["dg_kcal_mol"]
            )
            row["ddg_error_kcal_mol"] = math.sqrt(
                complex_record["dg_error_kcal_mol"] ** 2
                + solvent_record["dg_error_kcal_mol"] ** 2
            )
            row["min_leg_mbar_offdiag_min"] = min(
                complex_record["mbar_offdiag_min"], solvent_record["mbar_offdiag_min"]
            )
            row["min_leg_replica_exchange_offdiag_min"] = min(
                complex_record["replica_exchange_offdiag_min"],
                solvent_record["replica_exchange_offdiag_min"],
            )
            row["max_leg_hysteresis_50pct_kcal_mol"] = max(
                complex_record["hysteresis_50pct_kcal_mol"],
                solvent_record["hysteresis_50pct_kcal_mol"],
            )
        else:
            row["ddg_kcal_mol"] = np.nan
            row["ddg_error_kcal_mol"] = np.nan
            row["min_leg_mbar_offdiag_min"] = np.nan
            row["min_leg_replica_exchange_offdiag_min"] = np.nan
            row["max_leg_hysteresis_50pct_kcal_mol"] = np.nan

        for leg_name, record in sorted(by_leg.items()):
            for column, value in record.items():
                if column in base_columns:
                    continue
                row[f"{leg_name}_{column}"] = value

        repeat_records.append(row)

    return _records_to_frame(
        repeat_records,
        sort_columns=["repeat_index", "ligand_i", "ligand_j"],
    )


def build_leg_summary(leg_results: pd.DataFrame) -> pd.DataFrame:
    if leg_results.empty:
        return pd.DataFrame()

    summary_records: list[dict[str, Any]] = []
    grouped = leg_results.groupby(
        ["transformation_id", "transformation", "ligand_i", "ligand_j", "leg"],
        sort=True,
        dropna=False,
    )

    for (transformation_id, transformation, ligand_i, ligand_j, leg), frame in grouped:
        summary_records.append(
            {
                "transformation_id": transformation_id,
                "transformation": transformation,
                "ligand_i": ligand_i,
                "ligand_j": ligand_j,
                "leg": leg,
                "n_repeats": int(len(frame)),
                "dg_mean_kcal_mol": frame["dg_kcal_mol"].mean(),
                "dg_std_kcal_mol": frame["dg_kcal_mol"].std(ddof=1),
                "dg_sem_kcal_mol": frame["dg_kcal_mol"].sem(ddof=1),
                "dg_error_mean_kcal_mol": frame["dg_error_kcal_mol"].mean(),
                "mbar_overlap_scalar_mean": frame["mbar_overlap_scalar"].mean(),
                "mbar_offdiag_min_mean": frame["mbar_offdiag_min"].mean(),
                "mbar_offdiag_max_mean": frame["mbar_offdiag_max"].mean(),
                "replica_exchange_offdiag_min_mean": frame[
                    "replica_exchange_offdiag_min"
                ].mean(),
                "replica_exchange_offdiag_max_mean": frame[
                    "replica_exchange_offdiag_max"
                ].mean(),
                "hysteresis_50pct_mean_kcal_mol": frame[
                    "hysteresis_50pct_kcal_mol"
                ].mean(),
                "hysteresis_70pct_mean_kcal_mol": frame[
                    "hysteresis_70pct_kcal_mol"
                ].mean(),
                "hysteresis_90pct_mean_kcal_mol": frame[
                    "hysteresis_90pct_kcal_mol"
                ].mean(),
                "throughput_ns_per_day_mean": frame["throughput_ns_per_day"].mean(),
            }
        )

    return _records_to_frame(
        summary_records,
        sort_columns=["ligand_i", "ligand_j", "leg"],
    )


def build_transformation_summary(repeat_results: pd.DataFrame) -> pd.DataFrame:
    if repeat_results.empty:
        return pd.DataFrame()

    summary_records: list[dict[str, Any]] = []
    grouped = repeat_results.groupby(
        ["transformation_id", "transformation", "ligand_i", "ligand_j"],
        sort=True,
        dropna=False,
    )

    for (transformation_id, transformation, ligand_i, ligand_j), frame in grouped:
        summary_records.append(
            {
                "ligand_i": ligand_i,
                "ligand_j": ligand_j,
                "n_repeats": int(len(frame)),
                "ddg_mean_kcal_mol": frame["ddg_kcal_mol"].mean(),
                "ddg_std_kcal_mol": frame["ddg_kcal_mol"].std(ddof=1),
                "ddg_sem_kcal_mol": frame["ddg_kcal_mol"].sem(ddof=1),
                "ddg_error_mean_kcal_mol": frame["ddg_error_kcal_mol"].mean(),
                "complex_dg_mean_kcal_mol": frame["complex_dg_kcal_mol"].mean(),
                "solvent_dg_mean_kcal_mol": frame["solvent_dg_kcal_mol"].mean(),
                "perturbed_heavy_atoms_total": _first_non_null(
                    frame["perturbed_heavy_atoms_total"]
                ),
                "n_lambda_windows": _first_non_null(frame["complex_n_lambda_windows"]),
                "simulation_time_ns_per_window": _first_non_null(
                    frame["complex_simulation_time_ns_per_window"]
                ),
                "complex_mbar_offdiag_min": frame["complex_mbar_offdiag_min"].min(),
                "solvent_mbar_offdiag_min": frame["solvent_mbar_offdiag_min"].min(),
                "complex_replica_exchange_offdiag_min": frame[
                    "complex_replica_exchange_offdiag_min"
                ].min(),
                "solvent_replica_exchange_offdiag_min": frame[
                    "solvent_replica_exchange_offdiag_min"
                ].min(),
                "complex_hysteresis_50pct_kcal_mol": frame[
                    "complex_hysteresis_50pct_kcal_mol"
                ].max(),
                "solvent_hysteresis_50pct_kcal_mol": frame[
                    "solvent_hysteresis_50pct_kcal_mol"
                ].max(),
                "complex_hysteresis_70pct_kcal_mol": frame[
                    "complex_hysteresis_70pct_kcal_mol"
                ].max(),
                "solvent_hysteresis_70pct_kcal_mol": frame[
                    "solvent_hysteresis_70pct_kcal_mol"
                ].max(),
                "complex_hysteresis_90pct_kcal_mol": frame[
                    "complex_hysteresis_90pct_kcal_mol"
                ].max(),
                "solvent_hysteresis_90pct_kcal_mol": frame[
                    "solvent_hysteresis_90pct_kcal_mol"
                ].max(),
                "min_leg_mbar_offdiag_min": frame["min_leg_mbar_offdiag_min"].min(),
                "min_leg_replica_exchange_offdiag_min": frame[
                    "min_leg_replica_exchange_offdiag_min"
                ].min(),
                "max_leg_hysteresis_50pct_kcal_mol": frame[
                    "max_leg_hysteresis_50pct_kcal_mol"
                ].max(),
            }
        )

    return _records_to_frame(
        summary_records,
        sort_columns=["ligand_i", "ligand_j"],
    )


def build_summary_energies(transformation_summary: pd.DataFrame) -> pd.DataFrame:
    if transformation_summary.empty:
        return pd.DataFrame()

    columns = [
        "ligand_i",
        "ligand_j",
        "ddg_mean_kcal_mol",
        "ddg_std_kcal_mol",
        "ddg_sem_kcal_mol",
        "ddg_error_mean_kcal_mol",
        "complex_dg_mean_kcal_mol",
        "solvent_dg_mean_kcal_mol",
    ]
    return transformation_summary.loc[:, columns].copy()


def build_summary_convergence(transformation_summary: pd.DataFrame) -> pd.DataFrame:
    if transformation_summary.empty:
        return pd.DataFrame()

    columns = [
        "ligand_i",
        "ligand_j",
        "n_repeats",
        "n_lambda_windows",
        "simulation_time_ns_per_window",
        "perturbed_heavy_atoms_total",
        "complex_mbar_offdiag_min",
        "solvent_mbar_offdiag_min",
        "complex_replica_exchange_offdiag_min",
        "solvent_replica_exchange_offdiag_min",
        "complex_hysteresis_50pct_kcal_mol",
        "complex_hysteresis_70pct_kcal_mol",
        "complex_hysteresis_90pct_kcal_mol",
        "solvent_hysteresis_50pct_kcal_mol",
        "solvent_hysteresis_70pct_kcal_mol",
        "solvent_hysteresis_90pct_kcal_mol",
        "min_leg_mbar_offdiag_min",
        "min_leg_replica_exchange_offdiag_min",
        "max_leg_hysteresis_50pct_kcal_mol",
    ]
    return transformation_summary.loc[:, columns].copy()


def plot_summary_convergence_heatmap(
    summary_convergence: pd.DataFrame,
    output_path: str | Path,
) -> None:
    if summary_convergence.empty:
        return

    plot_columns = [
        column
        for _, columns in QUALITY_COLUMN_GROUPS
        for column in columns
    ]
    missing_columns = [column for column in plot_columns if column not in summary_convergence]
    if missing_columns:
        raise ValueError(
            "Missing columns required for convergence heatmap: "
            + ", ".join(missing_columns)
        )

    row_labels = [
        f"{row.ligand_i} \u2192 {row.ligand_j}"
        for row in summary_convergence[["ligand_i", "ligand_j"]].itertuples(index=False)
    ]
    column_labels = [QUALITY_METRIC_SPECS[column]["label"] for column in plot_columns]

    color_grid = np.array(
        [
            [
                to_rgb(classify_quality(summary_convergence.iloc[row_idx][column], column))
                for column in plot_columns
            ]
            for row_idx in range(len(summary_convergence))
        ]
    )
    text_grid = [
        [
            format_quality_value(summary_convergence.iloc[row_idx][column], column)
            for column in plot_columns
        ]
        for row_idx in range(len(summary_convergence))
    ]

    fig_width = max(9.8, 0.36 * len(plot_columns) + 6.4)
    fig_height = max(4.2, 0.26 * len(row_labels) + 1.8)
    fig, (ax, info_ax) = plt.subplots(
        ncols=2,
        figsize=(fig_width, fig_height),
        gridspec_kw={"width_ratios": [3.7, 1.6]},
        constrained_layout=True,
    )
    ax.imshow(color_grid, aspect="auto")

    ax.set_xticks(np.arange(len(column_labels)))
    ax.set_xticklabels(column_labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=9)
    ax.set_title("Quality Summary", fontsize=12, pad=25)
    ax.set_xticks(np.arange(-0.5, len(column_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for row_idx, row_text in enumerate(text_grid):
        for col_idx, text in enumerate(row_text):
            column = plot_columns[col_idx]
            ax.text(
                col_idx,
                row_idx,
                text,
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold" if column in EMPHASIZED_QUALITY_COLUMNS else "normal",
            )

    group_boundaries: list[float] = []
    group_centers: list[tuple[float, str]] = []
    start_index = 0
    for group_name, group_columns in QUALITY_COLUMN_GROUPS:
        end_index = start_index + len(group_columns)
        if start_index > 0:
            group_boundaries.append(start_index - 0.5)
        center = (start_index + end_index - 1) / 2.0
        group_centers.append((center, group_name))
        y_line = -1.05
        ax.plot(
            [start_index - 0.35, end_index - 0.65],
            [y_line, y_line],
            color="#495057",
            linewidth=1.4,
            clip_on=False,
        )
        ax.plot(
            [start_index - 0.35, start_index - 0.35],
            [y_line, y_line + 0.18],
            color="#495057",
            linewidth=1.4,
            clip_on=False,
        )
        ax.plot(
            [end_index - 0.65, end_index - 0.65],
            [y_line, y_line + 0.18],
            color="#495057",
            linewidth=1.4,
            clip_on=False,
        )
        start_index = end_index

    for boundary in group_boundaries:
        ax.plot(
            [boundary, boundary],
            [-0.5, len(row_labels) - 0.5],
            color="#adb5bd",
            linewidth=2.2,
            clip_on=False,
        )

    for center, group_name in group_centers:
        ax.text(
            center,
            -1.45,
            group_name,
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
            clip_on=False,
        )

    info_ax.axis("off")
    reference_rows = [
        [
            metric_label,
            QUALITY_RULES[rule_key]["good_label"],
            QUALITY_RULES[rule_key]["ok_label"],
            QUALITY_RULES[rule_key]["bad_label"],
        ]
        for metric_label, rule_key in QUALITY_REFERENCE_RULES
    ]
    reference_colors = [
        [
            "#e9ecef",
            "#e9ecef",
            "#e9ecef",
            "#e9ecef",
        ]
        for row in reference_rows
    ]
    reference_table = info_ax.table(
        cellText=reference_rows,
        cellColours=reference_colors,
        colLabels=["Metric", "GOOD", "OK", "BAD"],
        colColours=["#dee2e6", QUALITY_COLORS["good"], QUALITY_COLORS["ok"], QUALITY_COLORS["bad"]],
        cellLoc="center",
        colLoc="center",
        bbox=[0.0, 0.4, 0.90, 0.5],
        colWidths=[0.40, 0.20, 0.20, 0.20],
    )
    reference_table.auto_set_font_size(False)
    reference_table.set_fontsize(6)
    reference_table.scale(1.0, 0.5)
    for (row_idx, col_idx), cell in reference_table.get_celld().items():
        cell.set_edgecolor("#adb5bd")
        cell.set_linewidth(0.5)
        cell.PAD = 0.1
        cell.get_text().set_wrap(True)
        if row_idx == 0:
            cell.get_text().set_fontweight("bold")
        if col_idx == 0:
            cell.get_text().set_ha("left")
        if row_idx > 0 and col_idx == 2 and not reference_rows[row_idx - 1][2]:
            cell.get_text().set_text("")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", format="png")
    plt.close(fig)

def classify_quality(value: Any, column: str) -> str:
    return QUALITY_COLORS[classify_quality_status(value, column)]


def classify_quality_status(value: Any, column: str) -> str:
    if pd.isna(value):
        return "missing"

    spec = QUALITY_METRIC_SPECS[column]
    rule = QUALITY_RULES[spec["rule"]]
    if rule is None:
        return "missing"

    direction = rule["direction"]
    good = rule["good"]
    ok = rule["ok"]

    if direction == "lower_is_better":
        if value <= good:
            return "good"
        if value <= ok:
            return "ok"
        return "bad"

    if direction == "higher_is_better":
        if value >= good:
            return "good"
        if value >= ok:
            return "ok"
        return "bad"

    raise ValueError(f"Unknown quality direction: {direction}")


def format_quality_value(value: Any, column: str) -> str:
    if pd.isna(value):
        return "NA"
    return QUALITY_METRIC_SPECS[column]["format"].format(value)


def decode_ndarray(node: dict[str, Any]) -> np.ndarray:
    raw = node["bytes"]["latin-1"].encode("latin-1")
    raw = _maybe_zstd_decompress(raw)
    dtype = np.dtype(node["dtype"])
    expected_nbytes = dtype.itemsize * math.prod(node["shape"])

    if len(raw) == expected_nbytes:
        array = np.frombuffer(raw, dtype=dtype)
        return array.reshape(tuple(node["shape"]))

    loaded = np.load(io.BytesIO(raw), allow_pickle=False)
    return np.asarray(loaded).reshape(tuple(node["shape"]))


def quantity_to_float(quantity: dict[str, Any]) -> float:
    magnitude = quantity["magnitude"]
    if isinstance(magnitude, dict):
        if magnitude.get("__class__") == "ndarray":
            raise TypeError("Expected scalar quantity but received an ndarray quantity")
        raise TypeError("Expected scalar quantity but received a mapping magnitude")
    return float(magnitude)


def quantity_to_ns(quantity: dict[str, Any]) -> float:
    unit = quantity["unit"]
    if unit not in TIME_UNIT_TO_NS:
        raise ValueError(f"Unsupported time unit for conversion to ns: {unit}")
    return quantity_to_float(quantity) * TIME_UNIT_TO_NS[unit]


def extract_hysteresis_metrics(
    fractions: np.ndarray,
    hysteresis_profile: np.ndarray,
) -> dict[str, float]:
    metrics: dict[str, float] = {}

    for target in TARGET_HYSTERESIS_FRACTIONS:
        index = int(np.argmin(np.abs(fractions - target)))
        percent = int(round(target * 100))
        metrics[f"hysteresis_{percent}pct_kcal_mol"] = float(hysteresis_profile[index])
        metrics[f"hysteresis_{percent}pct_fraction_actual"] = float(fractions[index])

    metrics["final_hysteresis_kcal_mol"] = float(hysteresis_profile[-1])
    return metrics


def extract_adjacent_diagonal_metrics(
    prefix: str,
    matrix: np.ndarray,
) -> dict[str, float | str]:
    upper = np.diag(matrix, k=1)
    lower = np.diag(matrix, k=-1)
    adjacent = np.concatenate([upper, lower])

    return {
        f"{prefix}_offdiag_values": _serialize_vector(adjacent),
        f"{prefix}_offdiag_min": float(np.min(adjacent)),
        f"{prefix}_offdiag_max": float(np.max(adjacent)),
        f"{prefix}_offdiag_mean": float(np.mean(adjacent)),
        f"{prefix}_offdiag_upper_values": _serialize_vector(upper),
        f"{prefix}_offdiag_lower_values": _serialize_vector(lower),
    }


def load_last_real_time_analysis(transformation_dir: Path) -> dict[str, Any]:
    yaml_paths = sorted(
        transformation_dir.glob(
            "shared_HybridTopologyMultiStateSimulationUnit-*/simulation_real_time_analysis.yaml"
        )
    )
    if not yaml_paths:
        return {}

    entries = yaml.safe_load(yaml_paths[0].read_text())
    if not entries:
        return {}

    return entries[-1]


def load_structural_metrics(
    results_root: Path,
    analysis_result: dict[str, Any],
) -> dict[str, float]:
    structural_path = Path(analysis_result["outputs"]["structural_analysis"]["path"])
    full_path = results_root.parent / structural_path
    if not full_path.exists():
        return {}

    metrics: dict[str, float] = {}
    with np.load(full_path) as archive:
        for key in ("protein_RMSD", "ligand_RMSD", "ligand_COM_drift", "protein_2D_RMSD"):
            if key not in archive.files:
                continue
            values = np.asarray(archive[key], dtype=float)
            if values.size == 0:
                continue
            metrics[f"{key}_mean"] = float(values.mean())
            metrics[f"{key}_max"] = float(values.max())

    return metrics


def load_perturbed_heavy_atom_metrics(
    result_path: Path,
    setup_result: dict[str, Any],
    results_root: Path,
) -> dict[str, int]:
    transformation_path = (
        results_root.parent / "network_setup" / "transformations" / result_path.name
    )
    if not transformation_path.exists():
        return {}

    atom_mapping = load_ligand_atom_mapping(str(transformation_path))

    ligand_a_atoms = setup_result["inputs"]["stateA"]["components"]["ligand"]["atoms"]
    ligand_b_atoms = setup_result["inputs"]["stateB"]["components"]["ligand"]["atoms"]
    heavy_a = {index for index, atom in enumerate(ligand_a_atoms) if atom[0] > 1}
    heavy_b = {index for index, atom in enumerate(ligand_b_atoms) if atom[0] > 1}

    reverse_mapping = {target: source for source, target in atom_mapping.items()}
    perturbed_a = {
        index
        for index in heavy_a
        if index not in atom_mapping
        or ligand_a_atoms[index][0] != ligand_b_atoms[atom_mapping[index]][0]
    }
    perturbed_b = {
        index
        for index in heavy_b
        if index not in reverse_mapping
        or ligand_b_atoms[index][0] != ligand_a_atoms[reverse_mapping[index]][0]
    }
    mapped_heavy_core = heavy_a - perturbed_a

    return {
        "ligand_i_heavy_atoms": len(heavy_a),
        "ligand_j_heavy_atoms": len(heavy_b),
        "mapped_heavy_atoms_core": len(mapped_heavy_core),
        "perturbed_heavy_atoms_ligand_i": len(perturbed_a),
        "perturbed_heavy_atoms_ligand_j": len(perturbed_b),
        "perturbed_heavy_atoms_total": len(perturbed_a) + len(perturbed_b),
    }


@lru_cache(maxsize=None)
def load_ligand_atom_mapping(transformation_path: str) -> dict[int, int]:
    items = json.loads(Path(transformation_path).read_text())
    lookup = {key: value for key, value in items}
    transformation = next(
        value for key, value in lookup.items() if key.startswith("Transformation-")
    )
    mapping_key = transformation["mapping"][":gufe-key:"]
    mapping = lookup[mapping_key]["componentA_to_componentB"]
    return {int(source): int(target) for source, target in mapping.items()}


def _find_analysis_result(payload: dict[str, Any]) -> dict[str, Any]:
    for result_list in payload["protocol_result"]["data"].values():
        for result in result_list:
            if result["name"].startswith("HybridTopology Analysis"):
                return result

    for result in payload["unit_results"].values():
        if result["name"].startswith("HybridTopology Analysis"):
            return result

    raise ValueError("Could not find HybridTopology Analysis result in payload")


def _split_state_name(state_name: str) -> tuple[str, str]:
    ligand, leg = state_name.rsplit("_", 1)
    return ligand, leg


def _parse_repeat_folder(repeat_name: str) -> tuple[str, int]:
    match = re.fullmatch(r"repeat(\d+)", repeat_name)
    if not match:
        raise ValueError(f"Unexpected repeat folder name: {repeat_name}")
    return repeat_name, int(match.group(1))


def _second_largest_eigenvalue(values: np.ndarray) -> float:
    sorted_values = np.sort(np.real(values))
    if sorted_values.size < 2:
        return float("nan")
    return float(sorted_values[-2])


def _maybe_zstd_decompress(raw: bytes) -> bytes:
    if not raw.startswith(ZSTD_MAGIC):
        return raw

    candidate_binaries = [
        shutil.which("zstd"),
        shutil.which("unzstd"),
        str((Path(sys.executable).resolve().parent / "zstd")),
        str((Path(sys.executable).resolve().parent / "unzstd")),
        str((Path(sys.prefix).resolve() / "bin" / "zstd")),
        str((Path(sys.prefix).resolve() / "bin" / "unzstd")),
    ]
    zstd_binary = next((path for path in candidate_binaries if path and Path(path).exists()), None)
    if zstd_binary is None:
        raise RuntimeError(
            "Encountered zstd-compressed openFE arrays but could not find a `zstd` binary."
        )

    return subprocess.run(
        [zstd_binary, "-d", "-c"],
        input=raw,
        capture_output=True,
        check=True,
    ).stdout


def _serialize_vector(values: np.ndarray) -> str:
    return json.dumps([round(float(value), 6) for value in values.tolist()])


def _parse_datetime(node: dict[str, Any]) -> datetime:
    return datetime.fromisoformat(node["isotime"])


def _elapsed_hours(start: datetime, end: datetime) -> float:
    return (end - start).total_seconds() / 3600.0


def _maybe_nested_get(mapping: dict[str, Any], *keys: str) -> float | None:
    current: Any = mapping
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return None
        current = current[key]
    return current


def _first_non_null(series: pd.Series) -> Any:
    non_null = series.dropna()
    if non_null.empty:
        return np.nan
    return non_null.iloc[0]


def _records_to_frame(
    records: list[dict[str, Any]],
    sort_columns: list[str],
) -> pd.DataFrame:
    frame = pd.DataFrame.from_records(records)
    if frame.empty:
        return frame
    return frame.sort_values(sort_columns).reset_index(drop=True)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract per-leg and per-transformation workup metrics from openFE "
            "relative free energy results."
        )
    )
    parser.add_argument(
        "results_dir",
        nargs="?",
        default="example/results",
        help="Directory containing repeat subdirectories with openFE result JSON files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory where CSV summaries should be written. Defaults to <results_dir>/../workup.",
    )
    parser.add_argument(
        "--no-export",
        action="store_true",
        help="Load and summarize results without writing CSV files.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Write CSV summaries but skip plot generation.",
    )
    return parser


def main() -> None:
    parser = build_argument_parser()
    args = parser.parse_args()

    workup = load_results(args.results_dir)

    print(
        "Loaded "
        f"{len(workup.leg_results)} leg rows, "
        f"{len(workup.repeat_results)} transformation-repeat rows, and "
        f"{len(workup.transformation_summary)} aggregated transformation rows."
    )

    if not workup.transformation_summary.empty:
        preview = workup.transformation_summary[
            ["ligand_i", "ligand_j", "ddg_mean_kcal_mol", "ddg_std_kcal_mol"]
        ]
        print()
        print("DDG summary preview:")
        print(preview.to_string(index=False))

    if args.no_export:
        return

    output_dir = args.output_dir or (Path(args.results_dir).resolve().parent / "workup")
    written_paths = workup.export_csv(output_dir)
    plot_paths = {} if args.no_plots else workup.export_plots(output_dir)

    print()
    print(f"Wrote CSV files to {output_dir}")
    for filename, path in written_paths.items():
        print(f"  {filename}: {path}")
    for filename, path in plot_paths.items():
        print(f"  {filename}: {path}")


if __name__ == "__main__":
    main()
