#!/usr/bin/env python3

from __future__ import annotations

import argparse
import io
import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence, cast

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, OffsetImage, TextArea, VPacker
import numpy as np
import pandas as pd


QUALITY_STATUS_ORDER = {
    "good": 0,
    "ok": 1,
    "bad": 2,
}
NETWORK_QUALITY_RULES = {
    "min_leg_mbar_offdiag_min": {
        "direction": "higher_is_better",
        "good": 0.10,
        "ok": 0.03,
    },
    "min_leg_replica_exchange_offdiag_min": {
        "direction": "higher_is_better",
        "good": 0.20,
        "ok": 0.10,
    },
}
DEFAULT_QUALITY_COLORS = {
    "good": "#b7e4c7",
    "ok": "#ffd6a5",
    "bad": "#ffadad",
    "missing": "#ffffff",
}
DEFAULT_EDGE_QUALITY_COLUMNS = (
    "min_leg_mbar_offdiag_min",
    "min_leg_replica_exchange_offdiag_min",
)
DEFAULT_NETWORK_EDGE_MISSING_COLOR = "#adb5bd"


def plot_summary_ligand_network(
    transformation_summary: pd.DataFrame,
    output_path: str | Path,
    results_root: Path | None = None,
    *,
    quality_colors: Mapping[str, str] | None = None,
    edge_quality_columns: Sequence[str] | None = None,
    classify_quality_status: Callable[[Any, str], str] | None = None,
    network_edge_missing_color: str = DEFAULT_NETWORK_EDGE_MISSING_COLOR,
) -> None:
    if transformation_summary.empty:
        return

    quality_colors = quality_colors or DEFAULT_QUALITY_COLORS
    edge_quality_columns = edge_quality_columns or DEFAULT_EDGE_QUALITY_COLUMNS
    classify_quality_status = classify_quality_status or default_classify_quality_status

    nodes = sorted(
        {
            cast(str, ligand_name)
            for ligand_name in set(transformation_summary["ligand_i"]).union(
                transformation_summary["ligand_j"]
            )
        }
    )
    edges = [
        (cast(str, row.ligand_i), cast(str, row.ligand_j))
        for row in transformation_summary[["ligand_i", "ligand_j"]].itertuples(index=False)
    ]
    positions = compute_network_layout(nodes, edges)
    ligand_components = load_ligand_components(results_root, transformation_summary)
    node_images = build_aligned_ligand_node_images(ligand_components)
    node_positions = [positions[ligand_name] for ligand_name in nodes]

    fig_width = max(10.8, 2.45 * math.sqrt(len(nodes)) + 6.3)
    fig_height = max(8.2, 2.25 * math.sqrt(len(nodes)) + 4.4)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    ax.set_facecolor("#f8f9fa")
    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.28, 1.12)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Ligand Network With Transformation ΔΔG", fontsize=14, pad=16)

    node_zoom = min(0.33, max(0.17, 0.74 / max(math.sqrt(len(nodes)), 1.0)))
    node_radius = 0.12 + 0.07 * node_zoom
    label_positions: list[np.ndarray] = []
    edge_rows = list(transformation_summary.itertuples(index=False))
    edge_segments = [
        (positions[cast(str, row.ligand_i)], positions[cast(str, row.ligand_j)])
        for row in edge_rows
    ]
    edge_label_specs: list[tuple[np.ndarray, np.ndarray]] = []
    for edge_index, row in enumerate(edge_rows):
        ligand_i = cast(str, row.ligand_i)
        ligand_j = cast(str, row.ligand_j)
        edge_label_position, edge_anchor = choose_edge_label_position(
            positions[ligand_i],
            positions[ligand_j],
            node_positions=node_positions,
            existing_positions=label_positions,
            node_radius=node_radius,
            edge_segments=edge_segments,
            edge_index=edge_index,
        )
        edge_label_specs.append((edge_label_position, edge_anchor))
        label_positions.append(edge_label_position)
    edge_label_specs = relax_edge_label_positions(
        edge_label_specs,
        node_positions=node_positions,
        edge_segments=edge_segments,
        node_radius=node_radius,
    )

    for row, (edge_label, edge_anchor) in zip(edge_rows, edge_label_specs):
        ligand_i = cast(str, row.ligand_i)
        ligand_j = cast(str, row.ligand_j)
        start = positions[ligand_i]
        end = positions[ligand_j]
        start_xy = (float(start[0]), float(start[1]))
        end_xy = (float(end[0]), float(end[1]))
        quality_status = summarize_edge_quality(
            row,
            edge_quality_columns=edge_quality_columns,
            classify_quality_status=classify_quality_status,
        )
        edge_color = (
            quality_colors.get(quality_status, network_edge_missing_color)
            if quality_status != "missing"
            else network_edge_missing_color
        )
        ax.annotate(
            "",
            xy=end_xy,
            xytext=start_xy,
            arrowprops={
                "arrowstyle": "-|>",
                "color": edge_color,
                "linewidth": 2.5,
                "alpha": 0.96,
                "mutation_scale": 16,
                "shrinkA": 36,
                "shrinkB": 36,
                "linestyle": "--" if quality_status == "missing" else "-",
                "connectionstyle": "arc3,rad=0.04",
            },
            zorder=1,
        )
        ax.plot(
            [float(edge_anchor[0]), float(edge_label[0])],
            [float(edge_anchor[1]), float(edge_label[1])],
            color=edge_color,
            linewidth=0.9,
            alpha=0.72,
            zorder=2,
        )

        ax.text(
            float(edge_label[0]),
            float(edge_label[1]),
            format_ddg_label(row.ddg_mean_kcal_mol, row.ddg_std_kcal_mol),
            ha="center",
            va="center",
            fontsize=7.6,
            linespacing=0.96,
            bbox={
                "boxstyle": "round,pad=0.20",
                "facecolor": "#ffffff",
                "edgecolor": edge_color,
                "linewidth": 1.1,
                "alpha": 0.97,
            },
            zorder=3,
        )

    for ligand_name in nodes:
        position = positions[ligand_name]
        position_xy = (float(position[0]), float(position[1]))
        node_image = node_images.get(ligand_name)

        if node_image is not None:
            node_box = VPacker(
                children=[
                    OffsetImage(node_image, zoom=node_zoom),
                    TextArea(
                        format_ligand_display_name(ligand_name),
                        textprops={
                            "fontsize": 7.8,
                            "fontweight": "bold",
                            "color": "#212529",
                            "multialignment": "center",
                        },
                    ),
                ],
                align="center",
                pad=0,
                sep=2,
            )
            annotation = AnnotationBbox(
                node_box,
                position_xy,
                frameon=True,
                bboxprops={
                    "edgecolor": "#495057",
                    "facecolor": "#ffffff",
                    "linewidth": 1.2,
                    "boxstyle": "round,pad=0.32",
                },
                zorder=4,
            )
            ax.add_artist(annotation)
            continue

        ax.text(
            position_xy[0],
            position_xy[1],
            ligand_name,
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            bbox={
                "boxstyle": "round,pad=0.42",
                "facecolor": "#ffffff",
                "edgecolor": "#495057",
                "linewidth": 1.2,
            },
            zorder=4,
        )

    legend_handles = [
        Line2D([0], [0], color=quality_colors["good"], lw=3, label="Good"),
        Line2D([0], [0], color=quality_colors["ok"], lw=3, label="OK"),
        Line2D([0], [0], color=quality_colors["bad"], lw=3, label="Bad"),
        Line2D([0], [0], color=network_edge_missing_color, lw=3, ls="--", label="Missing"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper left",
        frameon=True,
        facecolor="#ffffff",
        edgecolor="#ced4da",
        title="Edge Quality",
        title_fontsize=10,
    )
    ax.text(
        0.01,
        0.02,
        (
            "Edge color reflects the worst classified overlap metric across "
            "MBAR overlap and HREX overlap."
        ),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.2,
        bbox={
            "boxstyle": "round,pad=0.38",
            "facecolor": "#ffffff",
            "edgecolor": "#ced4da",
        },
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", format="png")
    plt.close(fig)


def compute_network_layout(
    nodes: list[str],
    directed_edges: list[tuple[str, str]],
    iterations: int = 260,
) -> dict[str, np.ndarray]:
    if not nodes:
        return {}

    if len(nodes) == 1:
        return {nodes[0]: np.array([0.5, 0.5], dtype=float)}

    if len(nodes) == 2:
        return {
            nodes[0]: np.array([0.30, 0.50], dtype=float),
            nodes[1]: np.array([0.70, 0.50], dtype=float),
        }

    node_to_index = {node: index for index, node in enumerate(nodes)}
    undirected_edges = sorted(
        {
            tuple(sorted((source, target)))
            for source, target in directed_edges
            if source in node_to_index and target in node_to_index and source != target
        }
    )
    edge_indices = np.array(
        [(node_to_index[source], node_to_index[target]) for source, target in undirected_edges],
        dtype=int,
    )
    best_positions: np.ndarray | None = None
    best_score = float("-inf")
    base_angles = np.linspace(0.0, 2.0 * math.pi, len(nodes), endpoint=False)
    rng = np.random.default_rng(20260329)

    for restart in range(18):
        angle_offset = (2.0 * math.pi / max(len(nodes), 1)) * (restart / 18.0)
        jitter = rng.normal(scale=0.03, size=(len(nodes), 2))
        radius_x = 0.34 + 0.03 * math.sin(restart)
        radius_y = 0.30 + 0.03 * math.cos(restart)
        initial_positions = np.column_stack(
            (
                np.cos(base_angles + angle_offset) * radius_x + 0.5,
                np.sin(base_angles + angle_offset) * radius_y + 0.5,
            )
        )
        initial_positions = np.clip(initial_positions + jitter, 0.08, 0.92)
        candidate_positions = _run_force_layout(
            initial_positions,
            edge_indices=edge_indices,
            iterations=iterations,
        )
        candidate_positions = rescale_layout_positions(candidate_positions)
        candidate_score = score_network_layout(candidate_positions, edge_indices)
        if candidate_score > best_score:
            best_score = candidate_score
            best_positions = candidate_positions

    assert best_positions is not None
    return {node: point for node, point in zip(nodes, best_positions)}


def _run_force_layout(
    positions: np.ndarray,
    *,
    edge_indices: np.ndarray,
    iterations: int,
) -> np.ndarray:
    positions = positions.copy()
    k = math.sqrt(1.15 / len(positions))
    temperature = 0.16

    for _ in range(iterations):
        displacement = np.zeros_like(positions)

        delta = positions[:, None, :] - positions[None, :, :]
        distance = np.linalg.norm(delta, axis=2)
        distance = np.maximum(distance, 1e-4)
        np.fill_diagonal(distance, np.inf)
        repulsive_force = (k * k) / (distance * distance)
        displacement += np.sum(
            (delta / distance[:, :, None]) * repulsive_force[:, :, None],
            axis=1,
        )

        for source_index, target_index in edge_indices:
            edge_delta = positions[source_index] - positions[target_index]
            edge_distance = max(np.linalg.norm(edge_delta), 1e-4)
            direction = edge_delta / edge_distance
            attractive_force = ((edge_distance - 0.28) * abs(edge_distance - 0.28)) / k
            displacement[source_index] -= direction * attractive_force
            displacement[target_index] += direction * attractive_force

        displacement -= (positions - 0.5) * 0.05
        step = np.linalg.norm(displacement, axis=1)
        step[step == 0.0] = 1.0
        positions += displacement / step[:, None] * np.minimum(step, temperature)[:, None]
        positions = np.clip(positions, 0.06, 0.94)
        temperature *= 0.985

    return positions


def score_network_layout(positions: np.ndarray, edge_indices: np.ndarray) -> float:
    score = 0.0

    pair_distances = []
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            distance = float(np.linalg.norm(positions[i] - positions[j]))
            pair_distances.append(distance)
            score += 2.8 * distance
            if distance < 0.23:
                score -= 20.0 * (0.23 - distance)

    for source_index, target_index in edge_indices:
        edge_length = float(np.linalg.norm(positions[source_index] - positions[target_index]))
        score -= 0.55 * abs(edge_length - 0.30)

    segments = [
        (positions[source_index], positions[target_index])
        for source_index, target_index in edge_indices
    ]
    for i in range(len(segments)):
        for j in range(i + 1, len(segments)):
            if _segments_share_endpoint(edge_indices[i], edge_indices[j]):
                continue
            if segments_cross(segments[i][0], segments[i][1], segments[j][0], segments[j][1]):
                score -= 6.0

    for point in positions:
        margin = min(point[0], point[1], 1.0 - point[0], 1.0 - point[1])
        score += 0.8 * margin

    if pair_distances:
        score += 1.8 * min(pair_distances)

    return score


def rescale_layout_positions(positions: np.ndarray) -> np.ndarray:
    minimum = positions.min(axis=0)
    maximum = positions.max(axis=0)
    span = maximum - minimum

    normalized = np.empty_like(positions)
    for axis in range(positions.shape[1]):
        if span[axis] < 1e-8:
            normalized[:, axis] = 0.5
        else:
            normalized[:, axis] = (positions[:, axis] - minimum[axis]) / span[axis]

    return normalized * 0.76 + 0.12


def edge_label_position(start: np.ndarray, end: np.ndarray) -> tuple[float, float]:
    midpoint = (start + end) / 2.0
    direction = end - start
    distance = np.linalg.norm(direction)
    if distance < 1e-8:
        return float(midpoint[0]), float(midpoint[1])

    perpendicular = np.array([-direction[1], direction[0]], dtype=float) / distance
    if perpendicular[1] < 0.0:
        perpendicular *= -1.0
    label_position = midpoint + perpendicular * 0.035
    return float(label_position[0]), float(label_position[1])


def choose_edge_label_position(
    start: np.ndarray,
    end: np.ndarray,
    *,
    node_positions: Sequence[np.ndarray],
    existing_positions: Sequence[np.ndarray],
    node_radius: float,
    edge_segments: Sequence[tuple[np.ndarray, np.ndarray]],
    edge_index: int,
) -> tuple[np.ndarray, np.ndarray]:
    midpoint = (start + end) / 2.0
    direction = end - start
    distance = np.linalg.norm(direction)
    if distance < 1e-8:
        return midpoint.copy(), midpoint.copy()

    direction /= distance
    perpendicular = np.array([-direction[1], direction[0]], dtype=float)
    preferred_sign = 1.0 if perpendicular[1] >= 0.0 else -1.0

    candidates: list[tuple[np.ndarray, np.ndarray]] = []
    for fraction in (0.34, 0.45, 0.56, 0.67):
        anchor = start + (end - start) * fraction
        for sign in (preferred_sign, -preferred_sign):
            for offset in (0.06, 0.09, 0.12, 0.15, 0.18):
                candidates.append((anchor + perpendicular * offset * sign, anchor))

    best_candidate = midpoint + perpendicular * 0.10 * preferred_sign
    best_anchor = midpoint
    best_score = float("-inf")
    for candidate, anchor in candidates:
        if not (0.06 <= candidate[0] <= 0.94 and 0.06 <= candidate[1] <= 0.94):
            continue

        min_node_distance = min(np.linalg.norm(candidate - node) for node in node_positions)
        min_label_distance = (
            min(np.linalg.norm(candidate - label) for label in existing_positions)
            if existing_positions
            else 1.0
        )
        endpoint_distance = min(
            np.linalg.norm(candidate - start),
            np.linalg.norm(candidate - end),
        )
        min_other_edge_distance = min(
            point_to_segment_distance(candidate, other_start, other_end)
            for idx, (other_start, other_end) in enumerate(edge_segments)
            if idx != edge_index
        )
        score = (
            3.2 * min_node_distance
            + 2.4 * min_label_distance
            + 1.8 * min_other_edge_distance
            + 0.5 * endpoint_distance
        )
        if min_node_distance < node_radius:
            score -= 10.0
        if min_label_distance < 0.06:
            score -= 8.0
        if min_other_edge_distance < 0.04:
            score -= 6.0
        if score > best_score:
            best_score = score
            best_candidate = candidate
            best_anchor = anchor

    return best_candidate, best_anchor


def relax_edge_label_positions(
    label_specs: Sequence[tuple[np.ndarray, np.ndarray]],
    *,
    node_positions: Sequence[np.ndarray],
    edge_segments: Sequence[tuple[np.ndarray, np.ndarray]],
    node_radius: float,
    iterations: int = 140,
) -> list[tuple[np.ndarray, np.ndarray]]:
    if not label_specs:
        return []

    positions = [label.copy() for label, _ in label_specs]
    anchors = [anchor.copy() for _, anchor in label_specs]
    home_positions = [label.copy() for label, _ in label_specs]

    for _ in range(iterations):
        displacements = [np.zeros(2, dtype=float) for _ in positions]

        for i, position in enumerate(positions):
            for j in range(i + 1, len(positions)):
                delta = position - positions[j]
                distance = float(np.linalg.norm(delta))
                if distance < 1e-6:
                    delta = np.array([0.001, 0.0], dtype=float)
                    distance = 0.001
                if distance < 0.11:
                    force = (0.11 - distance) * 0.22
                    direction = delta / distance
                    displacements[i] += direction * force
                    displacements[j] -= direction * force

            for node_position in node_positions:
                delta = position - node_position
                distance = float(np.linalg.norm(delta))
                minimum_distance = node_radius + 0.05
                if distance < 1e-6:
                    delta = np.array([0.001, 0.0], dtype=float)
                    distance = 0.001
                if distance < minimum_distance:
                    force = (minimum_distance - distance) * 0.18
                    displacements[i] += (delta / distance) * force

            for edge_index, (edge_start, edge_end) in enumerate(edge_segments):
                if edge_index == i:
                    continue
                distance = point_to_segment_distance(position, edge_start, edge_end)
                if distance < 0.035:
                    projection = closest_point_on_segment(position, edge_start, edge_end)
                    delta = position - projection
                    norm = float(np.linalg.norm(delta))
                    if norm < 1e-6:
                        continue
                    displacements[i] += (delta / norm) * (0.035 - distance) * 0.12

            displacements[i] += (home_positions[i] - position) * 0.045
            displacements[i] += (anchors[i] - position) * 0.01

        for i, displacement in enumerate(displacements):
            step = float(np.linalg.norm(displacement))
            if step > 0.018:
                displacement *= 0.018 / step
            positions[i] = np.clip(positions[i] + displacement, 0.05, 0.95)

    return list(zip(positions, anchors))


def point_to_segment_distance(point: np.ndarray, start: np.ndarray, end: np.ndarray) -> float:
    segment = end - start
    length_sq = float(np.dot(segment, segment))
    if length_sq < 1e-12:
        return float(np.linalg.norm(point - start))

    t = float(np.dot(point - start, segment) / length_sq)
    t = max(0.0, min(1.0, t))
    projection = start + t * segment
    return float(np.linalg.norm(point - projection))


def closest_point_on_segment(point: np.ndarray, start: np.ndarray, end: np.ndarray) -> np.ndarray:
    segment = end - start
    length_sq = float(np.dot(segment, segment))
    if length_sq < 1e-12:
        return start.copy()

    t = float(np.dot(point - start, segment) / length_sq)
    t = max(0.0, min(1.0, t))
    return start + t * segment


def _segments_share_endpoint(first: np.ndarray, second: np.ndarray) -> bool:
    return bool(set(first.tolist()) & set(second.tolist()))


def segments_cross(a0: np.ndarray, a1: np.ndarray, b0: np.ndarray, b1: np.ndarray) -> bool:
    def orientation(p: np.ndarray, q: np.ndarray, r: np.ndarray) -> float:
        return float((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))

    o1 = orientation(a0, a1, b0)
    o2 = orientation(a0, a1, b1)
    o3 = orientation(b0, b1, a0)
    o4 = orientation(b0, b1, a1)
    return (o1 * o2 < 0.0) and (o3 * o4 < 0.0)


def summarize_edge_quality(
    row: Any,
    *,
    edge_quality_columns: Sequence[str],
    classify_quality_status: Callable[[Any, str], str],
) -> str:
    statuses = [
        classify_quality_status(getattr(row, column, np.nan), column)
        for column in edge_quality_columns
    ]
    non_missing = [status for status in statuses if status != "missing"]
    if not non_missing:
        return "missing"
    return max(non_missing, key=lambda status: QUALITY_STATUS_ORDER[status])


def format_ddg_label(mean_value: Any, std_value: Any) -> str:
    if pd.isna(mean_value):
        return "ΔΔG NA"
    if pd.isna(std_value):
        return f"ΔΔG {mean_value:+.2f}"
    return f"ΔΔG {mean_value:+.2f}\n± {std_value:.2f}"


def format_ligand_display_name(ligand_name: str) -> str:
    if len(ligand_name) <= 14:
        return ligand_name

    parts = ligand_name.split("_")
    if len(parts) == 1:
        midpoint = len(ligand_name) // 2
        return f"{ligand_name[:midpoint]}\n{ligand_name[midpoint:]}"

    best_split = 1
    best_score = float("inf")
    for split_index in range(1, len(parts)):
        first = " ".join(parts[:split_index])
        second = " ".join(parts[split_index:])
        score = abs(len(first) - len(second))
        if score < best_score:
            best_score = score
            best_split = split_index

    return (
        f"{' '.join(parts[:best_split])}\n"
        f"{' '.join(parts[best_split:])}"
    )


def default_classify_quality_status(value: Any, column: str) -> str:
    if pd.isna(value):
        return "missing"

    rule = NETWORK_QUALITY_RULES.get(column)
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


def load_transformation_summary_from_workup_csvs(workup_dir: str | Path) -> pd.DataFrame:
    workup_path = Path(workup_dir)
    energies_path = workup_path / "summary_energies.csv"
    convergence_path = workup_path / "summary_convergence.csv"

    if not energies_path.exists():
        raise FileNotFoundError(f"Could not find summary energies CSV: {energies_path}")
    if not convergence_path.exists():
        raise FileNotFoundError(
            f"Could not find summary convergence CSV: {convergence_path}"
        )

    energies = pd.read_csv(energies_path)
    convergence = pd.read_csv(convergence_path)
    return energies.merge(convergence, on=["ligand_i", "ligand_j"], how="inner")


def load_ligand_components(
    results_root: Path | None,
    transformation_summary: pd.DataFrame,
) -> dict[str, dict[str, Any]]:
    if results_root is None:
        return {}

    ligand_components: dict[str, dict[str, Any]] = {}
    for row in transformation_summary[["ligand_i", "ligand_j"]].itertuples(index=False):
        ligand_i = cast(str, row.ligand_i)
        ligand_j = cast(str, row.ligand_j)
        transformation_path = find_transformation_definition_path(
            results_root,
            ligand_i,
            ligand_j,
        )
        if transformation_path is None:
            continue

        for ligand_name, component in load_ligand_components_from_transformation(
            str(transformation_path)
        ).items():
            ligand_components.setdefault(ligand_name, component)

    return ligand_components


def find_transformation_definition_path(
    results_root: Path,
    ligand_i: str,
    ligand_j: str,
) -> Path | None:
    transformations_dir = results_root.parent / "network_setup" / "transformations"
    if not transformations_dir.exists():
        return None

    explicit_candidates = [
        transformations_dir / f"rbfe_{ligand_i}_complex_{ligand_j}_complex.json",
        transformations_dir / f"rbfe_{ligand_i}_solvent_{ligand_j}_solvent.json",
    ]
    for candidate in explicit_candidates:
        if candidate.exists():
            return candidate

    fallback_matches = sorted(
        transformations_dir.glob(f"rbfe_{ligand_i}_*_{ligand_j}_*.json")
    )
    if fallback_matches:
        return fallback_matches[0]

    return None


@lru_cache(maxsize=None)
def load_transformation_lookup(transformation_path: str) -> dict[str, Any]:
    return {key: value for key, value in json.loads(Path(transformation_path).read_text())}


def load_ligand_components_from_transformation(
    transformation_path: str,
) -> dict[str, dict[str, Any]]:
    lookup = load_transformation_lookup(transformation_path)
    ligand_components: dict[str, dict[str, Any]] = {}

    for key, value in lookup.items():
        if not key.startswith("SmallMoleculeComponent-"):
            continue

        ligand_name = value.get("molprops", {}).get("ofe-name")
        if ligand_name:
            ligand_components[ligand_name] = value

    return ligand_components


def build_aligned_ligand_node_images(
    ligand_components: Mapping[str, dict[str, Any]],
) -> dict[str, np.ndarray | None]:
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor, rdFMCS
    except ImportError:
        return {
            ligand_name: build_ligand_node_image(component)
            for ligand_name, component in ligand_components.items()
        }

    molecules: dict[str, Any] = {}
    for ligand_name, component in ligand_components.items():
        molecule = build_rdkit_molecule(component)
        if molecule is not None:
            molecules[ligand_name] = Chem.Mol(molecule)

    if not molecules:
        return {}

    rdDepictor.SetPreferCoordGen(True)
    reference_name = choose_reference_ligand_name(molecules)
    reference = Chem.Mol(molecules[reference_name])
    rdDepictor.Compute2DCoords(reference)

    aligned_images: dict[str, np.ndarray | None] = {}
    for ligand_name, molecule in molecules.items():
        depicted = Chem.Mol(molecule)
        if ligand_name != reference_name:
            align_ligand_depiction_to_reference(depicted, reference, rdFMCS, rdDepictor)
        aligned_images[ligand_name] = render_ligand_node_image(depicted)

    return aligned_images


def choose_reference_ligand_name(molecules: Mapping[str, Any]) -> str:
    try:
        from rdkit.Chem import rdFMCS
    except ImportError:
        return next(iter(molecules))

    scores: dict[str, int] = {}
    for candidate_name, candidate in molecules.items():
        score = 0
        for other_name, other in molecules.items():
            if candidate_name == other_name:
                continue
            mcs_result = rdFMCS.FindMCS(
                [candidate, other],
                atomCompare=rdFMCS.AtomCompare.CompareElements,
                bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                ringMatchesRingOnly=True,
                completeRingsOnly=True,
                timeout=2,
            )
            if mcs_result.numAtoms > 0:
                score += int(mcs_result.numAtoms)
        scores[candidate_name] = score

    return max(scores.items(), key=lambda item: item[1])[0]


def align_ligand_depiction_to_reference(
    molecule: Any,
    reference: Any,
    rdFMCS: Any,
    rdDepictor: Any,
) -> None:
    mcs_result = rdFMCS.FindMCS(
        [reference, molecule],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        timeout=2,
    )
    if mcs_result.numAtoms < 4:
        rdDepictor.Compute2DCoords(molecule)
        return

    try:
        from rdkit import Chem
    except ImportError:
        rdDepictor.Compute2DCoords(molecule)
        return

    query = Chem.MolFromSmarts(mcs_result.smartsString)
    if query is None:
        rdDepictor.Compute2DCoords(molecule)
        return

    ref_match = reference.GetSubstructMatch(query)
    mol_match = molecule.GetSubstructMatch(query)
    if not ref_match or not mol_match:
        rdDepictor.Compute2DCoords(molecule)
        return

    atom_map = list(zip(ref_match, mol_match))
    try:
        rdDepictor.GenerateDepictionMatching2DStructure(
            molecule,
            reference,
            atom_map,
        )
    except Exception:
        rdDepictor.Compute2DCoords(molecule)


def build_ligand_node_image(component: dict[str, Any]) -> np.ndarray | None:
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
    except ImportError:
        return None

    molecule = build_rdkit_molecule(component)
    if molecule is None:
        return None

    molecule = Chem.Mol(molecule)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(molecule)

    return render_ligand_node_image(molecule)


def render_ligand_node_image(molecule: Any) -> np.ndarray | None:
    try:
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        return None

    drawer = rdMolDraw2D.MolDraw2DCairo(320, 240)
    drawer_options = drawer.drawOptions()
    drawer_options.padding = 0.04
    if hasattr(drawer_options, "clearBackground"):
        drawer_options.clearBackground = False
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, molecule)
    drawer.FinishDrawing()

    return plt.imread(io.BytesIO(drawer.GetDrawingText()), format="png")


def build_rdkit_molecule(component: dict[str, Any]) -> Any:
    try:
        from rdkit import Chem
    except ImportError:
        return None

    atoms = component.get("atoms", [])
    bonds = component.get("bonds", [])
    if not atoms:
        return None

    molecule = Chem.RWMol()
    for atom_data in atoms:
        atom = Chem.Atom(int(atom_data[0]))
        atom.SetFormalCharge(int(atom_data[2]))
        if int(atom_data[1]):
            atom.SetIsotope(int(atom_data[1]))
        if bool(atom_data[3]):
            atom.SetIsAromatic(True)
        molecule.AddAtom(atom)

    bond_type_map = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        12: Chem.BondType.AROMATIC,
    }
    for bond_data in bonds:
        begin_atom, end_atom, bond_order = bond_data[:3]
        bond_type = bond_type_map.get(int(bond_order), Chem.BondType.SINGLE)
        molecule.AddBond(int(begin_atom), int(end_atom), bond_type)
        if bond_type == Chem.BondType.AROMATIC:
            bond = molecule.GetBondBetweenAtoms(int(begin_atom), int(end_atom))
            if bond is not None:
                bond.SetIsAromatic(True)

    molecule = molecule.GetMol()
    sanitize_result = Chem.SanitizeMol(molecule, catchErrors=True)
    if sanitize_result != Chem.SanitizeFlags.SANITIZE_NONE:
        return None

    try:
        return Chem.RemoveHs(molecule)
    except Exception:
        return molecule


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Render summary_ligand_network.png from workup CSV exports without relying "
            "on workup.py."
        )
    )
    parser.add_argument(
        "workup_dir",
        nargs="?",
        default="example/workup",
        help="Directory containing summary_energies.csv and summary_convergence.csv.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help=(
            "Results root used to locate transformation definitions for ligand depictions. "
            "Defaults to <workup_dir>/../results when present."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output PNG path. Defaults to <workup_dir>/summary_ligand_network.png."
        ),
    )
    return parser


def main() -> None:
    parser = build_argument_parser()
    args = parser.parse_args()

    workup_dir = Path(args.workup_dir).resolve()
    transformation_summary = load_transformation_summary_from_workup_csvs(workup_dir)

    results_dir = args.results_dir
    if results_dir is None:
        default_results_dir = workup_dir.parent / "results"
        results_dir = default_results_dir if default_results_dir.exists() else None

    output_path = args.output or (workup_dir / "summary_ligand_network.png")
    plot_summary_ligand_network(
        transformation_summary,
        output_path,
        results_root=results_dir.resolve() if results_dir is not None else None,
    )
    print(f"Wrote ligand network plot to {output_path}")


if __name__ == "__main__":
    main()
