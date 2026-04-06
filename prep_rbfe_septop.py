#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import Counter
import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

LOGGER = logging.getLogger("prep_rbfe_septop")
LAMBDA_SETTING_FIELDS = (
    "lambda_elec_A",
    "lambda_elec_B",
    "lambda_vdw_A",
    "lambda_vdw_B",
    "lambda_restraints_A",
    "lambda_restraints_B",
)


@dataclass(frozen=True)
class CliConfig:
    receptor_path: Path
    ligands_path: Path
    mapper: str
    scorer: str
    network: str
    custom_network_path: Path | None
    central_ligand: str | None
    n_windows: int | None
    window_length_ns: float | None
    equilibration_length_ns: float | None
    protocol_repeats: int
    host_min_distance_nm: float
    host_max_distance_nm: float
    output_dir: Path


def configure_logging() -> None:
    LOGGER.handlers.clear()
    LOGGER.setLevel(logging.INFO)
    LOGGER.propagate = False

    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("%(message)s"))
    LOGGER.addHandler(handler)

    root_logger = logging.getLogger()
    root_logger.handlers.clear()
    root_logger.addHandler(logging.NullHandler())
    root_logger.setLevel(logging.CRITICAL + 1)

    for logger_name in (
        "openfe",
        "gufe",
        "kartograf",
        "lomap",
        "rdkit",
        "matplotlib",
        "openmm",
        "openmmtools",
        "pymbar",
    ):
        external_logger = logging.getLogger(logger_name)
        external_logger.handlers.clear()
        external_logger.propagate = False
        external_logger.disabled = True


def existing_file(path_value: str) -> Path:
    path = Path(path_value)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"File does not exist: {path}")
    return path


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("Value must be a positive integer.")
    return parsed


def positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("Value must be a positive number.")
    return parsed


def parse_args(argv: Sequence[str] | None = None) -> CliConfig:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare an OpenFE SepTop RBFE network and write the transformation "
            "JSON files to disk."
        ),
    )
    parser.add_argument(
        "--rec",
        type=existing_file,
        required=True,
        help="Path to the prepared receptor PDB file.",
    )
    parser.add_argument(
        "--ligs",
        type=existing_file,
        required=True,
        help="Path to the ligand SDF file.",
    )
    parser.add_argument(
        "--mapper",
        choices=("lomap", "kartograf"),
        default="kartograf",
        help="Atom mapper used to build the ligand network.",
    )
    parser.add_argument(
        "--scorer",
        choices=("lomap", "kartograf_rmsd"),
        default="lomap",
        help="Mapping scorer used during network planning.",
    )
    parser.add_argument(
        "--network",
        choices=("minimal_spanning", "minimal_redundant", "radial", "maximal", "custom"),
        default="minimal_spanning",
        help="Ligand network planner to use.",
    )
    parser.add_argument(
        "--custom-network",
        type=existing_file,
        default=None,
        help=(
            "Path to a YAML file defining custom ligand edges under a top-level "
            "'edges' key. Only used with --network custom."
        ),
    )
    parser.add_argument(
        "--central-ligand",
        default=None,
        help=(
            "Name of the central ligand for radial networks. If omitted, the first "
            "ligand in the SDF is used."
        ),
    )
    parser.add_argument(
        "--windows",
        type=positive_int,
        default=None,
        help=(
            "Number of lambda windows and OpenMM replicas to use for both solvent "
            "and complex SepTop legs. If omitted, the protocol default schedule is used."
        ),
    )
    parser.add_argument(
        "--windowtime",
        type=positive_float,
        default=None,
        help=(
            "Production length per lambda window in nanoseconds for both solvent "
            "and complex SepTop legs. If omitted, the protocol default is used."
        ),
    )
    parser.add_argument(
        "--equilibration-time",
        type=positive_float,
        default=2.0,
        help=(
            "Equilibration length in nanoseconds for both solvent and complex "
            "alchemical simulations. Default: 2.0 ns."
        ),
    )
    parser.add_argument(
        "--protocol-repeats",
        type=positive_int,
        default=1,
        help=(
            "Number of protocol repeats encoded into each SepTop transformation. "
            "Default: 1."
        ),
    )
    parser.add_argument(
        "--host-min-distance",
        type=positive_float,
        default=0.5,
        help=(
            "Minimum distance in nanometers between protein and ligand atoms when "
            "building complex restraints. Default: 0.5 nm."
        ),
    )
    parser.add_argument(
        "--host-max-distance",
        type=positive_float,
        default=1.5,
        help=(
            "Maximum distance in nanometers between protein and ligand atoms when "
            "building complex restraints. Default: 1.5 nm."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory where network artifacts and transformations will be written.",
    )
    args = parser.parse_args(argv)

    if args.network == "custom" and args.custom_network is None:
        parser.error("--custom-network is required when --network custom is selected.")
    if args.network != "custom" and args.custom_network is not None:
        parser.error("--custom-network can only be used with --network custom.")
    if args.windows is not None and args.windows < 2:
        parser.error("--windows must be at least 2 for a valid lambda schedule.")
    if args.host_min_distance >= args.host_max_distance:
        parser.error("--host-min-distance must be smaller than --host-max-distance.")

    return CliConfig(
        receptor_path=args.rec,
        ligands_path=args.ligs,
        mapper=args.mapper,
        scorer=args.scorer,
        network=args.network,
        custom_network_path=args.custom_network,
        central_ligand=args.central_ligand,
        n_windows=args.windows,
        window_length_ns=args.windowtime,
        equilibration_length_ns=args.equilibration_time,
        protocol_repeats=args.protocol_repeats,
        host_min_distance_nm=args.host_min_distance,
        host_max_distance_nm=args.host_max_distance,
        output_dir=args.output_dir,
    )


def load_ligands(ligands_path: Path) -> list[Any]:
    import openfe
    from rdkit import Chem

    LOGGER.info("Loading ligands from %s", ligands_path)
    supplier = Chem.SDMolSupplier(str(ligands_path), removeHs=False)

    ligands: list[Any] = []
    invalid_records: list[int] = []
    for index, molecule in enumerate(supplier):
        if molecule is None:
            invalid_records.append(index)
            continue
        ligands.append(openfe.SmallMoleculeComponent.from_rdkit(molecule))

    if invalid_records:
        invalid_text = ", ".join(str(index) for index in invalid_records)
        raise ValueError(
            "Encountered invalid molecules in the SDF file at record indices: "
            f"{invalid_text}"
        )
    if len(ligands) < 2:
        raise ValueError("The ligand SDF must contain at least two valid molecules.")

    validate_ligand_names(ligands)
    LOGGER.info("Loaded %d ligands", len(ligands))
    return ligands


def validate_ligand_names(ligands: Sequence[Any]) -> None:
    missing_name_indices = [
        index
        for index, ligand in enumerate(ligands)
        if getattr(ligand, "name", None) is None
        or not str(getattr(ligand, "name")).strip()
    ]
    if missing_name_indices:
        missing_text = ", ".join(str(index) for index in missing_name_indices)
        raise ValueError(
            "All ligands must have non-empty names so the script can create stable "
            f"transformation names and output files. Missing names at SDF record indices: {missing_text}"
        )

    ligand_name_counts = Counter(str(ligand.name) for ligand in ligands)
    duplicate_ligand_names = sorted(
        ligand_name for ligand_name, count in ligand_name_counts.items() if count > 1
    )
    if duplicate_ligand_names:
        duplicate_names_text = ", ".join(duplicate_ligand_names)
        raise ValueError(
            "All ligands must have unique names so transformation filenames do not "
            f"collide. Duplicate ligand names found: {duplicate_names_text}"
        )


def assign_partial_charges(ligands: Sequence[Any]) -> list[Any]:
    from openfe.protocols.openmm_utils.charge_generation import (
        bulk_assign_partial_charges,
    )
    from openfe.protocols.openmm_utils.omm_settings import (
        OpenFFPartialChargeSettings,
    )

    LOGGER.info("Assigning partial charges to ligands")
    charge_settings = OpenFFPartialChargeSettings(
        partial_charge_method="am1bcc",
        off_toolkit_backend="ambertools",
    )
    return bulk_assign_partial_charges(
        molecules=list(ligands),
        overwrite=False,
        method=charge_settings.partial_charge_method,
        toolkit_backend=charge_settings.off_toolkit_backend,
        generate_n_conformers=charge_settings.number_of_conformers,
        nagl_model=charge_settings.nagl_model,
        processors=1,
    )


def load_custom_network_edges(
    custom_network_path: Path,
    ligands: Sequence[Any],
) -> list[tuple[str, str]]:
    import yaml

    LOGGER.info("Loading custom network definition from %s", custom_network_path)
    custom_network = yaml.safe_load(custom_network_path.read_text(encoding="utf-8"))

    if not isinstance(custom_network, dict) or "edges" not in custom_network:
        raise ValueError("Custom network YAML must contain a top-level 'edges' key.")

    raw_edges = custom_network["edges"]
    if not isinstance(raw_edges, list):
        raise ValueError(
            "Custom network YAML 'edges' must be a list of 2-item ligand-name pairs."
        )
    if not raw_edges:
        raise ValueError("Custom network YAML must contain at least one edge.")

    valid_ligand_names = {str(ligand.name) for ligand in ligands}
    custom_edges: list[tuple[str, str]] = []
    seen_edges: set[tuple[str, str]] = set()
    for index, raw_edge in enumerate(raw_edges):
        if not isinstance(raw_edge, (list, tuple)) or len(raw_edge) != 2:
            raise ValueError(
                f"Custom network edge at index {index} must be a 2-item ligand-name pair."
            )

        ligand_a, ligand_b = raw_edge
        if not isinstance(ligand_a, str) or not isinstance(ligand_b, str):
            raise ValueError(
                f"Custom network edge at index {index} must contain string ligand names."
            )
        if not ligand_a or not ligand_b:
            raise ValueError(
                f"Custom network edge at index {index} must contain non-empty ligand names."
            )
        if ligand_a == ligand_b:
            raise ValueError(
                f"Custom network edge at index {index} must connect two different ligands."
            )

        missing_names = [
            ligand_name
            for ligand_name in (ligand_a, ligand_b)
            if ligand_name not in valid_ligand_names
        ]
        if missing_names:
            missing_names_text = ", ".join(missing_names)
            raise ValueError(
                f"Custom network edge at index {index} references unknown ligand "
                f"name(s): {missing_names_text}"
            )

        if ligand_a < ligand_b:
            edge_key = (ligand_a, ligand_b)
        else:
            edge_key = (ligand_b, ligand_a)
        if edge_key in seen_edges:
            raise ValueError(
                f"Duplicate custom network edge at index {index}: {ligand_a}, {ligand_b}"
            )

        seen_edges.add(edge_key)
        custom_edges.append((ligand_a, ligand_b))

    LOGGER.info("Loaded %d custom network edges", len(custom_edges))
    return custom_edges


def build_mapper(mapper_name: str) -> Any:
    import openfe

    if mapper_name == "kartograf":
        from kartograf import KartografAtomMapper

        return KartografAtomMapper(atom_map_hydrogens=True)
    return openfe.LomapAtomMapper(max3d=1.0, element_change=False)


def build_scorer(scorer_name: str) -> Any:
    import openfe

    if scorer_name == "kartograf_rmsd":
        from kartograf.atom_mapping_scorer import MappingRMSDScorer

        return MappingRMSDScorer()
    return openfe.lomap_scorers.default_lomap_score


def create_ligand_network(
    charged_ligands: Sequence[Any],
    mapper_name: str,
    scorer_name: str,
    network_name: str,
    custom_network_path: Path | None,
    central_ligand_name: str | None,
) -> Any:
    import openfe

    charged_ligands = list(charged_ligands)
    mapper = build_mapper(mapper_name)
    LOGGER.info(
        "Planning a SepTop ligand network. Atom mappings are used only to choose edges "
        "and generate diagnostic plots; the SepTop transformations themselves do not "
        "store atom mappings."
    )

    if network_name != "radial" and central_ligand_name is not None:
        LOGGER.info(
            "Ignoring --central-ligand because it is only used when --network radial is selected"
        )

    if network_name == "custom":
        if custom_network_path is None:
            raise ValueError(
                "--custom-network is required when --network custom is selected."
            )
        LOGGER.info("Ignoring --scorer because it is not used when --network custom is selected")
        custom_edges = load_custom_network_edges(custom_network_path, charged_ligands)
        LOGGER.info(
            "Creating ligand network with network=%s mapper=%s",
            network_name,
            mapper_name,
        )
        return openfe.ligand_network_planning.generate_network_from_names(
            ligands=charged_ligands,
            mapper=mapper,
            names=custom_edges,
        )

    scorer = build_scorer(scorer_name)
    network_planner = {
        "minimal_spanning": openfe.ligand_network_planning.generate_minimal_spanning_network,
        "minimal_redundant": openfe.ligand_network_planning.generate_minimal_redundant_network,
        "radial": openfe.ligand_network_planning.generate_radial_network,
        "maximal": openfe.ligand_network_planning.generate_maximal_network,
    }[network_name]

    planner_kwargs: dict[str, Any] = {
        "ligands": charged_ligands,
        "mappers": [mapper],
        "scorer": scorer,
    }
    if network_name == "radial":
        if central_ligand_name is None:
            central_ligand = charged_ligands[0]
            LOGGER.info(
                "Using radial network with central ligand %s (first ligand in input)",
                central_ligand.name,
            )
        else:
            matches = [
                ligand for ligand in charged_ligands if ligand.name == central_ligand_name
            ]
            if not matches:
                available_names = ", ".join(ligand.name for ligand in charged_ligands)
                raise ValueError(
                    f"central-ligand '{central_ligand_name}' was not found in the input SDF. "
                    f"Available ligand names: {available_names}"
                )
            if len(matches) > 1:
                raise ValueError(
                    f"central-ligand '{central_ligand_name}' is ambiguous because multiple "
                    "ligands share that name."
                )
            central_ligand = matches[0]
            LOGGER.info(
                "Using radial network with central ligand %s (selected by --central-ligand)",
                central_ligand.name,
            )
        planner_kwargs["central_ligand"] = central_ligand

    LOGGER.info(
        "Creating ligand network with network=%s mapper=%s scorer=%s",
        network_name,
        mapper_name,
        scorer_name,
    )
    return network_planner(**planner_kwargs)


def write_ligand_network_artifacts(ligand_network: Any, output_dir: Path) -> tuple[Path, Path]:
    from openfe.utils.atommapping_network_plotting import plot_atommapping_network

    output_dir.mkdir(parents=True, exist_ok=True)
    network_setup_dir = output_dir / "network_setup"
    network_setup_dir.mkdir(parents=True, exist_ok=True)
    mappings_dir = output_dir / "mappings"
    mappings_dir.mkdir(parents=True, exist_ok=True)
    transformation_dir = network_setup_dir / "transformations"
    transformation_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info("Writing ligand network plot")
    network_figure = plot_atommapping_network(ligand_network)
    network_figure.savefig(
        network_setup_dir / "ligand_network.png",
        dpi=300,
        bbox_inches="tight",
    )

    LOGGER.info("Writing ligand network graph")
    (network_setup_dir / "ligand_network.graphml").write_text(
        ligand_network.to_graphml(),
        encoding="utf-8",
    )

    return mappings_dir, transformation_dir


def get_lambda_schedule_length(lambda_settings: Any, settings_name: str) -> int:
    schedule_lengths = {
        field_name: len(getattr(lambda_settings, field_name))
        for field_name in LAMBDA_SETTING_FIELDS
    }
    unique_lengths = set(schedule_lengths.values())
    if len(unique_lengths) != 1:
        details = ", ".join(
            f"{field_name}={field_length}"
            for field_name, field_length in sorted(schedule_lengths.items())
        )
        raise ValueError(
            f"The {settings_name} lambda settings are internally inconsistent: {details}"
        )
    return next(iter(unique_lengths))


def resample_lambda_schedule(values: Sequence[float], n_windows: int) -> list[float]:
    if n_windows < 2:
        raise ValueError("SepTop lambda schedules require at least two windows.")

    original_values = [float(value) for value in values]
    if len(original_values) < 2:
        raise ValueError(
            "Cannot resample a lambda schedule with fewer than two original points."
        )
    if len(original_values) == n_windows:
        return original_values

    resampled: list[float] = []
    scale = (len(original_values) - 1) / (n_windows - 1)
    for output_index in range(n_windows):
        source_position = output_index * scale
        lower_index = int(math.floor(source_position))
        upper_index = int(math.ceil(source_position))
        if lower_index == upper_index:
            resampled.append(original_values[lower_index])
            continue

        upper_weight = source_position - lower_index
        lower_weight = 1.0 - upper_weight
        resampled.append(
            (lower_weight * original_values[lower_index])
            + (upper_weight * original_values[upper_index])
        )

    resampled[0] = original_values[0]
    resampled[-1] = original_values[-1]
    return resampled


def resize_lambda_settings(lambda_settings: Any, settings_name: str, n_windows: int) -> None:
    current_windows = get_lambda_schedule_length(lambda_settings, settings_name)
    if current_windows == n_windows:
        LOGGER.info(
            "Using %s lambda schedule with %d windows",
            settings_name,
            current_windows,
        )
        return

    for field_name in LAMBDA_SETTING_FIELDS:
        values = getattr(lambda_settings, field_name)
        setattr(
            lambda_settings,
            field_name,
            resample_lambda_schedule(values, n_windows),
        )

    LOGGER.info(
        "Resampled %s lambda schedule from %d to %d windows",
        settings_name,
        current_windows,
        n_windows,
    )


def build_protocol(
    *,
    n_windows: int | None,
    window_length_ns: float | None,
    equilibration_length_ns: float | None,
    protocol_repeats: int,
    host_min_distance_nm: float,
    host_max_distance_nm: float,
) -> Any:
    from openfe.protocols.openmm_septop import SepTopProtocol
    from openff.units import unit

    settings: Any = SepTopProtocol.default_settings()
    settings.protocol_repeats = protocol_repeats
    settings.complex_restraint_settings.host_min_distance = (
        host_min_distance_nm * unit.nanometer
    )
    settings.complex_restraint_settings.host_max_distance = (
        host_max_distance_nm * unit.nanometer
    )

    if window_length_ns is not None:
        settings.solvent_simulation_settings.production_length = (
            window_length_ns * unit.nanosecond
        )
        settings.complex_simulation_settings.production_length = (
            window_length_ns * unit.nanosecond
        )

    if equilibration_length_ns is not None:
        settings.solvent_simulation_settings.equilibration_length = (
            equilibration_length_ns * unit.nanosecond
        )
        settings.complex_simulation_settings.equilibration_length = (
            equilibration_length_ns * unit.nanosecond
        )

    solvent_default_windows = get_lambda_schedule_length(
        settings.solvent_lambda_settings,
        "solvent",
    )
    complex_default_windows = get_lambda_schedule_length(
        settings.complex_lambda_settings,
        "complex",
    )
    solvent_active_windows = solvent_default_windows
    complex_active_windows = complex_default_windows
    if n_windows is None:
        LOGGER.info(
            "Using SepTop default lambda schedules: solvent=%d windows, complex=%d windows",
            solvent_active_windows,
            complex_active_windows,
        )
    else:
        resize_lambda_settings(settings.solvent_lambda_settings, "solvent", n_windows)
        resize_lambda_settings(settings.complex_lambda_settings, "complex", n_windows)
        solvent_active_windows = n_windows
        complex_active_windows = n_windows

    settings.solvent_simulation_settings.n_replicas = solvent_active_windows
    settings.complex_simulation_settings.n_replicas = complex_active_windows

    if solvent_active_windows == complex_active_windows:
        LOGGER.info(
            "Configured SepTop protocol with repeats=%d and %d lambda windows per leg",
            protocol_repeats,
            solvent_active_windows,
        )
    else:
        LOGGER.info(
            "Configured SepTop protocol with repeats=%d, solvent windows=%d, and complex windows=%d",
            protocol_repeats,
            solvent_active_windows,
            complex_active_windows,
        )
    return SepTopProtocol(settings)


def create_alchemical_network(
    ligand_network: Any,
    receptor_path: Path,
    mappings_dir: Path,
    protocol: Any,
) -> Any:
    import openfe

    LOGGER.info("Loading receptor from %s", receptor_path)
    solvent = openfe.SolventComponent()
    protein = openfe.ProteinComponent.from_pdb_file(str(receptor_path))

    LOGGER.info("Creating SepTop alchemical network")
    transformations: list[Any] = []
    for mapping in ligand_network.edges:
        LOGGER.info("%s -> %s", mapping.componentA.name, mapping.componentB.name)
        LOGGER.info("Score: %s", mapping.annotations)
        mapping.draw_to_file(
            str(mappings_dir / f"mapping_{mapping.componentA.name}_{mapping.componentB.name}.png")
        )

        state_a = openfe.ChemicalSystem(
            {
                "ligand": mapping.componentA,
                "protein": protein,
                "solvent": solvent,
            },
            name=f"{mapping.componentA.name}",
        )
        state_b = openfe.ChemicalSystem(
            {
                "ligand": mapping.componentB,
                "protein": protein,
                "solvent": solvent,
            },
            name=f"{mapping.componentB.name}",
        )
        transformations.append(
            openfe.Transformation(
                stateA=state_a,
                stateB=state_b,
                mapping=None,
                protocol=protocol,
                name=f"rbfe_{state_a.name}_{state_b.name}",
            )
        )

    if not transformations:
        raise ValueError(
            "The ligand network did not contain any edges, so no SepTop transformations "
            "could be created."
        )

    LOGGER.info("Created %d SepTop transformations", len(transformations))
    return openfe.AlchemicalNetwork(transformations)


def write_transformations(network: Any, transformation_dir: Path) -> int:
    LOGGER.info("Writing transformations to %s", transformation_dir)

    written = 0
    for transformation in network.edges:
        transformation.to_json(transformation_dir / f"{transformation.name}.json")
        written += 1
    return written


def main(argv: Sequence[str] | None = None) -> int:
    configure_logging()
    config = parse_args(argv)

    try:
        ligands = load_ligands(config.ligands_path)
        charged_ligands = assign_partial_charges(ligands)
        ligand_network = create_ligand_network(
            charged_ligands,
            mapper_name=config.mapper,
            scorer_name=config.scorer,
            network_name=config.network,
            custom_network_path=config.custom_network_path,
            central_ligand_name=config.central_ligand,
        )
        protocol = build_protocol(
            n_windows=config.n_windows,
            window_length_ns=config.window_length_ns,
            equilibration_length_ns=config.equilibration_length_ns,
            protocol_repeats=config.protocol_repeats,
            host_min_distance_nm=config.host_min_distance_nm,
            host_max_distance_nm=config.host_max_distance_nm,
        )
        mappings_dir, transformation_dir = write_ligand_network_artifacts(
            ligand_network,
            config.output_dir,
        )
        alchemical_network = create_alchemical_network(
            ligand_network,
            config.receptor_path,
            mappings_dir,
            protocol,
        )
        written_transformations = write_transformations(
            alchemical_network,
            transformation_dir,
        )
    except ImportError as exc:
        LOGGER.error(
            "Missing dependency: %s. Activate the micromamba environment that contains "
            "OpenFE 1.7.0 or newer with SepTop support, RDKit, OpenFF, and Kartograf "
            "before running this script.",
            exc,
        )
        return 1
    except ValueError as exc:
        LOGGER.error(str(exc))
        return 1

    LOGGER.info("Finished. Wrote %d transformation files.", written_transformations)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
