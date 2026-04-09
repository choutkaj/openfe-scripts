#!/usr/bin/env python3

from __future__ import annotations

import argparse
from collections import Counter
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

LOGGER = logging.getLogger("prep_rbfe_hybridtop")


@dataclass(frozen=True)
class CliConfig:
    receptor_path: Path
    ligands_path: Path
    mapper: str
    scorer: str
    network: str
    custom_network_path: Path | None
    central_ligand: str | None
    n_windows: int
    window_length_ns: float
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

    # Suppress external library loggers that may attach their own handlers.
    for logger_name in (
        "openfe",
        "gufe",
        "kartograf",
        "lomap",
        "rdkit",
        "matplotlib",
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
            "Prepare an OpenFE hybrid-topology RBFE network and write the "
            "transformation JSON files to disk."
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
        default=11,
        help="Number of lambda windows and OpenMM replicas to use.",
    )
    parser.add_argument(
        "--windowtime",
        type=positive_float,
        default=5.0,
        help="Production length per window in nanoseconds.",
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

    LOGGER.info("Loaded %d ligands", len(ligands))
    return ligands


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

    ligand_name_counts = Counter(ligand.name for ligand in ligands)
    duplicate_ligand_names = sorted(
        ligand_name
        for ligand_name, count in ligand_name_counts.items()
        if count > 1
    )
    if duplicate_ligand_names:
        duplicate_names_text = ", ".join(duplicate_ligand_names)
        raise ValueError(
            "Custom network mode requires unique ligand names in the SDF. "
            f"Duplicate ligand names found: {duplicate_names_text}"
        )

    valid_ligand_names = set(ligand_name_counts)
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


def build_protocol(
    *,
    window_length_ns: float,
    n_windows: int,
    solvent_padding_nm: float,
) -> Any:
    from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
    from openff.units import unit

    settings: Any = RelativeHybridTopologyProtocol.default_settings()
    settings.solvation_settings.solvent_padding = solvent_padding_nm * unit.nanometer
    settings.simulation_settings.production_length = window_length_ns * unit.nanosecond
    settings.simulation_settings.n_replicas = n_windows
    settings.lambda_settings.lambda_windows = n_windows
    settings.protocol_repeats = 1
    return RelativeHybridTopologyProtocol(settings)


def create_alchemical_network(
    ligand_network: Any,
    receptor_path: Path,
    mappings_dir: Path,
    *,
    window_length_ns: float,
    n_windows: int,
) -> Any:
    import openfe

    LOGGER.info("Loading receptor from %s", receptor_path)
    solvent = openfe.SolventComponent()
    protein = openfe.ProteinComponent.from_pdb_file(str(receptor_path))

    solvent_protocol = build_protocol(
        window_length_ns=window_length_ns,
        n_windows=n_windows,
        solvent_padding_nm=2.0,
    )
    complex_protocol = build_protocol(
        window_length_ns=window_length_ns,
        n_windows=n_windows,
        solvent_padding_nm=1.0,
    )

    LOGGER.info("Creating alchemical network")
    transformations: list[Any] = []
    for mapping in ligand_network.edges:
        LOGGER.info("%s -> %s", mapping.componentA.name, mapping.componentB.name)
        LOGGER.info("Score: %s", mapping.annotations)
        mapping.draw_to_file(
            str(mappings_dir / f"mapping_{mapping.componentA.name}_{mapping.componentB.name}.png")
        )

        for leg, protocol in (
            ("solvent", solvent_protocol),
            ("complex", complex_protocol),
        ):
            system_a = {
                "ligand": mapping.componentA,
                "solvent": solvent,
            }
            system_b = {
                "ligand": mapping.componentB,
                "solvent": solvent,
            }
            if leg == "complex":
                system_a["protein"] = protein
                system_b["protein"] = protein

            state_a = openfe.ChemicalSystem(
                system_a,
                name=f"{mapping.componentA.name}_{leg}",
            )
            state_b = openfe.ChemicalSystem(
                system_b,
                name=f"{mapping.componentB.name}_{leg}",
            )
            transformations.append(
                openfe.Transformation(
                    stateA=state_a,
                    stateB=state_b,
                    mapping=mapping,
                    protocol=protocol,
                    name=f"rbfe_{state_a.name}_{state_b.name}",
                )
            )

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
        mappings_dir, transformation_dir = write_ligand_network_artifacts(
            ligand_network,
            config.output_dir,
        )
        alchemical_network = create_alchemical_network(
            ligand_network,
            config.receptor_path,
            mappings_dir,
            window_length_ns=config.window_length_ns,
            n_windows=config.n_windows,
        )
        written_transformations = write_transformations(
            alchemical_network,
            transformation_dir,
        )
    except ImportError as exc:
        LOGGER.error(
            "Missing dependency: %s. Activate the micromamba environment that contains "
            "OpenFE, RDKit, OpenFF, and Kartograf before running this script.",
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
