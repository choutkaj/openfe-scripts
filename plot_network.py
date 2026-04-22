#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
from pathlib import Path


_SCRIPT_PATH = Path(__file__).with_name("plot-network.py")
_SPEC = importlib.util.spec_from_file_location("plot_network_cli", _SCRIPT_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError(f"Could not load plot helper from {_SCRIPT_PATH}")

_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

plot_summary_ligand_network = _MODULE.plot_summary_ligand_network
main = _MODULE.main


if __name__ == "__main__":
    main()
