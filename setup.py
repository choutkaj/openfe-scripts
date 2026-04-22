from pathlib import Path

from setuptools import setup


README_PATH = Path(__file__).with_name("README.md")


setup(
    name="openfe-scripts",
    version="0.1.0",
    description="Helper scripts that extend the OpenFE CLI workflow.",
    long_description=README_PATH.read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    python_requires=">=3.10",
    py_modules=["plot_network"],
    scripts=[
        "prep-rbfe-hybridtop.py",
        "prep-rbfe-septop.py",
        "workup-hybridtop.py",
        "workup-septop.py",
        "plot-network.py",
    ],
)
