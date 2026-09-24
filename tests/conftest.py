from __future__ import annotations

from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "tests" / "fixtures"
DEMO_FILES = [
    ROOT / "ref_1kgPhase3eur_LDr2p1_DEMO.mat",
    ROOT / "CTG_COG_2018_DEMO.mat",
    ROOT / "SSGAC_EDU_2016_DEMO.mat",
]

requires_demo = pytest.mark.skipif(
    not all(p.exists() for p in DEMO_FILES),
    reason="demo data missing: tar -xzvf pleioFDR_demo_data.tar.gz in the repository root",
)


@pytest.fixture(scope="session")
def demo_reference():
    from pleiofdr.io import load_reference

    if not DEMO_FILES[0].exists():
        pytest.skip("demo data missing")
    return load_reference(DEMO_FILES[0])
