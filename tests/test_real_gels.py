"""Regression test using the user-supplied 50 bp Mycoplasma ARDRA gel."""

from pathlib import Path

import numpy as np
import pytest

from ardronics.gel import LADDERS, analyze_gel

GEL = Path(__file__).parents[1] / "data" / "ardra_example_1.jpg"
NAMES = [
    "M_gallisepticum", "M_pneumoniae", "M_glycophilum", "M_bovirhinis",
    "M_columborale", "M_neurolyticum", "M_flocculare", "M_dispar",
    "M_hyopneumoniae", "M_hyorhinis", "M_arginini", "M_hyosynoviae",
    "M_orale", "M_columbinasale", "M_agalactiae", "M_bovis",
    "M_lipofaciens", "M_gallinarum",
]


@pytest.mark.skipif(not GEL.exists(), reason="user-supplied data/ardra_example_1.jpg is unavailable")
def test_mycoplasma_gel_with_two_ladders(tmp_path: Path):
    result = analyze_gel(
        GEL, "unknown", tmp_path, n_lanes=20, crop=(20, 138, 535, 334),
        ladder_lane=1, secondary_ladder_lane=20, ladder_sizes=LADDERS["50bp"],
        sample_names=NAMES,
    )
    assert set(result.bands["sample"]) == set(NAMES) | {"ladder"}
    assert min(result.lane_centers) >= 0
    assert max(result.lane_centers) < 535
    for lane in (1, 20):
        ladder = result.bands.loc[result.bands["lane"] == lane]
        assert len(ladder) == 13
        assert np.allclose(ladder["fragment_bp"], LADDERS["50bp"], rtol=1e-6)
    assert all((result.bands["sample"] == name).any() for name in NAMES)
    assert (tmp_path / "qc_overlay.png").stat().st_size > 1000


ALUI_GEL = Path(__file__).parents[1] / "data" / "Vol_8_no2_Typi_Afaf_fig5.jpg"
ALUI_Y = [59.8245, 71.3201, 82.2045, 87.8058, 94.2063, 101.1093, 109.5424, 119.8382, 132.6827, 146.0004, 163.7537, 186.1212, 216.2830]
ALUI_SIZES = [2072, 1500, 1000, 900, 800, 700, 600, 500, 400, 300, 250, 200, 100]


@pytest.mark.skipif(not ALUI_GEL.exists(), reason="user-supplied AluI gel is unavailable")
def test_alui_gel_with_reviewed_ladder(tmp_path: Path):
    names = ["41", "42", "44", "45", "46", "47", "48", "51", "52", "53", "54", "55", "MW_right_uncalibrated"]
    result = analyze_gel(
        ALUI_GEL, "AluI", tmp_path, n_lanes=14, crop=(88, 42, 339, 285),
        ladder_lane=1, ladder_sizes=ALUI_SIZES, ladder_y=ALUI_Y, sample_names=names,
    )
    ladder = result.bands.loc[result.bands["lane"] == 1]
    assert len(ladder) == 13
    assert np.allclose(ladder["fragment_bp"], ALUI_SIZES, rtol=1e-4)
    assert set(result.bands["sample"]) == set(names) | {"ladder"}
    sample_41 = result.bands.loc[result.bands["sample"] == "41", "fragment_bp"]
    assert len(sample_41) == 3
    assert sample_41.between(100, 2072).all()
