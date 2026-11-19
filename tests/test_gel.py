from pathlib import Path

import numpy as np
import pytest
from conftest import LADDER_SIZES, SAMPLE_SIZES, migration_y

from ardronics.gel import analyze_gel, calibrate_ladder


def test_synthetic_gel_sizes_and_intensities(synthetic_gel: Path, tmp_path: Path):
    result = analyze_gel(
        synthetic_gel, "HaeIII", tmp_path / "out", n_lanes=4,
        ladder_lane=1, ladder_sizes=LADDER_SIZES,
        sample_names=["a", "b", "c"],
    )
    assert len(result.lane_centers) == 4
    assert np.max(np.abs(np.array(result.lane_centers) - [120, 340, 560, 780])) <= 10
    samples = result.bands.loc[result.bands["sample"] != "ladder"]
    for name, expected in zip(["a", "b", "c"], SAMPLE_SIZES):
        observed = samples.loc[samples["sample"] == name, "fragment_bp"].to_numpy()
        assert len(observed) == len(expected)
        assert np.allclose(observed, expected, rtol=0.02)
        weights = samples.loc[samples["sample"] == name, "relative_intensity"]
        assert weights.sum() == pytest.approx(1.0)
    assert (tmp_path / "out" / "qc_overlay.png").exists()
    assert (tmp_path / "out" / "analysis.json").exists()


def test_dark_gel_polarity(synthetic_gel: Path, tmp_path: Path):
    import cv2
    inverted = cv2.bitwise_not(cv2.imread(str(synthetic_gel)))
    path = tmp_path / "inverted.png"
    cv2.imwrite(str(path), inverted)
    result = analyze_gel(
        path, "HaeIII", tmp_path / "dark", 4, ladder_lane=1,
        ladder_sizes=LADDER_SIZES, polarity="dark",
    )
    assert len(result.bands) == 14


def test_manual_ladder_positions(synthetic_gel: Path, tmp_path: Path):
    positions = [migration_y(size) for size in LADDER_SIZES]
    result = analyze_gel(
        synthetic_gel, "HaeIII", tmp_path / "manual", 4,
        ladder_lane=1, ladder_sizes=LADDER_SIZES, ladder_y=positions,
    )
    ladder = result.bands.loc[result.bands["sample"] == "ladder"]
    assert np.allclose(ladder["fragment_bp"], LADDER_SIZES, rtol=0.02)


def test_ladder_count_mismatch_fails_loudly():
    with pytest.raises(ValueError, match="Detected 2 ladder bands"):
        calibrate_ladder([100, 200], [1000, 500, 100])


def test_ladder_must_be_monotonic():
    with pytest.raises(ValueError, match="strictly decrease"):
        calibrate_ladder([100, 200], [100, 1000])


def test_uncalibrated_gel_keeps_pixel_positions(synthetic_gel: Path, tmp_path: Path):
    result = analyze_gel(synthetic_gel, "HaeIII", tmp_path / "pixels", 4)
    assert result.bands["fragment_bp"].isna().all()
    assert result.bands["pixel_y"].notna().all()


def test_dual_ladder_corrects_linear_smile(tmp_path: Path):
    import cv2
    height, width = 720, 680
    centers = [120, 340, 560]
    expected = [1000, 500, 200]
    yy, xx = np.mgrid[:height, :width]
    image = np.full((height, width), 12.0)
    for lane, center in enumerate(centers):
        sizes = LADDER_SIZES if lane != 1 else expected
        offset = 30 * lane / 2
        for size in sizes:
            y = migration_y(size) + offset
            image += 220 * np.exp(-0.5 * ((xx - center) / 28) ** 2 - 0.5 * ((yy - y) / 3) ** 2)
    path = tmp_path / "smile.png"
    cv2.imwrite(str(path), np.clip(image, 0, 255).astype(np.uint8))
    result = analyze_gel(
        path, "HaeIII", tmp_path / "out", 3, lane_centers=centers,
        ladder_lane=1, secondary_ladder_lane=3, ladder_sizes=LADDER_SIZES,
        ladder_y=[migration_y(size) for size in LADDER_SIZES],
        secondary_ladder_y=[migration_y(size) + 30 for size in LADDER_SIZES],
        sample_names=["middle"],
    )
    observed = result.bands.loc[result.bands["sample"] == "middle", "fragment_bp"]
    assert np.allclose(observed, expected, rtol=0.02)
