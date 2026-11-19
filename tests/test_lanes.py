from pathlib import Path

import pytest

from ardronics.lanes import read_lane_map


def test_example_lane_map_provides_all_lane_inputs():
    path = Path(__file__).parents[1] / "examples" / "lane_map.csv"
    result = read_lane_map(path)
    assert result.n_lanes == 6
    assert result.ladder_lane == 1
    assert result.secondary_ladder_lane == 6
    assert result.sample_names == ["control_a", "unknown_b", "empty_lane_4", "MW_right"]
    assert result.roles[2] == "control"
    assert result.sample_types[3] == "isolate"
    assert result.excluded_samples == ["empty_lane_4", "MW_right"]
    assert result.expected_matches == {"control_a": "reference_a"}


def test_lane_map_requires_consecutive_lanes(tmp_path: Path):
    path = tmp_path / "bad.csv"
    path.write_text("lane,role,sample,sample_type\n1,ladder,,,\n3,sample,x,isolate\n")
    with pytest.raises(ValueError, match="consecutively"):
        read_lane_map(path)


def test_lane_map_requires_names_for_samples(tmp_path: Path):
    path = tmp_path / "bad.csv"
    path.write_text("lane,role,sample,sample_type\n1,sample,,isolate\n")
    with pytest.raises(ValueError, match="requires a sample name"):
        read_lane_map(path)
