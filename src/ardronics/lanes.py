"""Validated lane manifests and positive-control metadata."""

import csv
from dataclasses import dataclass
from pathlib import Path

ROLES = {"sample", "control", "ladder", "ladder_secondary", "empty", "failed", "unused"}
SAMPLE_TYPES = {"isolate", "mixed", "unknown"}


@dataclass
class LaneMap:
    n_lanes: int
    ladder_lane: int | None
    secondary_ladder_lane: int | None
    sample_names: list[str]
    roles: dict[int, str]
    sample_types: dict[int, str]
    excluded_samples: list[str]
    expected_matches: dict[str, str]


def read_lane_map(path: str | Path) -> LaneMap:
    with Path(path).open(newline="") as source:
        rows = list(csv.DictReader(source))
    if not rows:
        raise ValueError("Lane map is empty")
    required = {"lane", "role", "sample", "sample_type"}
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"Lane map missing columns: {sorted(missing)}")

    entries = {}
    for row in rows:
        try:
            lane = int(row["lane"])
        except (TypeError, ValueError) as error:
            raise ValueError(f"Invalid lane number: {row['lane']}") from error
        role = (row["role"] or "sample").strip()
        sample_type = (row["sample_type"] or "unknown").strip()
        sample = (row["sample"] or "").strip()
        if lane in entries:
            raise ValueError(f"Duplicate lane-map row for lane {lane}")
        if role not in ROLES:
            raise ValueError(f"Lane {lane}: role must be one of {sorted(ROLES)}")
        if sample_type not in SAMPLE_TYPES:
            raise ValueError(f"Lane {lane}: sample_type must be isolate, mixed, or unknown")
        if role in {"sample", "control"} and not sample:
            raise ValueError(f"Lane {lane}: {role} requires a sample name")
        if not sample:
            sample = f"{role}_lane_{lane}"
        expected = (row.get("expected_match") or "").strip()
        entries[lane] = {"role": role, "sample": sample, "sample_type": sample_type, "expected": expected}

    lanes = sorted(entries)
    if lanes != list(range(1, len(lanes) + 1)):
        raise ValueError("Lane map must contain every lane consecutively starting at 1")
    ladders = [lane for lane in lanes if entries[lane]["role"] == "ladder"]
    secondary = [lane for lane in lanes if entries[lane]["role"] == "ladder_secondary"]
    if len(ladders) > 1 or len(secondary) > 1:
        raise ValueError("Lane map permits one primary and one secondary ladder")
    non_ladders = [lane for lane in lanes if entries[lane]["role"] not in {"ladder", "ladder_secondary"}]
    names = [entries[lane]["sample"] for lane in non_ladders]
    if len(set(names)) != len(names):
        raise ValueError("Sample names in the lane map must be unique")
    return LaneMap(
        n_lanes=len(lanes), ladder_lane=ladders[0] if ladders else None,
        secondary_ladder_lane=secondary[0] if secondary else None,
        sample_names=names,
        roles={lane: entries[lane]["role"] for lane in lanes},
        sample_types={lane: entries[lane]["sample_type"] for lane in lanes},
        excluded_samples=[entries[lane]["sample"] for lane in lanes if entries[lane]["role"] in {"empty", "failed", "unused"}],
        expected_matches={
            entries[lane]["sample"]: entries[lane]["expected"] for lane in lanes
            if entries[lane]["role"] == "control" and entries[lane]["expected"]
        },
    )
