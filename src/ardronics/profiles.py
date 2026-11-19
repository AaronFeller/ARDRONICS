"""Shared, strictly numeric fragment-profile I/O."""

from pathlib import Path

import numpy as np
import pandas as pd


def read_profiles(path: str | Path) -> pd.DataFrame:
    """Read tidy fragment_bp profiles or legacy Sample/enznames/1..25 CSVs."""
    frame = pd.read_csv(path)
    if "fragment_bp" in frame.columns:
        required = {"sample", "enzyme", "fragment_bp"}
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"Missing profile columns: {sorted(missing)}")
        frame = frame.copy()
        frame["fragment_bp"] = pd.to_numeric(frame["fragment_bp"], errors="raise")
    else:
        sample_col = "sample" if "sample" in frame else "Sample"
        enzyme_col = "enzyme" if "enzyme" in frame else "enznames"
        if sample_col not in frame or enzyme_col not in frame:
            raise ValueError("Profiles need sample/Sample and enzyme/enznames columns")
        fragment_cols = [c for c in frame.columns if str(c).isdigit()]
        if not fragment_cols:
            raise ValueError("No numeric fragment columns found")
        long = frame.assign(_row=np.arange(len(frame))).melt(
            id_vars=["_row", sample_col, enzyme_col], value_vars=fragment_cols,
            var_name="band_number", value_name="fragment_bp",
        ).sort_values(["_row", "band_number"])
        frame = long.rename(columns={sample_col: "sample", enzyme_col: "enzyme"})
        frame["fragment_bp"] = pd.to_numeric(frame["fragment_bp"], errors="coerce")
        frame = frame.dropna(subset=["fragment_bp"])
    frame = frame.dropna(subset=["sample", "enzyme", "fragment_bp"])
    if frame.empty:
        raise ValueError(f"No fragments found in {path}")
    if (frame["fragment_bp"] <= 0).any():
        raise ValueError("Fragment lengths must be positive numbers")
    frame["sample"] = frame["sample"].astype(str)
    frame["enzyme"] = frame["enzyme"].astype(str)
    columns = ["sample", "enzyme", "fragment_bp"]
    if "description" in frame:
        frame["description"] = frame["description"].fillna("").astype(str)
        columns.append("description")
    return frame[columns].reset_index(drop=True)


def write_profiles(frame: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, float_format="%.3f")
