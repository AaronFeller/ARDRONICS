from pathlib import Path

import pandas as pd

from ardronics.profiles import read_profiles


def test_legacy_wide_profiles_are_converted_to_numeric(tmp_path: Path):
    path = tmp_path / "legacy.csv"
    pd.DataFrame({
        "Sample": ["a", "b"], "enznames": ["HaeIII", "HaeIII"],
        "1": [99, 1200], "2": [500, None], "3": [None, None],
    }).to_csv(path, index=False)
    result = read_profiles(path)
    assert result["fragment_bp"].tolist() == [99, 500, 1200]
    assert pd.api.types.is_numeric_dtype(result["fragment_bp"])
