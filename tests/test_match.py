import pandas as pd
import pytest

from ardronics.match import (
    identify_samples,
    identify_samples_consensus,
    jaccard_distance,
    match_band_count,
)


def profile(sample, fragments, enzyme="HaeIII"):
    return pd.DataFrame({"sample": sample, "enzyme": enzyme, "fragment_bp": fragments})


def test_numeric_matching_is_not_lexicographic():
    observed = profile("unknown", [99, 500, 1200])
    reference = profile("known", [100, 490, 1210])
    result = identify_samples(observed, reference, tolerance_bp=20, min_bp=100)
    assert result.iloc[0]["matched_bands"] == 2
    assert result.iloc[0]["sample_bands"] == 2


def test_one_to_one_matching_prevents_double_counting():
    assert match_band_count([100, 102], [110], tolerance_bp=20) == 1
    distance, matched = jaccard_distance([100, 102], [110], tolerance_bp=20)
    assert matched == 1
    assert distance == pytest.approx(0.5)


def test_identical_references_are_grouped_and_flagged():
    observed = profile("unknown", [200, 500])
    references = pd.concat([profile("taxon_a", [200, 500]), profile("taxon_b", [200, 500])])
    result = identify_samples(observed, references, tolerance_bp=10)
    best = result.iloc[0]
    assert best["jaccard_distance"] == 0
    assert best["references_in_group"] == 2
    assert best["ambiguous"]


def test_enzymes_are_not_cross_matched():
    observed = profile("unknown", [200, 500], enzyme="HhaI")
    references = profile("known", [200, 500], enzyme="HaeIII")
    assert identify_samples(observed, references).empty


def test_out_of_tolerance_has_distance_one():
    distance, matched = jaccard_distance([200, 500], [800], tolerance_bp=20)
    assert matched == 0
    assert distance == 1


def test_library_descriptions_are_retained():
    references = profile("accession_1", [200, 500])
    references["description"] = "accession_1 Example organism strain A"
    result = identify_samples(profile("unknown", [200, 500]), references)
    assert result.iloc[0]["library_description"] == "accession_1 Example organism strain A"


def test_multi_enzyme_consensus_resolves_single_enzyme_ties():
    observed = pd.concat([
        profile("unknown", [200, 500], "HaeIII"),
        profile("unknown", [300, 700], "AluI"),
    ])
    references = pd.concat([
        profile("taxon_a", [200, 500], "HaeIII"), profile("taxon_a", [900], "AluI"),
        profile("taxon_b", [850], "HaeIII"), profile("taxon_b", [300, 700], "AluI"),
        profile("taxon_c", [200, 500], "HaeIII"), profile("taxon_c", [300, 700], "AluI"),
    ])
    result = identify_samples_consensus(observed, references, tolerance_bp=20)
    best = result.iloc[0]
    assert best["library_match"] == "taxon_c"
    assert best["jaccard_distance"] == 0
    assert best["enzyme"] == "AluI+HaeIII"
    assert best["enzyme_coverage"] == 1


def test_consensus_rejects_references_missing_an_enzyme():
    observed = pd.concat([profile("unknown", [200], "HaeIII"), profile("unknown", [300], "AluI")])
    references = profile("incomplete", [200], "HaeIII")
    assert identify_samples_consensus(observed, references).empty
