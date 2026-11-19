from pathlib import Path

import pandas as pd
import pytest

from ardronics.digest import (
    FORWARD_27F,
    REVERSE_1492R,
    digest_sequence,
    extract_amplicon,
    generate_library,
    primer_matches,
)


def test_iupac_primer_matches_degenerate_base():
    sequence = "TTTAGAGTTTGATCATGGCTCAGAAA"
    matches = primer_matches(sequence, FORWARD_27F, max_mismatches=0)
    assert matches == [(3, 0)]


def test_amplicon_includes_primers_and_haeiii_cut_offset():
    sequence = "TTAACGAAGGCCAAACGTTT"
    amplicon = extract_amplicon(sequence, "AACG", "ACGT", max_mismatches=0)
    assert amplicon == "AACGAAGGCCAAACGT"
    assert len(amplicon) == 16
    # HaeIII recognizes GGCC and cuts GG^CC, producing two 8 bp pieces.
    assert digest_sequence(amplicon, ["HaeIII"]) == [8, 8]


def test_double_digest_unites_cut_sites():
    sequence = "AACGAAGGCCGCGCAAACGT"
    single = digest_sequence(sequence, ["HaeIII"])
    double = digest_sequence(sequence, ["HaeIII", "HhaI"])
    assert len(double) == len(single) + 1
    assert sum(double) == len(sequence)


def test_library_is_numeric_and_reports_primer_failures(tmp_path: Path):
    fasta = tmp_path / "references.fa"
    fasta.write_text(
        ">known description here\nTTAACGAAGGCCAAACGTTT\n"
        ">no_primers\nAAAAAAAAAAAAAAAAAAAA\n"
    )
    output = tmp_path / "library.csv"
    library = generate_library(
        fasta, output, enzymes=["HaeIII"], forward_primer="AACG",
        reverse_primer="ACGT", max_primer_mismatches=0,
    )
    assert library["fragment_bp"].tolist() == [8, 8]
    assert pd.api.types.is_numeric_dtype(library["fragment_bp"])
    qc = pd.read_csv(tmp_path / "library_qc.csv")
    assert dict(zip(qc["sample"], qc["status"])) == {"known": "ok", "no_primers": "failed"}


def test_no_valid_sequences_raises_with_qc(tmp_path: Path):
    fasta = tmp_path / "bad.fa"
    fasta.write_text(">bad\nAAAA\n")
    with pytest.raises(ValueError, match="No library profiles"):
        generate_library(fasta, tmp_path / "bad.csv", forward_primer="AACG", reverse_primer="ACGT")
    assert (tmp_path / "bad_qc.csv").exists()


def test_unclipped_mode_is_explicit_in_qc(tmp_path: Path):
    fasta = tmp_path / "partial.fa"
    fasta.write_text(">partial\nAACGAAGGCCAAACGT\n")
    library = generate_library(
        fasta, tmp_path / "partial.csv", enzymes=["HaeIII"],
        forward_primer=FORWARD_27F, reverse_primer=REVERSE_1492R,
        allow_unclipped=True,
    )
    assert library["fragment_bp"].tolist() == [8, 8]
    qc = pd.read_csv(tmp_path / "partial_qc.csv")
    assert qc.iloc[0]["status"] == "ok_unclipped"
