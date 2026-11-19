"""IUPAC-aware PCR extraction and in-silico restriction digestion."""

from collections.abc import Iterable
from itertools import pairwise
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Restriction import AllEnzymes
from Bio.Seq import Seq

FORWARD_27F = "AGAGTTTGATCMTGGCTCAG"
REVERSE_1492R = "GGTTACCTTGTTACGACTT"
IUPAC = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "U": set("T"), "R": set("AG"), "Y": set("CT"), "S": set("GC"),
    "W": set("AT"), "K": set("GT"), "M": set("AC"), "B": set("CGT"),
    "D": set("AGT"), "H": set("ACT"), "V": set("ACG"), "N": set("ACGT"),
}


def primer_matches(sequence: str, primer: str, max_mismatches: int = 2) -> list[tuple[int, int]]:
    """Return (start, mismatches), honoring degenerate bases in both inputs."""
    sequence = "".join(sequence.upper().split())
    primer = "".join(primer.upper().split())
    if not primer or len(primer) > len(sequence):
        return []
    unknown = set(sequence + primer) - set(IUPAC)
    if unknown:
        raise ValueError(f"Unsupported DNA characters: {sorted(unknown)}")
    matches = []
    for start in range(len(sequence) - len(primer) + 1):
        mismatches = sum(
            not IUPAC[sequence[start + offset]].intersection(IUPAC[base])
            for offset, base in enumerate(primer)
        )
        if mismatches <= max_mismatches:
            matches.append((start, mismatches))
    return matches


def extract_amplicon(
    sequence: str,
    forward_primer: str = FORWARD_27F,
    reverse_primer: str = REVERSE_1492R,
    max_mismatches: int = 2,
) -> str:
    """Extract the PCR product, including both primer sequences."""
    forward = primer_matches(sequence, forward_primer, max_mismatches)
    reverse = primer_matches(sequence, str(Seq(reverse_primer).reverse_complement()), max_mismatches)
    if not forward or not reverse:
        raise ValueError("Forward or reverse primer was not found")
    forward_start = min(forward, key=lambda item: (item[1], item[0]))[0]
    reverse_start = min(reverse, key=lambda item: (item[1], -item[0]))[0]
    reverse_end = reverse_start + len(reverse_primer)
    if reverse_end <= forward_start + len(forward_primer):
        raise ValueError("Primers do not define a valid amplicon")
    return sequence[forward_start:reverse_end]


def restriction_sites(sequence: str, enzyme_name: str) -> list[int]:
    """Return zero-based cut boundaries using Biopython's enzyme definitions."""
    enzyme_class = next((item for item in AllEnzymes if str(item) == enzyme_name), None)
    if enzyme_class is None:
        raise ValueError(f"Unknown restriction enzyme: {enzyme_name}")
    # Bio.Restriction reports biological (one-based) cut coordinates.
    return [position - 1 for position in enzyme_class.search(Seq(sequence), linear=True)]


def digest_sequence(sequence: str, enzyme_names: Iterable[str]) -> list[int]:
    cuts = sorted({site for name in enzyme_names for site in restriction_sites(sequence, name)})
    boundaries = [0, *cuts, len(sequence)]
    return [right - left for left, right in pairwise(boundaries) if right > left]


def generate_library(
    fasta_path: str | Path,
    output_csv: str | Path,
    enzymes: Iterable[str] = ("HaeIII", "HhaI"),
    forward_primer: str = FORWARD_27F,
    reverse_primer: str = REVERSE_1492R,
    max_primer_mismatches: int = 2,
    double_digests: Iterable[tuple[str, str]] = (),
    allow_unclipped: bool = False,
) -> pd.DataFrame:
    """Build tidy sample/enzyme/fragment_bp profiles and a sidecar QC report."""
    fasta_path = Path(fasta_path)
    output_csv = Path(output_csv)
    treatments: list[tuple[str, tuple[str, ...]]] = [(name, (name,)) for name in enzymes]
    treatments.extend(("+".join(pair), tuple(pair)) for pair in double_digests)
    if not treatments:
        raise ValueError("At least one restriction enzyme is required")

    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    identifiers = [record.id for record in records]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("FASTA record identifiers must be unique")

    rows: list[dict] = []
    qc_rows: list[dict] = []
    for record in records:
        try:
            status = "ok"
            try:
                amplicon = extract_amplicon(
                    str(record.seq), forward_primer, reverse_primer, max_primer_mismatches
                )
            except ValueError:
                if not allow_unclipped:
                    raise
                if not len(record.seq):
                    raise
                amplicon = str(record.seq)
                status = "ok_unclipped"
            for treatment, enzyme_pair in treatments:
                for fragment in digest_sequence(amplicon, enzyme_pair):
                    rows.append({
                        "sample": record.id, "enzyme": treatment,
                        "fragment_bp": fragment, "amplicon_bp": len(amplicon),
                        "description": record.description,
                    })
            qc_rows.append({"sample": record.id, "status": status, "amplicon_bp": len(amplicon), "message": ""})
        except (ValueError, TypeError) as error:
            qc_rows.append({"sample": record.id, "status": "failed", "amplicon_bp": pd.NA, "message": str(error)})

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    library = pd.DataFrame(rows, columns=["sample", "enzyme", "fragment_bp", "amplicon_bp", "description"])
    library.to_csv(output_csv, index=False)
    qc_path = output_csv.with_name(output_csv.stem + "_qc.csv")
    pd.DataFrame(qc_rows).to_csv(qc_path, index=False)
    if library.empty:
        raise ValueError(f"No library profiles generated; inspect {qc_path}")
    return library
