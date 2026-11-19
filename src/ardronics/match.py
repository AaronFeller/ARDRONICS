"""Numerically safe ARDRA profile matching and clustering."""

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import squareform


def match_band_count(observed: np.ndarray, reference: np.ndarray, tolerance_bp: float) -> int:
    """Maximum one-to-one band assignment within tolerance, minimizing total error."""
    observed = np.asarray(observed, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if not len(observed) or not len(reference):
        return 0
    differences = np.abs(observed[:, None] - reference[None, :])
    rows, columns = linear_sum_assignment(differences)
    return int(np.sum(differences[rows, columns] <= tolerance_bp))


def jaccard_distance(observed: np.ndarray, reference: np.ndarray, tolerance_bp: float) -> tuple[float, int]:
    matched = match_band_count(observed, reference, tolerance_bp)
    union = len(observed) + len(reference) - matched
    return (1.0 - matched / union if union else 1.0), matched


def _filtered_profiles(frame: pd.DataFrame, min_bp: float) -> dict[tuple[str, str], np.ndarray]:
    clean = frame.loc[frame["fragment_bp"] >= min_bp]
    return {
        (sample, enzyme): group["fragment_bp"].to_numpy(dtype=float)
        for (sample, enzyme), group in clean.groupby(["sample", "enzyme"], sort=False)
    }


def identify_samples(
    samples: pd.DataFrame,
    library: pd.DataFrame,
    tolerance_bp: float = 20,
    min_bp: float = 100,
    top_n: int = 3,
) -> pd.DataFrame:
    """Rank reference profiles, collapsing references with identical fragments."""
    if tolerance_bp < 0 or min_bp < 0 or top_n < 1:
        raise ValueError("Tolerance/min_bp must be nonnegative and top_n positive")
    observed = _filtered_profiles(samples, min_bp)
    references = _filtered_profiles(library, min_bp)
    groups: dict[tuple[str, tuple[float, ...]], list[str]] = {}
    descriptions = (
        library.groupby("sample")["description"].first().to_dict()
        if "description" in library else {}
    )
    for (sample, enzyme), fragments in references.items():
        key = (enzyme, tuple(sorted(np.round(fragments, 3))))
        groups.setdefault(key, []).append(sample)

    rows = []
    for (sample, enzyme), fragments in observed.items():
        candidates = []
        for (reference_enzyme, reference_key), members in groups.items():
            if enzyme != reference_enzyme:
                continue
            reference = np.asarray(reference_key)
            distance, matched = jaccard_distance(fragments, reference, tolerance_bp)
            candidates.append({
                "sample": sample, "enzyme": enzyme, "library_match": ";".join(sorted(members)),
                "library_description": ";".join(descriptions.get(member, "") for member in sorted(members)),
                "jaccard_distance": distance, "similarity_percent": (1 - distance) * 100,
                "matched_bands": matched, "sample_bands": len(fragments),
                "reference_bands": len(reference), "references_in_group": len(members),
            })
        candidates.sort(key=lambda row: (row["jaccard_distance"], row["library_match"]))
        for rank, candidate in enumerate(candidates[:top_n], 1):
            next_distance = candidates[rank]["jaccard_distance"] if rank < len(candidates) else None
            candidate["rank"] = rank
            candidate["ambiguous"] = bool(
                candidate["references_in_group"] > 1
                or (next_distance is not None and abs(next_distance - candidate["jaccard_distance"]) < 1e-12)
            )
            candidate["margin_to_next"] = (
                next_distance - candidate["jaccard_distance"] if next_distance is not None else np.nan
            )
            rows.append(candidate)
    columns = [
        "sample", "enzyme", "rank", "library_match", "library_description", "jaccard_distance", "similarity_percent",
        "matched_bands", "sample_bands", "reference_bands", "references_in_group", "ambiguous", "margin_to_next",
    ]
    return pd.DataFrame(rows, columns=columns)


def identify_samples_consensus(
    samples: pd.DataFrame,
    library: pd.DataFrame,
    tolerance_bp: float = 20,
    min_bp: float = 100,
    top_n: int = 3,
) -> pd.DataFrame:
    """Combine evidence across enzymes, requiring complete reference coverage."""
    reference_count = library["sample"].nunique()
    detailed = identify_samples(samples, library, tolerance_bp, min_bp, max(1, reference_count))
    if detailed.empty:
        return detailed.assign(enzyme_coverage=pd.Series(dtype=float))
    detailed = detailed.assign(library_match=detailed["library_match"].str.split(";")).explode("library_match")
    descriptions = library.groupby("sample")["description"].first().to_dict() if "description" in library else {}
    output = []
    for sample, sample_rows in detailed.groupby("sample", sort=False):
        enzymes = sorted(samples.loc[samples["sample"] == sample, "enzyme"].unique())
        signatures = {}
        evidence = {}
        for member, member_rows in sample_rows.groupby("library_match"):
            if set(member_rows["enzyme"]) != set(enzymes) or len(member_rows) != len(enzymes):
                continue
            ordered = member_rows.set_index("enzyme").loc[enzymes]
            signature = tuple(
                (round(row.jaccard_distance, 12), int(row.matched_bands), int(row.reference_bands))
                for row in ordered.itertuples()
            )
            signatures.setdefault(signature, []).append(member)
            evidence[signature] = ordered
        candidates = []
        for signature, members in signatures.items():
            ordered = evidence[signature]
            distance = float(ordered["jaccard_distance"].mean())
            candidates.append({
                "sample": sample, "enzyme": "+".join(enzymes),
                "library_match": ";".join(sorted(members)),
                "library_description": ";".join(descriptions.get(member, "") for member in sorted(members)),
                "jaccard_distance": distance, "similarity_percent": (1 - distance) * 100,
                "matched_bands": int(ordered["matched_bands"].sum()),
                "sample_bands": int(ordered["sample_bands"].sum()),
                "reference_bands": int(ordered["reference_bands"].sum()),
                "references_in_group": len(members), "enzyme_coverage": 1.0,
            })
        candidates.sort(key=lambda row: (row["jaccard_distance"], row["library_match"]))
        for rank, candidate in enumerate(candidates[:top_n], 1):
            next_distance = candidates[rank]["jaccard_distance"] if rank < len(candidates) else None
            candidate["rank"] = rank
            candidate["ambiguous"] = bool(candidate["references_in_group"] > 1 or (next_distance is not None and abs(next_distance - candidate["jaccard_distance"]) < 1e-12))
            candidate["margin_to_next"] = next_distance - candidate["jaccard_distance"] if next_distance is not None else np.nan
            output.append(candidate)
    ordered_columns = ["sample", "enzyme", "rank", "library_match", "library_description", "jaccard_distance", "similarity_percent", "matched_bands", "sample_bands", "reference_bands", "references_in_group", "ambiguous", "margin_to_next", "enzyme_coverage"]
    return pd.DataFrame(output, columns=ordered_columns)


def validate_controls(report: pd.DataFrame, expectations: dict[str, str]) -> pd.DataFrame:
    """Check rank-one control calls against accession or description text."""
    rows = []
    best = report.loc[report["rank"] == 1].set_index("sample") if not report.empty else report
    for sample, expected in expectations.items():
        observed, distance, passed = "", np.nan, False
        if not best.empty and sample in best.index:
            match = best.loc[sample]
            observed = str(match["library_match"])
            distance = float(match["jaccard_distance"])
            passed = expected.casefold() in (observed + " " + str(match["library_description"])).casefold()
        rows.append({"sample": sample, "expected": expected, "observed": observed, "jaccard_distance": distance, "passed": passed})
    return pd.DataFrame(rows)


def write_dendrograms(samples: pd.DataFrame, library: pd.DataFrame, output_dir: str | Path, tolerance_bp: float, min_bp: float) -> list[Path]:
    """Write one UPGMA dendrogram per enzyme represented in the unknowns."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    observed = _filtered_profiles(samples, min_bp)
    references = _filtered_profiles(library, min_bp)
    overlap = set(observed).intersection(references)
    unknown_enzymes = set(samples["enzyme"])
    paths = []
    for enzyme in sorted(unknown_enzymes):
        selected = []
        for source, profiles in (("sample", observed), ("reference", references)):
            for (name, enz), fragments in profiles.items():
                if enz == enzyme:
                    label = f"{source}:{name}" if (name, enz) in overlap else name
                    selected.append((label, fragments))
        if len(selected) < 2:
            continue
        labels = [name for name, _ in selected]
        matrix = np.zeros((len(selected), len(selected)))
        for row in range(len(selected)):
            for column in range(row + 1, len(selected)):
                distance, _ = jaccard_distance(selected[row][1], selected[column][1], tolerance_bp)
                matrix[row, column] = matrix[column, row] = distance
        condensed = squareform(matrix, checks=False)
        tree = linkage(condensed, method="average")
        fig, axis = plt.subplots(figsize=(max(8, len(labels) * 0.3), 6))
        dendrogram(tree, labels=labels, ax=axis, leaf_rotation=90, leaf_font_size=8)
        axis.set_title(f"ARDRONICS UPGMA — {enzyme}")
        axis.set_ylabel("Jaccard distance")
        fig.tight_layout()
        safe_enzyme = re.sub(r"[^A-Za-z0-9_.-]+", "_", enzyme)
        path = output_dir / f"dendrogram_{safe_enzyme}.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(path)
    return paths
