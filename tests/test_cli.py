from pathlib import Path

import pandas as pd
from conftest import LADDER_SIZES, SAMPLE_SIZES

from ardronics.cli import main


def amplicon_for_fragments(fragments):
    sequence = list("A" * sum(fragments))
    boundary = 0
    for fragment in fragments[:-1]:
        boundary += fragment
        sequence[boundary - 2:boundary + 2] = list("GGCC")
    sequence[:4] = list("AACG")
    sequence[-4:] = list("ACGT")
    return "".join(sequence)


def test_complete_cli_run(synthetic_gel: Path, tmp_path: Path):
    fasta = tmp_path / "references.fa"
    with fasta.open("w") as output:
        for name, fragments in zip(["reference_a", "reference_b", "reference_c"], SAMPLE_SIZES):
            output.write(f">{name}\n{amplicon_for_fragments(fragments)}\n")
    output_dir = tmp_path / "run"
    status = main([
        "run", "--image", str(synthetic_gel), "--fasta", str(fasta),
        "--enzyme", "HaeIII", "--enzymes", "HaeIII", "--n-lanes", "4",
        "--ladder-lane", "1", "--ladder-sizes", ",".join(map(str, LADDER_SIZES)),
        "--sample-names", "a,b,c", "--forward-primer", "AACG", "--reverse-primer", "ACGT",
        "--primer-mismatches", "0", "--tolerance", "30", "--output-dir", str(output_dir),
    ])
    assert status == 0
    result = pd.read_csv(output_dir / "identification.csv")
    best = result.loc[result["rank"] == 1].set_index("sample")
    assert best.loc["a", "library_match"] == "reference_a"
    assert best.loc["b", "library_match"] == "reference_b"
    assert best.loc["c", "library_match"] == "reference_c"
    assert (best["jaccard_distance"] == 0).all()
    for filename in ["report.html", "library.csv", "library_qc.csv", "dendrogram_HaeIII.png"]:
        assert (output_dir / filename).exists()
    assert (output_dir / "gel" / "qc_overlay.png").exists()


def test_complete_run_with_lane_map_and_controls(synthetic_gel: Path, tmp_path: Path):
    fasta = tmp_path / "references.fa"
    with fasta.open("w") as output:
        for name, fragments in zip(["reference_a", "reference_b", "reference_c"], SAMPLE_SIZES):
            output.write(f">{name}\n{amplicon_for_fragments(fragments)}\n")
    lane_map = tmp_path / "lanes.csv"
    lane_map.write_text(
        "lane,role,sample,sample_type,expected_match\n"
        "1,ladder,,,\n"
        "2,control,a,isolate,reference_a\n"
        "3,sample,b,mixed,\n"
        "4,control,c,isolate,reference_c\n"
    )
    metadata = tmp_path / "metadata.json"
    metadata.write_text('{"gel_percent": 1.5, "voltage": 100}')
    output_dir = tmp_path / "mapped"
    import json
    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "fasta": str(fasta), "output_dir": str(output_dir),
        "gel": {
            "image": str(synthetic_gel), "enzyme": "HaeIII", "lane_map": str(lane_map),
            "ladder_sizes": LADDER_SIZES, "experiment_metadata": str(metadata),
        },
        "digest": {
            "enzymes": "HaeIII", "forward_primer": "AACG", "reverse_primer": "ACGT",
            "primer_mismatches": 0,
        },
        "matching": {"tolerance": 30},
    }))
    status = main(["run", "--config", str(config)])
    assert status == 0
    controls = pd.read_csv(output_dir / "control_validation.csv")
    assert controls["passed"].all()
    bands = pd.read_csv(output_dir / "gel" / "bands.csv")
    assert bands.loc[bands["sample"] == "a", "role"].unique().tolist() == ["control"]
    assert bands.loc[bands["sample"] == "b", "sample_type"].unique().tolist() == ["mixed"]
    import json
    analysis = json.loads((output_dir / "gel" / "analysis.json").read_text())
    assert analysis["experiment"]["gel_percent"] == 1.5
