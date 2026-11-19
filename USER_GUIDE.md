# ARDRONICS User Guide

ARDRONICS estimates DNA fragment sizes from gel images, simulates 16S restriction digests from reference sequences, and ranks matching organisms. It supports single/double digests, two-ladder smile correction, positive controls, and multi-enzyme consensus.

**Important:** an in-silico organism library replaces a known-organism control, not a molecular-weight ladder. Without sizing information, gel pixels cannot be reliably converted to base pairs.

## 1. Install and check

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -e '.[dev]'
.venv/bin/ardronics --help
.venv/bin/python -m pytest -q
```

Python 3.11 or newer is required. A standard `venv` plus `.venv/bin/pip install -e '.[dev]'` also works.

## 2. Recommended workflow

Prepare:

1. `gel.jpg` — original gel image.
2. `config.json` — image, ladder, enzyme, primer, FASTA, and matching settings.
3. `lanes.csv` — lane assignments and controls.
4. `references.fa` — reference 16S sequences.
5. Optional `metadata.json` for experimental notes.

Copy and edit the templates:

```bash
cp config.example.json my_config.json
cp examples/lane_map.csv my_lanes.csv
cp examples/experiment_metadata.json my_metadata.json
.venv/bin/ardronics run --config my_config.json
```

Explicit CLI options override JSON values. Paths are relative to your current working directory.

## 3. Input checklist

| Input | JSON location | CLI equivalent | Required? |
|---|---|---|---|
| Gel image | `gel.image` | `--image` | Yes |
| Lane map | `gel.lane_map` | `--lane-map` | Recommended |
| Lane count | `gel.n_lanes` | `--n-lanes` | If no lane map |
| Sample names | `gel.sample_names` | `--sample-names` | If no lane map |
| Ladder preset | `gel.ladder` | `--ladder` | For bp sizing |
| Custom ladder sizes | `gel.ladder_sizes` | `--ladder-sizes` | Alternative to preset |
| Primary ladder lane | `gel.ladder_lane` or lane map | `--ladder-lane` | For bp sizing |
| Secondary ladder | `gel.secondary_ladder_lane` or lane map | `--secondary-ladder-lane` | Optional |
| Gel digest enzyme | `gel.enzyme` | `--enzyme` | For matching |
| Library enzymes | `digest.enzymes` | `--enzymes` | Yes |
| Primers | `digest.forward_primer`, `digest.reverse_primer` | Corresponding flags | Yes |
| Reference FASTA | top-level `fasta` | `--fasta` | For matching |
| Matching tolerance | `matching.tolerance` | `--tolerance` | Default 20 bp |
| Minimum fragment | `matching.min_bp` | `--min-bp` | Default 100 bp |
| Number of matches | `matching.top_n` | `--top-n` | Default 3 |
| Multi-enzyme consensus | `matching.combine_enzymes` | `--combine-enzymes` | Optional |
| Experiment metadata | `gel.experiment_metadata` | `--experiment-metadata` | Optional |

A lane map is authoritative. Do not also supply `n_lanes`, `sample_names`, or ladder-lane options for that run.

## 4. Lane map

```csv
lane,role,sample,sample_type,expected_match
1,ladder,,,
2,control,control_a,isolate,reference_a
3,sample,unknown_b,isolate,
4,empty,,,
5,unused,MW_right,,
6,ladder_secondary,,,
```

List every lane consecutively from left to right, starting at 1, including empty and marker lanes. `sample_type` is `isolate`, `mixed`, or `unknown`. `expected_match` is optional and only used for controls.

| Role | Meaning | Included in identification? |
|---|---|---|
| `sample` | Unknown experimental sample | Yes |
| `control` | Known positive control | Yes, plus validation |
| `ladder` | Primary molecular-weight ladder | No |
| `ladder_secondary` | Opposite ladder for smile correction | No |
| `empty` | Intentionally empty lane | No |
| `failed` | Failed or unusable reaction | No |
| `unused` | Retain in QC but exclude from matching | No |

Roles and sample types are retained in `bands.csv` and `analysis.json`.

### Positive controls

Set `expected_match` to an accession or distinctive text from the FASTA description. Matching is case-insensitive. A complete run writes `control_validation.csv` with expected/observed rank-one matches, distance, and pass/fail.

A failed control means you should investigate ladder assignment, enzyme/primer metadata, reference choice, and tolerance before trusting unknowns.


## 5. Ladder and sizing

Presets are `50bp`, `100bp`, and `1kb`. For another ladder:

```json
"ladder_sizes": [2072, 1500, 1000, 900, 800, 700, 600, 500, 400, 300, 250, 200, 100]
```

Or use `--ladder-sizes 2072,1500,1000,900,800,700,600,500,400,300,250,200,100`.

Provide sizes in **top-to-bottom migration order**, normally largest to smallest. The detected ladder-band count must equal the size count. ARDRONICS stops on a mismatch instead of assigning sizes silently.

Use two ladders when available. Each lane then receives calibration interpolated between the references. If automatic detection includes artifacts, inspect `qc_overlay.png` and supply reviewed positions after cropping:

```bash
--ladder-y 60,71,82,88,94,101,110,120,133,146,164,186,216
```

For a second ladder, use `--secondary-ladder-y` with the same number of positions.

### No usable ladder

Run `ardronics gel` without ladder sizes to obtain pixel positions and intensities. Marker roles in a lane map are allowed without sizes. Do not use these pixel-only results for organism identification against bp profiles.

A known-organism lane can sometimes substitute for a ladder if its expected fragments and band assignments are independently convincing. Use those fragments as `ladder_sizes` and reviewed positions as `ladder_y`. This is circular if the same assignment is then used as evidence for that organism; retain an independent control whenever possible.

## 6. Enzymes, primers, and FASTA references

The gel enzyme describes the wet-lab digest. The library enzyme list determines which theoretical profiles are generated:

```json
"digest": {
  "enzymes": "HaeIII,AluI,RsaI",
  "forward_primer": "AGAGTTTGATCMTGGCTCAG",
  "reverse_primer": "GGTTACCTTGTTACGACTT",
  "primer_mismatches": 2,
  "double_digest": ["HaeIII+HhaI"]
}
```

Use the actual PCR primers. They determine amplicon boundaries and terminal fragments. Restriction cut positions come from Biopython's enzyme definitions.

FASTA identifiers must be unique. Put the organism and strain in the description:

```text
>NR_117677.1 Acinetobacter baumannii DSM 30007
ACGT...
```

Prefer sequences containing the complete primer-defined amplicon. By default, records without both primers fail QC. `--allow-unclipped` instead digests the available sequence and marks it `ok_unclipped`; internal fragments remain useful, but terminal fragments are approximate.

For a public example library:

```bash
.venv/bin/python scripts/build_organism_library.py
```

## 7. Image corrections and band detection

| Option | Purpose |
|---|---|
| `--crop x,y,width,height` | Remove labels, captions, and unrelated panels |
| `--polarity dark` | Analyze dark bands on a light background |
| `--lane-centers 50,150,250` | Manually set post-crop lane x positions |
| `--prominence 2.0` | Detect weaker bands; higher is more conservative |
| `--min-band-distance 3` | Minimum peak spacing in pixels; adjust close doublets |
| `--exclude-sample NAME` | Exclude a marker/sample during standalone matching; repeatable |

Always inspect `qc_overlay.png`. A regular-grid lane fallback is a warning to verify lane alignment, not confirmation that alignment is correct.

Optional `metadata.json` can record gel percentage, matrix, voltage, run time, operator, and notes. These values are saved in `analysis.json` but do not alter calibration.

## 8. Matching and multiple enzymes

Standalone matching:

```bash
.venv/bin/ardronics identify \
  --samples results/gel/profiles.csv \
  --library results/library.csv \
  --tolerance 20 --min-bp 100 --top-n 3 \
  --output-dir results/identification
```

The input is tidy `sample,enzyme,fragment_bp` CSV. Legacy wide `Sample,enznames,1,...,25` CSV is also accepted. Each observed band can match only one reference band within tolerance.

### Combine independent digests

Use the same sample names in every gel's `profiles.csv`, then combine the files:

```bash
.venv/bin/python -c "import pandas as pd; pd.concat([pd.read_csv('results/haeiii/profiles.csv'), pd.read_csv('results/alui/profiles.csv')]).to_csv('results/combined_profiles.csv', index=False)"

.venv/bin/ardronics identify \
  --samples results/combined_profiles.csv \
  --library results/organism_library/library.csv \
  --combine-enzymes --tolerance 20 \
  --exclude-sample MW_left --exclude-sample MW_right \
  --output-dir results/consensus
```

Consensus averages the per-enzyme Jaccard distances with equal weight per enzyme. Every reference must have a profile for every observed enzyme; incomplete references are excluded. Dendrograms remain per enzyme.

Multiple enzymes can resolve single-enzyme ties, but organisms identical across all tested enzymes remain ambiguous. Definitive separation may require additional enzymes or sequencing.

Mixed-community lanes are labeled as `mixed`, but matching a library of isolates is not a community-deconvolution method.



## 9. Outputs and interpretation

| File | Contents |
|---|---|
| `gel/bands.csv` | Lane, role, type, pixel position, bp, intensity, relative intensity, prominence |
| `gel/profiles.csv` | Tidy profiles used for matching |
| `gel/qc_overlay.png` | Detected lanes and bands for visual review |
| `gel/analysis.json` | Settings, metadata, lane roles, ladder coordinates, warnings |
| `library.csv` | Theoretical fragments and FASTA descriptions |
| `library_qc.csv` | `ok`, `ok_unclipped`, or `failed` per reference |
| `identification.csv` | Ranked matches, distances, band counts, ties, and margins |
| `control_validation.csv` | Positive-control pass/fail, when expectations were supplied |
| `dendrogram_ENZYME.png` | UPGMA clustering for each enzyme |
| `report.html` | Combined report for a complete run |

Jaccard distance is 0 for identical matched profiles and 1 for no matching bands. Similarity is `1 - distance`; it is not a probability of correct identification.

`ambiguous=true` means references are tied or the next candidate has the same distance. `margin_to_next` shows the separation from the next candidate. Do not treat a small or zero margin as a confident unique identification.

## 10. Troubleshooting

**“Detected N ladder bands but M sizes were supplied”**  
Inspect the ladder overlay. Correct the size list, crop, sensitivity, or provide reviewed `ladder_y` coordinates. Do not truncate the list merely to force a match.

**Lane-grid fallback or shifted markings**  
Use a tighter crop or explicit `lane_centers`. Verify both outer markers and several central lanes.

**Primers not found / failed library records**  
Check sequence orientation, completeness, primer sequences, and mismatch allowance. Use exact amplicon references when possible. Use `allow_unclipped` only with awareness of terminal-fragment uncertainty.

**No candidates in multi-enzyme consensus**  
Ensure every reference has every required enzyme profile and that enzyme names match exactly, including capitalization.

**Unexpected organism matches**  
Check the enzyme, ladder identity/order, primers, control validation, reference coverage, and tolerance. Do not tune settings solely to obtain the expected species.

**Everything is ambiguous**  
Add an independent enzyme, improve reference coverage, or use sequencing. Image processing cannot separate organisms with identical restriction profiles.

## 11. Supplied worked examples

```bash
# Mycoplasma example: see README.md for its full command
scripts/run_more_examples.sh

# Public NCBI reference library
.venv/bin/python scripts/build_organism_library.py

# Rebuild from downloaded FASTA without network access
.venv/bin/python scripts/build_organism_library.py --skip-download
```

The supplied examples demonstrate mechanics and QC. Their matches are not a substitute for validation on your own positive controls.

