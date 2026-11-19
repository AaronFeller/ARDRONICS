# ARDRONICS

ARDRONICS is a Python-native Amplified Ribosomal DNA Restriction Analysis (ARDRA) tool. It reads gel bands, estimates fragment sizes and densitometry weights, generates theoretical 16S restriction profiles from FASTA references, and ranks likely organism matches.

The active implementation is in `src/ardronics/`. Legacy R scripts and the original untrained U-Net scaffold are preserved in `.old/` and are not used.

**Start with the [complete User Guide](USER_GUIDE.md).** It covers the JSON configuration, lane-map CSV, ladders, controls, metadata, reference FASTA, and multi-enzyme consensus.

## Features

- Native-resolution analysis of bright or dark bands.
- Automatic lanes with regular-grid fallback and manual overrides.
- Sub-pixel band centers, integrated intensity, and relative intensity.
- Strict log10(bp) ladder interpolation and optional two-ladder smile correction.
- IUPAC-aware primer matching, with both primers included in the amplicon.
- Biopython enzyme-specific cut offsets and single/double digests.
- Numeric one-to-one band matching, Jaccard distances, and ambiguity flags.
- CSV results, QC overlays, UPGMA dendrograms, and an HTML report.

No trained vision model is required; the pipeline is deterministic and inspectable.

## Install

Python 3.11 or newer is required:

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -e '.[dev]'
```

Or use `python3 -m venv .venv` followed by `.venv/bin/pip install -e '.[dev]'`.

```bash
.venv/bin/ardronics --help
.venv/bin/python -m pytest -q
```

## Supplied Mycoplasma gel

`data/ardra_example_1.jpg` has 18 samples and 50 bp ladders in lanes 1 and 20:

```bash
.venv/bin/ardronics gel \
  --image data/ardra_example_1.jpg \
  --enzyme unknown --n-lanes 20 --crop 20,138,535,334 \
  --ladder 50bp --ladder-lane 1 --secondary-ladder-lane 20 \
  --sample-names M_gallisepticum,M_pneumoniae,M_glycophilum,M_bovirhinis,M_columborale,M_neurolyticum,M_flocculare,M_dispar,M_hyopneumoniae,M_hyorhinis,M_arginini,M_hyosynoviae,M_orale,M_columbinasale,M_agalactiae,M_bovis,M_lipofaciens,M_gallinarum \
  --output-dir results/example_1
```

Review `results/example_1/qc_overlay.png` before using the values. Dense-gel sample bands are marked without text; all values are in `bands.csv`. The image does not document its digest enzyme, so this example uses `unknown`; substitute the actual wet-lab enzyme before matching against a library.

### More supplied examples

Run the AluI figure and all four HSCR panels:

```bash
scripts/run_more_examples.sh
```

This writes `results/aluI_example/` and `results/hscr/{HaeIII,MspI,RsaI,AluI}/`. The AluI example uses 13 manually reviewed left-ladder positions; its faint right marker is not used for calibration. HSCR results retain pixel migration positions because the caption does not identify the visible 123 bp ladder multiples.

## Commands

### Read a gel

```bash
.venv/bin/ardronics gel \
  --image data/gel.jpg --enzyme HaeIII --n-lanes 4 \
  --ladder 100bp --ladder-lane 1 \
  --sample-names isolate_a,isolate_b,isolate_c \
  --output-dir results/gel
```

Presets are `50bp`, `100bp`, and `1kb`. For another ladder, use `--ladder-sizes 1500,1000,700,500,300,100` in top-to-bottom order (largest to smallest).

Correction options:

- `--crop x,y,width,height`: exclude labels and unrelated panels.
- `--polarity dark`: dark bands on a light background.
- `--lane-centers 50,150,250,350`: explicit lane x positions after cropping.
- `--ladder-y 100,140,190,250`: explicit ladder y positions after cropping.
- `--secondary-ladder-lane 4`: second ladder for smile correction.
- `--prominence 2.0`: weaker bands; higher values are more conservative.

The ladder band count must exactly match the supplied sizes. A mismatch stops analysis rather than silently assigning incorrect sizes.

Outputs are `bands.csv` (positions, sizes, intensities), `profiles.csv` (tidy matching input), `qc_overlay.png`, and `analysis.json` (settings/warnings).

### Generate a 16S reference library

```bash
.venv/bin/ardronics library \
  --fasta data/references.fa --enzymes HaeIII,HhaI \
  --double-digest HaeIII+HhaI --output results/library.csv
```

Default primers are 27F (`AGAGTTTGATCMTGGCTCAG`) and 1492R (`GGTTACCTTGTTACGACTT`). Override with `--forward-primer`, `--reverse-primer`, and `--primer-mismatches`.

`library.csv` retains full FASTA descriptions. `library_qc.csv` identifies records whose primers could not be matched.

To download a small public NCBI reference set for the supplied Mycoplasma and Acinetobacter examples and generate HaeIII, HhaI, MspI, RsaI, and AluI profiles:

```bash
.venv/bin/python scripts/build_organism_library.py
```

The script writes `data/references/example_16s.fasta`, accession-level provenance in `results/organism_library/reference_manifest.json`, and the digest library/QC CSVs. Use `--organism "Species name"` to subset it or `--skip-download` to rebuild from the saved FASTA.

The example builder enables `--allow-unclipped`: records lacking a primer-defined amplicon are digested as supplied and marked `ok_unclipped` in QC. Their internal fragments are usable, but terminal fragments are approximate. Strict primer-defined amplicons remain the default for normal `ardronics library` runs.

To compare the calibrated AluI example while excluding its uncalibrated marker:

```bash
.venv/bin/ardronics identify \
  --samples results/aluI_example/profiles.csv \
  --library results/organism_library/library.csv \
  --exclude-sample MW_right_uncalibrated \
  --tolerance 30 --output-dir results/aluI_identification
```

This is a mechanics demonstration, not a validated taxonomic result; the example's observed and theoretical AluI profiles have substantial disagreement.

An in-silico library replaces a known-organism control, but it does **not** replace a molecular-weight ladder: pixel positions cannot be compared reliably to bp profiles without sizing information.

### Match existing profiles

```bash
.venv/bin/ardronics identify \
  --samples results/gel/profiles.csv --library results/library.csv \
  --tolerance 20 --min-bp 100 --top-n 3 \
  --output-dir results/identification
```

Legacy wide CSVs with `Sample,enznames,1,2,...,25` are accepted and converted to numeric values.

### Complete pipeline

Edit `config.example.json` for your gel and FASTA, then run:

```bash
.venv/bin/ardronics run --config config.example.json
```

Explicit CLI options override JSON values. Open `results/run/report.html` for the report, identifications, library QC, gel overlay, and dendrograms.

## Interpretation and limitations

- Matches are hypotheses, not definitive identifications. Organisms sharing a profile are grouped and flagged ambiguous.
- Tune `--tolerance` using known controls; gel sizing is less precise than capillary electrophoresis.
- Inspect QC, especially after a regular-grid lane fallback.
- Dual ladders correct linear smile. Strong curvature, tilt, overlap, and saturation may require manual coordinates or further processing.
- Biological identification accuracy still needs validation on your own known controls.

## Development

```bash
.venv/bin/python -m pytest -q
.venv/bin/python -m compileall -q src
```

Tests cover cut coordinates, IUPAC primers, numeric conversion, matching, synthetic and inverted gels, dual-ladder correction, the supplied real gel, and a complete image-to-identification run.
