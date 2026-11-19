#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
ARDRONICS=.venv/bin/ardronics

# AluI: the left ladder has 13 usable bands. Top/bottom edge artifacts are
# excluded by the explicit coordinates. The faint right marker is not used
# for calibration or treated as an unknown organism.
"$ARDRONICS" gel \
  --image data/Vol_8_no2_Typi_Afaf_fig5.jpg \
  --enzyme AluI --n-lanes 14 --crop 88,42,339,285 \
  --ladder-lane 1 \
  --ladder-sizes 2072,1500,1000,900,800,700,600,500,400,300,250,200,100 \
  --ladder-y 59.8245,71.3201,82.2045,87.8058,94.2063,101.1093,109.5424,119.8382,132.6827,146.0004,163.7537,186.1212,216.2830 \
  --sample-names 41,42,44,45,46,47,48,51,52,53,54,55,MW_right_uncalibrated \
  --output-dir results/aluI_example

# HSCR: each panel has two markers and samples T0-T15. The caption specifies
# a 123 bp ladder but not the visible multiples, so bp calibration is withheld.
HSCR_IMAGE=data/Gel-image-of-ARDRA-analysis-of-bacterial-16S-rDNA-retrieved-from-the-feces-of-an-HSCR.png
HSCR_NAMES=MW_left,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,MW_right
for panel in a:HaeIII:22:180 b:MspI:228:178 c:RsaI:435:174 d:AluI:638:130; do
  IFS=: read -r label enzyme y height <<< "$panel"
  "$ARDRONICS" gel \
    --image "$HSCR_IMAGE" --enzyme "$enzyme" \
    --n-lanes 18 --crop "0,$y,476,$height" --polarity dark \
    --sample-names "$HSCR_NAMES" \
    --output-dir "results/hscr/$enzyme"
done
