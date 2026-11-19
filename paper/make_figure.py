"""Generate the original ARDRONICS workflow figure for the manuscript."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = Path(__file__).parent / "figures"
OUT.mkdir(parents=True, exist_ok=True)
fig, ax = plt.subplots(figsize=(12, 6.5))
ax.set_xlim(0, 12)
ax.set_ylim(0, 7)
ax.axis("off")


def box(x, y, w, h, title, lines, color):
    patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.08", facecolor=color, edgecolor="#263746", linewidth=1.2)
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h - 0.28, title, ha="center", va="top", fontsize=11, fontweight="bold", color="#152431")
    ax.text(x + w / 2, y + h - 0.68, lines, ha="center", va="top", fontsize=9.2, linespacing=1.45, color="#263746")


def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=15, linewidth=1.4, color="#486581"))


box(0.3, 4.25, 2.6, 1.8, "Gel inputs", "Original image\nLane-map CSV\nLadder sizes and locations", "#dbeafe")
box(0.3, 1.0, 2.6, 1.8, "Reference inputs", "16S FASTA\nPCR primers\nRestriction enzymes", "#e0e7ff")
box(3.65, 4.25, 2.8, 1.8, "Image analysis", "Lane traces and band peaks\nDensitometry and QC overlay\nSingle/dual-ladder calibration", "#bfdbfe")
box(3.65, 1.0, 2.8, 1.8, "In-silico library", "IUPAC-aware amplicon clipping\nEnzyme-specific cut positions\nReference QC and provenance", "#c7d2fe")
box(7.25, 2.65, 2.15, 1.8, "Profile matching", "One-to-one band\nassignment\nJaccard distances\nMulti-enzyme consensus", "#bbf7d0")
box(10.0, 2.65, 1.75, 1.8, "Outputs", "CSV profiles\nQC overlays and\ndendrograms\nHTML/control reports", "#fef3c7")
arrow(2.9, 5.15, 3.65, 5.15)
arrow(2.9, 1.9, 3.65, 1.9)
arrow(6.45, 5.0, 7.55, 4.45)
arrow(6.45, 2.05, 7.55, 2.65)
arrow(9.4, 3.55, 10.0, 3.55)
ax.text(6, 6.65, "ARDRONICS: reproducible gel-to-16S ARDRA analysis", ha="center", va="center", fontsize=16, fontweight="bold", color="#152431")
ax.text(6, 0.35, "Explicit metadata, ambiguity reporting, and positive-control validation connect the complete workflow.", ha="center", fontsize=10, color="#486581")
fig.savefig(OUT / "workflow.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(OUT / "workflow.pdf", bbox_inches="tight", facecolor="white")
