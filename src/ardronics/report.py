"""Self-contained, human-readable run report."""

from html import escape
from pathlib import Path

import pandas as pd


def write_html_report(output_dir: str | Path, identifications: pd.DataFrame, warnings: list[str]) -> Path:
    output_dir = Path(output_dir)
    warning_html = "".join(f"<li>{escape(item)}</li>" for item in warnings) or "<li>None</li>"
    table = identifications.to_html(index=False, float_format=lambda value: f"{value:.3f}")
    dendrograms = sorted(output_dir.glob("dendrogram_*.png"))
    figures = "".join(
        f'<h3>{path.stem}</h3><img src="{path.name}" alt="{path.stem}">' for path in dendrograms
    )
    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>ARDRONICS report</title>
<style>body{{font:15px system-ui;max-width:1200px;margin:2rem auto;padding:0 1rem}}
table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #ccc;padding:.4rem;text-align:left}}
th{{background:#eee}}img{{max-width:100%;height:auto}}.warning{{color:#8a4b00}}</style></head>
<body><h1>ARDRONICS identification report</h1>
<p>Matches are hypotheses based on restriction profiles. Review ambiguous calls and the gel QC overlay.</p>
<h2>Gel QC</h2><img src="gel/qc_overlay.png" alt="Detected gel lanes and bands">
<h2 class="warning">Warnings</h2><ul>{warning_html}</ul>
<h2>Closest library profiles</h2>{table}{figures}</body></html>"""
    path = output_dir / "report.html"
    path.write_text(html)
    return path
