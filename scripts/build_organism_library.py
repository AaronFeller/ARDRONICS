#!/usr/bin/env python3
"""Download public 16S references from NCBI and build an ARDRA library."""

import argparse
import json
import time
import urllib.parse
import urllib.request
from io import StringIO
from pathlib import Path

from Bio import SeqIO

from ardronics.digest import generate_library

ORGANISMS = [
    "Mycoplasma gallisepticum", "Mycoplasma pneumoniae", "Mycoplasma glycophilum",
    "Mycoplasma bovirhinis", "Mycoplasma columborale", "Mycoplasma neurolyticum",
    "Mycoplasma flocculare", "Mycoplasma dispar", "Mycoplasma hyopneumoniae",
    "Mycoplasma hyorhinis", "Mycoplasma arginini", "Mycoplasma hyosynoviae",
    "Mycoplasma orale", "Mycoplasma columbinasale", "Mycoplasma agalactiae",
    "Mycoplasma bovis", "Mycoplasma lipofaciens", "Mycoplasma gallinarum",
    "Acinetobacter baumannii", "Acinetobacter nosocomialis",
    "Acinetobacter pittii", "Acinetobacter calcoaceticus",
]
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def request(endpoint, parameters):
    parameters = {**parameters, "tool": "ardronics"}
    url = BASE + endpoint + "?" + urllib.parse.urlencode(parameters)
    req = urllib.request.Request(url, headers={"User-Agent": "ARDRONICS/0.1"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as response:
                return response.read().decode()
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("unreachable")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--organism", action="append", help="Repeat to subset the default organism list")
    parser.add_argument("--records-per-organism", type=int, default=2)
    parser.add_argument("--fasta", default="data/references/example_16s.fasta")
    parser.add_argument("--output-dir", default="results/organism_library")
    parser.add_argument("--skip-download", action="store_true")
    args = parser.parse_args()
    fasta = Path(args.fasta)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    organisms = args.organism or ORGANISMS
    manifest = []

    if not args.skip_download:
        fasta.parent.mkdir(parents=True, exist_ok=True)
        all_records = []
        for organism in organisms:
            base_query = f'"{organism}"[Organism] AND 16S[Title] AND 1200:1800[SLEN]'
            query = base_query + " AND srcdb_refseq[PROP]"
            reference_source = "RefSeq"
            search = json.loads(request("esearch.fcgi", {
                "db": "nuccore", "term": query, "retmax": args.records_per_organism,
                "retmode": "json",
            }))
            ids = search["esearchresult"]["idlist"]
            time.sleep(0.4)
            if not ids:
                query = base_query
                reference_source = "GenBank"
                search = json.loads(request("esearch.fcgi", {
                    "db": "nuccore", "term": query, "retmax": args.records_per_organism,
                    "retmode": "json",
                }))
                ids = search["esearchresult"]["idlist"]
                time.sleep(0.4)
            records = []
            if ids:
                downloaded = request("efetch.fcgi", {
                    "db": "nuccore", "id": ",".join(ids), "rettype": "fasta", "retmode": "text",
                })
                records = list(SeqIO.parse(StringIO(downloaded), "fasta"))
            for record in records:
                record.description = f"{organism} | {record.description}"
            all_records.extend(records)
            manifest.append({
                "organism": organism, "query": query, "reference_source": reference_source,
                "ncbi_ids": ids, "records": [r.id for r in records],
            })
            print(f"{organism}: {len(records)} record(s) from {reference_source}", flush=True)
            time.sleep(0.4)
        written = SeqIO.write(all_records, fasta, "fasta")
        if not written:
            raise ValueError("No FASTA records were downloaded")
        (output / "reference_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    library = generate_library(
        fasta, output / "library.csv",
        enzymes=["HaeIII", "HhaI", "MspI", "RsaI", "AluI"],
        allow_unclipped=True,
    )
    print(f"Wrote {len(library)} fragments for {library['sample'].nunique()} references to {output / 'library.csv'}")


if __name__ == "__main__":
    main()
