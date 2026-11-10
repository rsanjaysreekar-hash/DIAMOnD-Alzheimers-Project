#!/usr/bin/env python3
"""
Download the DIAMOnD supplementary datasets and prepare Alzheimer's DIAMOnD inputs.

Run from the repository root:
    python download_data.py

Outputs:
    data/raw/original_interactome.tsv
    data/raw/original_disease_genes.tsv
    data/processed/PPI.txt
    data/processed/seed_genes.txt

The script uses the same PLOS supplementary-file sources referenced in the
project notebook. It keeps only PPI interaction columns 1 and 2, and extracts
the "Alzheimer disease" seed-gene identifiers that are present in the network.
"""

from pathlib import Path
import re
import sys
import urllib.request

BASE_DIR = Path(__file__).resolve().parent
RAW_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"

SOURCES = {
    "original_interactome.tsv": (
        "https://journals.plos.org/ploscompbiol/article/file?"
        "id=10.1371/journal.pcbi.1004120.s003&type=supplementary"
    ),
    "original_disease_genes.tsv": (
        "https://journals.plos.org/ploscompbiol/article/file?"
        "id=10.1371/journal.pcbi.1004120.s004&type=supplementary"
    ),
}


def download_file(filename: str, url: str) -> Path:
    destination = RAW_DIR / filename
    request = urllib.request.Request(
        url, headers={"User-Agent": "Mozilla/5.0 (compatible; DIAMOnD-project/1.0)"}
    )
    print(f"Downloading {filename} ...")
    with urllib.request.urlopen(request, timeout=90) as response:
        content = response.read()

    if len(content) < 1000 or b"<html" in content[:1000].lower():
        raise RuntimeError(
            f"{filename} does not look like the expected supplementary data. "
            "The source may be temporarily unavailable; check the PLOS URL."
        )

    destination.write_bytes(content)
    print(f"Saved {destination} ({len(content):,} bytes)")
    return destination


def read_interactome(path: Path):
    """Read the first two tab-separated gene-ID columns, ignoring comments."""
    interactions = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            gene_a, gene_b = fields[0].strip(), fields[1].strip()
            # Skip a possible non-comment header row.
            if gene_a.lower().startswith("gene_id") or gene_a.startswith("#gene"):
                continue
            if gene_a and gene_b:
                interactions.append((gene_a, gene_b))

    # Preserve order while removing duplicate undirected interactions.
    unique = []
    seen = set()
    for a, b in interactions:
        key = tuple(sorted((a, b)))
        if key not in seen:
            seen.add(key)
            unique.append((a, b))
    if not unique:
        raise RuntimeError("No PPI interactions were parsed from the interactome file.")
    return unique


def read_alzheimer_seed_ids(path: Path):
    """Parse disease rows ending in slash-separated numeric gene IDs."""
    diseases = {}
    pattern = re.compile(r"^(.*?)\s+([0-9]+(?:/[0-9]+)*)$")
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            match = pattern.match(line)
            if match:
                disease_name = match.group(1).strip().lower()
                gene_ids = set(match.group(2).split("/"))
                diseases[disease_name] = gene_ids

    if "alzheimer disease" not in diseases:
        raise RuntimeError(
            "Could not find an 'Alzheimer disease' row in the disease-gene file. "
            "Inspect the downloaded file format before proceeding."
        )
    return diseases["alzheimer disease"]


def main():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    for filename, url in SOURCES.items():
        target = RAW_DIR / filename
        if not target.exists() or target.stat().st_size < 1000:
            download_file(filename, url)
        else:
            print(f"Using existing raw file: {target}")

    interactions = read_interactome(RAW_DIR / "original_interactome.tsv")
    network_genes = {gene for pair in interactions for gene in pair}

    all_seeds = read_alzheimer_seed_ids(RAW_DIR / "original_disease_genes.tsv")
    mapped_seeds = sorted(all_seeds.intersection(network_genes))

    if not mapped_seeds:
        raise RuntimeError("No Alzheimer's seed genes were found in the PPI network.")

    ppi_path = PROCESSED_DIR / "PPI.txt"
    with ppi_path.open("w", encoding="utf-8") as handle:
        for gene_a, gene_b in interactions:
            handle.write(f"{gene_a}\t{gene_b}\n")

    seed_path = PROCESSED_DIR / "seed_genes.txt"
    seed_path.write_text("\n".join(mapped_seeds) + "\n", encoding="utf-8")

    print("\nData preparation complete.")
    print(f"PPI interactions: {len(interactions):,}")
    print(f"Unique network proteins: {len(network_genes):,}")
    print(f"Alzheimer's seed genes in source: {len(all_seeds)}")
    print(f"Alzheimer's seed genes present in network: {len(mapped_seeds)}")
    print(f"Prepared network: {ppi_path}")
    print(f"Prepared seed list: {seed_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
