#!/usr/bin/env python3
"""Download GenBank flat files for all accessions listed in ictv_table1_all.csv."""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import requests

DEFAULT_CSV = Path(__file__).with_name("ictv_table1_all.csv")
DEFAULT_OUTPUT_DIR = Path(__file__).with_name("subtype_files")
EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download GenBank files for all accession numbers in ictv_table1_all.csv."
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=DEFAULT_CSV,
        help=f"CSV file containing accession numbers (default: {DEFAULT_CSV})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory where .gb files will be saved (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--email",
        help="Optional email to send with NCBI E-utilities requests.",
    )
    parser.add_argument(
        "--api-key",
        help="Optional NCBI API key.",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.34,
        help="Delay in seconds between requests (default: 0.34).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Redownload files even if they already exist.",
    )
    return parser.parse_args(argv)


def load_accessions(csv_path: Path) -> list[str]:
    if not csv_path.exists():
        raise SystemExit(f"Input CSV not found: {csv_path}")

    accessions: list[str] = []
    seen: set[str] = set()

    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            field = row.get("Accession number(s)", "")
            for accession in (part.strip() for part in field.split(",")):
                if not accession:
                    continue
                if accession not in seen:
                    seen.add(accession)
                    accessions.append(accession)

    if not accessions:
        raise SystemExit(f"No accession numbers found in {csv_path}")

    return accessions


def fetch_genbank(accession: str, email: str | None, api_key: str | None) -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gbwithparts",
        "retmode": "text",
        "tool": "ictv-subtype-downloader",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    try:
        response = requests.get(EUTILS_URL, params=params, timeout=60)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise RuntimeError(f"request failed for {accession}: {exc}") from exc

    text = response.text.strip()
    if not text:
        raise RuntimeError(f"empty response for {accession}")
    if text.startswith("Error:") or "Failed to understand id" in text:
        raise RuntimeError(f"NCBI returned an error for {accession}: {text}")
    if not text.startswith("LOCUS"):
        raise RuntimeError(f"unexpected response for {accession}: {text[:160]}")
    return text + "\n"


def download_all(
    accessions: list[str],
    output_dir: Path,
    email: str | None,
    api_key: str | None,
    delay: float,
    overwrite: bool,
) -> tuple[int, int]:
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = 0
    skipped = 0

    for accession in accessions:
        output_path = output_dir / f"{accession}.gb"
        if output_path.exists() and not overwrite:
            skipped += 1
            continue

        record = fetch_genbank(accession, email, api_key)
        output_path.write_text(record, encoding="utf-8")
        downloaded += 1
        time.sleep(delay)

    return downloaded, skipped


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    accessions = load_accessions(args.input_csv)

    downloaded, skipped = download_all(
        accessions=accessions,
        output_dir=args.output_dir,
        email=args.email,
        api_key=args.api_key,
        delay=args.delay,
        overwrite=args.overwrite,
    )
    print(
        f"Processed {len(accessions)} accessions: downloaded {downloaded}, skipped {skipped}, "
        f"output_dir={args.output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
