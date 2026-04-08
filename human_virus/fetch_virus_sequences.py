#!/usr/bin/env python3
"""Fetch reference sequences from accessions listed in a virus CSV."""

from __future__ import annotations

import argparse
import csv
import re
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

BASE_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
INVALID_CHARS = re.compile(r"[^A-Za-z0-9._-]+")


def sanitize_name(value: str, fallback: str) -> str:
    cleaned = INVALID_CHARS.sub("_", value.strip()).strip("._-")
    return cleaned or fallback


def fetch_fasta(
    accession: str,
    email: str | None,
    api_key: str | None,
    retries: int,
    timeout: float,
) -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    url = f"{BASE_EFETCH_URL}?{urlencode(params)}"

    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urlopen(url, timeout=timeout) as response:
                body = response.read().decode("utf-8", errors="replace").strip()
            if not body or not body.startswith(">"):
                raise ValueError(f"Unexpected response for {accession}: {body[:120]!r}")
            return body + "\n"
        except (HTTPError, URLError, TimeoutError, ValueError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(min(2**attempt, 10))

    raise RuntimeError(f"Failed to fetch {accession}: {last_error}")


def iter_accessions(raw: str, delimiter: str) -> list[str]:
    return [part.strip() for part in raw.split(delimiter) if part.strip()]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Read a virus CSV and fetch FASTA by reference accession."
    )
    parser.add_argument("--csv", required=True, help="Path to input CSV file")
    parser.add_argument(
        "--output-dir",
        default="virus_sequences",
        help="Directory to write virus folders into (default: virus_sequences)",
    )
    parser.add_argument(
        "--name-column",
        default="name",
        help="CSV column used for virus folder names (default: name)",
    )
    parser.add_argument(
        "--accession-column",
        default="reference_accessions",
        help="CSV column containing accessions (default: reference_accessions)",
    )
    parser.add_argument(
        "--accession-delimiter",
        default=";",
        help="Delimiter between accessions in a cell (default: ';')",
    )
    parser.add_argument(
        "--email",
        default=None,
        help="Optional email to pass to NCBI E-utilities",
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="Optional NCBI API key for higher rate limits",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of retries per accession (default: 3)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="HTTP timeout in seconds (default: 30)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing FASTA files if present",
    )
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=0.34,
        help="Delay between requests (default: 0.34)",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()

    csv_path = Path(args.csv)
    out_root = Path(args.output_dir)

    if not csv_path.exists():
        print(f"CSV not found: {csv_path}", file=sys.stderr)
        return 1

    out_root.mkdir(parents=True, exist_ok=True)

    fetched = 0
    skipped = 0
    failed = 0

    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            print("CSV appears to have no header.", file=sys.stderr)
            return 1
        if args.name_column not in reader.fieldnames:
            print(
                f"Missing name column '{args.name_column}'. Available: {reader.fieldnames}",
                file=sys.stderr,
            )
            return 1
        if args.accession_column not in reader.fieldnames:
            print(
                "Missing accession column "
                f"'{args.accession_column}'. Available: {reader.fieldnames}",
                file=sys.stderr,
            )
            return 1

        for row_idx, row in enumerate(reader, start=2):
            virus_name = (row.get(args.name_column) or "").strip()
            accession_raw = (row.get(args.accession_column) or "").strip()
            if not accession_raw:
                continue

            folder_name = sanitize_name(virus_name, f"virus_{row_idx}")
            virus_dir = out_root / folder_name
            virus_dir.mkdir(parents=True, exist_ok=True)

            accessions = iter_accessions(accession_raw, args.accession_delimiter)
            for accession in accessions:
                out_file = virus_dir / f"{sanitize_name(accession, accession)}.fasta"
                if out_file.exists() and not args.overwrite:
                    skipped += 1
                    continue

                try:
                    fasta = fetch_fasta(
                        accession=accession,
                        email=args.email,
                        api_key=args.api_key,
                        retries=args.retries,
                        timeout=args.timeout,
                    )
                except Exception as exc:  # noqa: BLE001
                    failed += 1
                    print(
                        f"[row {row_idx}] Failed {accession} for '{virus_name}': {exc}",
                        file=sys.stderr,
                    )
                    continue

                out_file.write_text(fasta, encoding="utf-8")
                fetched += 1
                if args.sleep_seconds > 0:
                    time.sleep(args.sleep_seconds)

    print(
        f"Done. fetched={fetched} skipped={skipped} failed={failed} "
        f"output_dir={out_root}"
    )
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
