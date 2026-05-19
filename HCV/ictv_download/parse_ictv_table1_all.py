#!/usr/bin/env python3
"""Download ICTV Table 1 and export all genotype/subtype rows to CSV."""

from __future__ import annotations

import argparse
import csv
import html
import re
import sys
from pathlib import Path

import requests
from bs4 import BeautifulSoup

DEFAULT_URL = "https://ictv.global/sg_wiki/flaviviridae/hepacivirus/table1"
DEFAULT_OUTPUT = "ictv_table1_all.csv"
GENOTYPE_HEADER_RE = re.compile(r"^Genotype\s+(\d+)$", re.I)
SUBTYPE_RE = re.compile(r"^(\d+[a-z]{1,3})\b")
COLUMN_SPLIT_RE = re.compile(r"\s{2,}")
TAG_RE = re.compile(r"<[^>]+>")
SUP_RE = re.compile(r"<sup[^>]*>.*?</sup>", re.S | re.I)
ACCESSION_TOKEN_RE = re.compile(r"[A-Z]{1,4}\d{5,8}")


def fetch_html(url: str) -> str:
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.text
    except requests.RequestException as exc:
        raise SystemExit(f"Failed to fetch {url}: {exc}") from exc


def normalize_text(text: str) -> str:
    text = text.replace("\xa0", " ").strip()
    text = re.sub(r"\s+([,;)])", r"\1", text)
    text = re.sub(r"(\()\s+", r"\1", text)
    return text


def paragraph_text(paragraph) -> str:
    raw_html = paragraph.decode_contents()
    raw_html = SUP_RE.sub("", raw_html)
    text = TAG_RE.sub("", raw_html)
    return normalize_text(html.unescape(text))


def normalize_accessions(text: str) -> str:
    matches = ACCESSION_TOKEN_RE.findall(text)
    if len(matches) >= 2 and "," not in text:
        return ", ".join(matches)
    return text


def parse_row(text: str, current_genotype: str) -> dict[str, str] | None:
    reference_index = text.find("(")
    if reference_index == -1:
        return None

    left_side = text[:reference_index].rstrip()
    references = text[reference_index:].strip()
    columns = [part.strip() for part in COLUMN_SPLIT_RE.split(left_side) if part.strip()]
    if len(columns) != 3:
        return None

    subtype = columns[0]
    if not SUBTYPE_RE.fullmatch(subtype):
        return None

    return {
        "Genotype": current_genotype,
        "Subtype": subtype,
        "Locus/Isolate(s)": columns[1],
        "Accession number(s)": normalize_accessions(columns[2]),
        "Reference(s)": references,
    }


def extract_rows(page_html: str) -> list[dict[str, str]]:
    soup = BeautifulSoup(page_html, "html.parser")
    heading = soup.find("h2", string=re.compile(r"Table 1", re.I))
    if heading is None:
        raise SystemExit("Could not find ICTV Table 1 in the page.")

    current_genotype = ""
    rows: list[dict[str, str]] = []

    for paragraph in heading.find_all_next("p"):
        classes = paragraph.get("class", [])
        if "EndNoteBibliography" in classes:
            break

        text = paragraph_text(paragraph)
        if not text or "Locus/Isolate(s)" in text:
            continue

        header_match = GENOTYPE_HEADER_RE.fullmatch(text)
        if header_match:
            current_genotype = header_match.group(1)
            continue

        if not current_genotype or not SUBTYPE_RE.match(text):
            continue

        row = parse_row(text, current_genotype)
        if row is not None:
            rows.append(row)

    if not rows:
        raise SystemExit("No table rows were extracted from ICTV Table 1.")

    return rows


def write_csv(rows: list[dict[str, str]], output_path: Path) -> None:
    fieldnames = [
        "Genotype",
        "Subtype",
        "Locus/Isolate(s)",
        "Accession number(s)",
        "Reference(s)",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract all rows from ICTV Hepacivirus Table 1 into a CSV file."
    )
    parser.add_argument("--url", default=DEFAULT_URL, help=f"Source URL (default: {DEFAULT_URL})")
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help=f"CSV output path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    page_html = fetch_html(args.url)
    rows = extract_rows(page_html)
    write_csv(rows, Path(args.output))
    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
