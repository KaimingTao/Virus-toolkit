#!/usr/bin/env python3
"""Download ICTV Table 1 page HTML and print its visible update date."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import requests
from bs4 import BeautifulSoup

DEFAULT_URL = "https://ictv.global/sg_wiki/flaviviridae/hepacivirus/table1"
DEFAULT_OUTPUT = "ictv_table1.html"
HEADING_PATTERN = re.compile(r"Table 1\s+[–-]\s+Confirmed HCV genotypes/subtypes\s+\(([^)]+)\)", re.I)


def fetch_page(url: str) -> str:
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.text
    except requests.RequestException as exc:
        raise SystemExit(f"Failed to download {url}: {exc}") from exc


def extract_update_date(page_html: str) -> str:
    soup = BeautifulSoup(page_html, "html.parser")
    heading = soup.find("h2", string=re.compile(r"Table 1", re.I))
    if heading is None:
        raise SystemExit("Could not find the Table 1 heading in the downloaded page.")

    match = HEADING_PATTERN.search(heading.get_text(" ", strip=True))
    if match is None:
        raise SystemExit("Could not extract the update date from the Table 1 heading.")
    return match.group(1)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download ICTV Table 1 and print the update date shown on the page."
    )
    parser.add_argument("--url", default=DEFAULT_URL, help=f"Source URL (default: {DEFAULT_URL})")
    parser.add_argument(
        "--output-html",
        default=DEFAULT_OUTPUT,
        help=f"Where to save the downloaded HTML (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    page_html = fetch_page(args.url)
    Path(args.output_html).write_text(page_html, encoding="utf-8")
    update_date = extract_update_date(page_html)
    print(update_date)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
