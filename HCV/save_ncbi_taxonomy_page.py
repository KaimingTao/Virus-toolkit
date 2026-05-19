#!/usr/bin/env python3
"""Download and save an NCBI taxonomy browser page."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import requests

DEFAULT_URL = (
    "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/"
    "wwwtax.cgi?command=show&mode=tree&id=3052230&lvl=3"
)
DEFAULT_OUTPUT = "ncbi_taxonomy_3052230_lvl3.html"


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download and save an NCBI taxonomy browser page."
    )
    parser.add_argument("--url", default=DEFAULT_URL, help=f"Source URL (default: {DEFAULT_URL})")
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help=f"Output HTML path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Request timeout in seconds (default: 30).",
    )
    return parser.parse_args(argv)


def fetch_page(url: str, timeout: float) -> str:
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return response.text
    except requests.RequestException as exc:
        raise SystemExit(f"Failed to download {url}: {exc}") from exc


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    page_html = fetch_page(args.url, args.timeout)
    output_path = Path(args.output)
    output_path.write_text(page_html, encoding="utf-8")
    print(f"Saved {args.url} to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
