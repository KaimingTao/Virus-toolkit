#!/usr/bin/env python3
"""
Search Crossref by keyword and export results to CSV files.

The command-line interface intentionally mirrors literature_scripts/pubmed_search.py:
    uv run crossref_search.py --term HIV --out hiv_crossref.csv --email you@example.com
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Dict, Iterator, List, Optional, Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


CROSSREF_API_URL = "https://api.crossref.org/works"
ROWS_PER_REQUEST = 1000
REQUEST_DELAY_SECONDS = 0.2
CACHE_DIR = Path(__file__).with_name("crossref_search_cache")
OUTPUT_DIR_NAME = "crossref_search"
CSV_COLUMNS = [
    "DOI",
    "Title",
    "Authors",
    "Journal/Book",
    "Publication Year",
    "Published Date",
    "Publisher",
    "Type",
    "URL",
    "Abstract",
]


class ChunkedCsvWriter:
    def __init__(
        self,
        base_path: Path,
        rows_per_file: int,
        fieldnames: List[str],
    ) -> None:
        if rows_per_file < 1:
            raise ValueError("rows_per_file must be at least 1.")

        self.base_path = base_path
        self.rows_per_file = rows_per_file
        self.fieldnames = fieldnames
        self.output_dir = base_path.with_name(OUTPUT_DIR_NAME)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self._file_handle = None
        self._writer: Optional[csv.DictWriter] = None
        self._rows_in_chunk = 0
        self._chunk_start_row = 1
        self.files_written = 0

    def _chunk_path(self, start_row: int, end_row: int) -> Path:
        return self.output_dir / f"{self.base_path.stem}_{start_row:06d}_{end_row:06d}.csv"

    def _open_chunk(self) -> None:
        end_row = self._chunk_start_row + self.rows_per_file - 1
        chunk_path = self._chunk_path(self._chunk_start_row, end_row)
        self._file_handle = chunk_path.open("w", newline="", encoding="utf-8-sig")
        self._writer = csv.DictWriter(self._file_handle, fieldnames=self.fieldnames)
        self._writer.writeheader()
        self.files_written += 1

    def write_row(self, row: Dict[str, str], row_number: int) -> None:
        if self._writer is None:
            self._chunk_start_row = row_number
            self._rows_in_chunk = 0
            self._open_chunk()
        elif self._rows_in_chunk >= self.rows_per_file:
            self.close()
            self._chunk_start_row = row_number
            self._rows_in_chunk = 0
            self._open_chunk()

        assert self._writer is not None
        self._writer.writerow(row)
        self._rows_in_chunk += 1

    def close(self) -> None:
        if self._file_handle is not None:
            self._file_handle.close()
            self._file_handle = None
            self._writer = None


def normalize_cache_params(params: Dict[str, str]) -> Dict[str, str]:
    return {key: params[key] for key in sorted(params) if key != "mailto"}


def build_cache_key(params: Dict[str, str]) -> str:
    payload = json.dumps(normalize_cache_params(params), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def load_cache_entry(params: Dict[str, str]) -> Optional[Dict[str, object]]:
    cache_path = CACHE_DIR / f"{build_cache_key(params)}.json"
    if not cache_path.exists():
        return None
    try:
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    if payload.get("params") != normalize_cache_params(params):
        return None
    return payload


def save_cache_entry(params: Dict[str, str], response: Dict[str, object]) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / f"{build_cache_key(params)}.json"
    cache_path.write_text(
        json.dumps({"params": normalize_cache_params(params), "response": response}, sort_keys=True),
        encoding="utf-8",
    )


def request_json(params: Dict[str, str], email: Optional[str]) -> Dict[str, object]:
    cached_entry = load_cache_entry(params)
    if isinstance(cached_entry, dict):
        cached_response = cached_entry.get("response")
        if isinstance(cached_response, dict):
            print(f"[cache] cursor={params.get('cursor', '*')} rows={params.get('rows')}")
            return cached_response

    query = urlencode(params)
    request = Request(f"{CROSSREF_API_URL}?{query}")
    if email:
        request.add_header("User-Agent", f"virus-toolkits-crossref-search/0.1 (mailto:{email})")
    print(f"[api] cursor={params.get('cursor', '*')} rows={params.get('rows')}")
    with urlopen(request, timeout=60) as response:
        payload = json.loads(response.read().decode("utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError("Unexpected Crossref response format.")
    save_cache_entry(params, payload)
    return payload


def first_text(values: object) -> str:
    if isinstance(values, list):
        for value in values:
            if isinstance(value, str) and value.strip():
                return value.strip()
        return ""
    if isinstance(values, str):
        return values.strip()
    return ""


def join_authors(authors: object) -> str:
    if not isinstance(authors, list):
        return ""
    names: List[str] = []
    for author in authors:
        if not isinstance(author, dict):
            continue
        given = str(author.get("given", "")).strip()
        family = str(author.get("family", "")).strip()
        literal = str(author.get("literal", "")).strip()
        if literal:
            names.append(literal)
            continue
        full_name = " ".join(part for part in [given, family] if part)
        if full_name:
            names.append(full_name)
    return "; ".join(names)


def extract_date_parts(value: object) -> tuple[str, str]:
    if not isinstance(value, dict):
        return "", ""
    date_parts = value.get("date-parts")
    if not isinstance(date_parts, list) or not date_parts:
        return "", ""
    first_part = date_parts[0]
    if not isinstance(first_part, list) or not first_part:
        return "", ""

    numeric_parts: List[int] = []
    for part in first_part[:3]:
        if isinstance(part, int):
            numeric_parts.append(part)
    if not numeric_parts:
        return "", ""

    year = str(numeric_parts[0])
    date_string = "-".join(f"{part:02d}" if index else str(part) for index, part in enumerate(numeric_parts))
    return year, date_string


def clean_abstract(value: object) -> str:
    if not isinstance(value, str):
        return ""
    return " ".join(value.replace("<jats:p>", " ").replace("</jats:p>", " ").split())


def item_to_row(item: Dict[str, object]) -> Dict[str, str]:
    year, published_date = extract_date_parts(item.get("published-print") or item.get("published-online") or item.get("issued"))
    return {
        "DOI": str(item.get("DOI", "")).strip(),
        "Title": first_text(item.get("title")),
        "Authors": join_authors(item.get("author")),
        "Journal/Book": first_text(item.get("container-title")),
        "Publication Year": year,
        "Published Date": published_date,
        "Publisher": str(item.get("publisher", "")).strip(),
        "Type": str(item.get("type", "")).strip(),
        "URL": str(item.get("URL", "")).strip(),
        "Abstract": clean_abstract(item.get("abstract")),
    }


def iter_crossref_rows(
    term: str,
    email: Optional[str],
    max_results: Optional[int],
) -> Iterator[Dict[str, str]]:
    cursor = "*"
    fetched = 0

    while True:
        remaining = None if max_results is None else max_results - fetched
        if remaining is not None and remaining <= 0:
            return

        rows = ROWS_PER_REQUEST if remaining is None else min(ROWS_PER_REQUEST, remaining)
        params = {
            "query.bibliographic": term,
            "rows": str(rows),
            "cursor": cursor,
        }
        if email:
            params["mailto"] = email

        response = request_json(params, email=email)
        message = response.get("message")
        if not isinstance(message, dict):
            raise RuntimeError("Unexpected Crossref response: missing message object.")

        items = message.get("items")
        if not isinstance(items, list):
            raise RuntimeError("Unexpected Crossref response: missing items list.")

        next_cursor = message.get("next-cursor")
        if not isinstance(next_cursor, str) or not next_cursor:
            raise RuntimeError("Unexpected Crossref response: missing next-cursor.")

        for item in items:
            if not isinstance(item, dict):
                continue
            yield item_to_row(item)
            fetched += 1
            if max_results is not None and fetched >= max_results:
                return

        if len(items) < rows:
            return

        cursor = next_cursor
        time.sleep(REQUEST_DELAY_SECONDS)


def write_csv(
    term: str,
    out_path: str,
    email: Optional[str],
    max_results: Optional[int],
    rows_per_file: int,
) -> tuple[int, int, Path]:
    base_path = Path(out_path)
    writer = ChunkedCsvWriter(base_path, rows_per_file, CSV_COLUMNS)

    written = 0
    try:
        for row_number, row in enumerate(iter_crossref_rows(term, email, max_results), start=1):
            writer.write_row(row, row_number)
            written += 1
    finally:
        writer.close()

    return written, writer.files_written, writer.output_dir


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Search Crossref and save results to CSV.")
    parser.add_argument("--term", default="HIV", help="Search query for Crossref.")
    parser.add_argument(
        "--out",
        default="hiv_crossref.csv",
        help="Base CSV file path. Results are written to ./crossref_search/*.csv.",
    )
    parser.add_argument("--email", help="Your email address for the Crossref polite pool.")
    parser.add_argument(
        "--cache-dir",
        default=str(CACHE_DIR),
        help="Directory for API cache files (default: ./crossref_search_cache next to the script).",
    )
    parser.add_argument("--max-results", type=int, help="Optional cap on number of records to fetch.")
    parser.add_argument(
        "--rows-per-file",
        type=int,
        default=10000,
        help="Number of data rows per output CSV file (default: 10000).",
    )
    return parser.parse_args(argv)


def configure_cache_dir(cache_dir: str) -> None:
    global CACHE_DIR
    CACHE_DIR = Path(cache_dir)


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    configure_cache_dir(args.cache_dir)
    try:
        written, files_written, output_dir = write_csv(
            args.term,
            args.out,
            args.email,
            args.max_results,
            args.rows_per_file,
        )
    except (HTTPError, URLError, RuntimeError) as exc:
        print(f"Failed: {exc}", file=sys.stderr)
        return 1

    if written == 0:
        print(f"No Crossref records were written for term {args.term!r}.")
        return 0

    print(f"Wrote {written} Crossref records across {files_written} file(s) in {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
