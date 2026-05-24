#!/usr/bin/env python3
"""Write compact PMC identifier and package URL CSVs from PubMed CSV exports."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET


IDCONV_URL = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"
OA_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"
IDCONV_BATCH_SIZE = 200
REQUEST_DELAY_SECONDS = 0.34
CACHE_DIR = Path(__file__).with_name("fetch_pmc_files_cache")
PMID_COLUMN = "PMID"
PMCID_COLUMN = "PMCID"
PMC_PDF_URL_COLUMN = "PMC PDF URL"
PMC_TGZ_URL_COLUMN = "PMC TGZ URL"
SKIP_SUFFIX = "_pmc.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan a folder of PubMed CSV files and write one compact sibling "
            "_pmc.csv per input file with PMC identifiers and package URLs."
        )
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Folder containing CSV files exported by literature_scripts/pubmed_search.py.",
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Contact email passed to NCBI APIs.",
    )
    parser.add_argument(
        "--tool",
        default="virus_toolkit_pmc_fetch",
        help="Tool name passed to NCBI APIs.",
    )
    parser.add_argument(
        "--cache-dir",
        default=str(CACHE_DIR),
        help="Directory for API cache files (default: ./fetch_pmc_files_cache next to the script).",
    )
    return parser.parse_args()


def iter_input_csvs(input_dir: Path) -> Iterator[Path]:
    for path in sorted(input_dir.glob("*.csv")):
        if not path.name.endswith(SKIP_SUFFIX):
            yield path


def chunked(values: Sequence[str], size: int) -> Iterator[List[str]]:
    for start in range(0, len(values), size):
        yield list(values[start:start + size])


def normalize_cache_params(params: Dict[str, str]) -> Dict[str, str]:
    return {key: params[key] for key in sorted(params) if key not in {"email"}}


def build_cache_key(endpoint: str, params: Dict[str, str]) -> str:
    normalized = {"endpoint": endpoint, "params": normalize_cache_params(params)}
    serialized = json.dumps(normalized, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def build_cache_path(endpoint: str, cache_key: str) -> Path:
    return CACHE_DIR / endpoint / f"{cache_key}.json"


def configure_cache_dir(cache_dir: str) -> None:
    global CACHE_DIR
    CACHE_DIR = Path(cache_dir)


def load_cache_entry(endpoint: str, params: Dict[str, str], cache_key: str) -> Optional[Dict[str, object]]:
    cache_path = build_cache_path(endpoint, cache_key)
    if not cache_path.exists():
        return None
    try:
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    if payload.get("endpoint") != endpoint or payload.get("params") != normalize_cache_params(params):
        return None
    return payload


def save_cache_entry(endpoint: str, params: Dict[str, str], cache_key: str, response: Dict[str, object]) -> None:
    cache_path = build_cache_path(endpoint, cache_key)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = cache_path.with_suffix(".tmp")
    temp_path.write_text(
        json.dumps(
            {
                "endpoint": endpoint,
                "params": normalize_cache_params(params),
                "response": response,
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    temp_path.replace(cache_path)


def request_json(
    endpoint: str,
    url: str,
    params: Dict[str, str],
    retries: int = 3,
) -> Tuple[Dict[str, object], bool]:
    cache_key = build_cache_key(endpoint, params)
    cached_entry = load_cache_entry(endpoint, params, cache_key)
    if isinstance(cached_entry, dict):
        cached_response = cached_entry.get("response")
        if isinstance(cached_response, dict):
            print(f"[cache] {endpoint}")
            return cached_response, True

    query = urlencode(params)
    full_url = f"{url}?{query}"
    last_error: Optional[Exception] = None
    print(f"[api] {endpoint}")
    for attempt in range(1, retries + 1):
        try:
            with urlopen(full_url, timeout=60) as response:
                payload = json.loads(response.read().decode("utf-8"))
            if not isinstance(payload, dict):
                raise RuntimeError("Unexpected JSON response format.")
            save_cache_entry(endpoint, params, cache_key, payload)
            return payload, False
        except (HTTPError, URLError, json.JSONDecodeError, RuntimeError) as exc:
            last_error = exc
            if attempt == retries:
                break
            print(f"[retry] {endpoint} attempt {attempt}/{retries} error={exc}")
            time.sleep(attempt)
    raise RuntimeError(f"Request failed for {endpoint}: {last_error}") from last_error


def request_xml(
    endpoint: str,
    url: str,
    params: Dict[str, str],
    retries: int = 3,
) -> Tuple[ET.Element, bool]:
    cache_key = build_cache_key(endpoint, params)
    cached_entry = load_cache_entry(endpoint, params, cache_key)
    if isinstance(cached_entry, dict):
        cached_response = cached_entry.get("response")
        if isinstance(cached_response, dict):
            xml_text = cached_response.get("xml")
            if isinstance(xml_text, str):
                print(f"[cache] {endpoint}")
                return ET.fromstring(xml_text), True

    query = urlencode(params)
    full_url = f"{url}?{query}"
    last_error: Optional[Exception] = None
    print(f"[api] {endpoint}")
    for attempt in range(1, retries + 1):
        try:
            with urlopen(full_url, timeout=60) as response:
                xml_text = response.read().decode("utf-8")
            root = ET.fromstring(xml_text)
            save_cache_entry(endpoint, params, cache_key, {"xml": xml_text})
            return root, False
        except (HTTPError, URLError, ET.ParseError) as exc:
            last_error = exc
            if attempt == retries:
                break
            print(f"[retry] {endpoint} attempt {attempt}/{retries} error={exc}")
            time.sleep(attempt)
    raise RuntimeError(f"Request failed for {endpoint}: {last_error}") from last_error


def sanitize_pmid(value: str) -> str:
    return value.strip()


def fetch_pmcid_map(pmids: Sequence[str], email: str, tool: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    unique_pmids = [pmid for pmid in dict.fromkeys(pmids) if pmid]
    total_batches = (len(unique_pmids) + IDCONV_BATCH_SIZE - 1) // IDCONV_BATCH_SIZE

    for batch_index, batch in enumerate(chunked(unique_pmids, IDCONV_BATCH_SIZE), start=1):
        print(
            f"[pmcid] batch {batch_index}/{total_batches} "
            f"pmids={len(batch)} resolved_so_far={len(mapping)}"
        )
        response, from_cache = request_json(
            "idconv",
            IDCONV_URL,
            {
                "ids": ",".join(batch),
                "idtype": "pmid",
                "format": "json",
                "email": email,
                "tool": tool,
            },
        )
        records = response.get("records", [])
        if not isinstance(records, list):
            raise RuntimeError("Unexpected ID converter response: missing records list.")

        for record in records:
            if not isinstance(record, dict):
                continue
            pmid = str(record.get("pmid", "")).strip()
            pmcid = str(record.get("pmcid", "")).strip()
            if pmid and pmcid:
                mapping[pmid] = pmcid

        if not from_cache:
            time.sleep(REQUEST_DELAY_SECONDS)

    return mapping


def fetch_oa_links(pmcid: str) -> Tuple[Dict[str, str], bool]:
    root, from_cache = request_xml("oa", OA_URL, {"id": pmcid})
    links = {
        PMC_PDF_URL_COLUMN: "",
        PMC_TGZ_URL_COLUMN: "",
    }
    for record in root.findall(".//record"):
        if record.attrib.get("id") != pmcid:
            continue
        for link in record.findall("link"):
            href = link.attrib.get("href", "").strip()
            fmt = link.attrib.get("format", "").strip().lower()
            if not href:
                continue
            if fmt == "pdf":
                links[PMC_PDF_URL_COLUMN] = href
            elif fmt == "tgz":
                links[PMC_TGZ_URL_COLUMN] = href
        break
    return links, from_cache


def process_csv(csv_path: Path, email: str, tool: str) -> Path:
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        if PMID_COLUMN not in fieldnames:
            raise RuntimeError(f"{csv_path.name} is missing the required {PMID_COLUMN!r} column.")
        rows = list(reader)

    print(f"[rows] {csv_path.name} total_rows={len(rows)}")
    pmids = [sanitize_pmid(row.get(PMID_COLUMN, "")) for row in rows]
    pmcid_map = fetch_pmcid_map(pmids, email=email, tool=tool)
    print(f"[pmcid] {csv_path.name} matched_pmcids={len(pmcid_map)}")

    output_fieldnames = [
        PMID_COLUMN,
        PMCID_COLUMN,
        PMC_PDF_URL_COLUMN,
        PMC_TGZ_URL_COLUMN,
    ]

    oa_links_by_pmcid: Dict[str, Dict[str, str]] = {}
    output_rows: List[Dict[str, str]] = []
    total_rows = len(rows)
    for row_index, row in enumerate(rows, start=1):
        pmid = sanitize_pmid(row.get(PMID_COLUMN, ""))
        pmcid = pmcid_map.get(pmid, "")
        output_row = {
            PMID_COLUMN: pmid,
            PMCID_COLUMN: pmcid,
            PMC_PDF_URL_COLUMN: "",
            PMC_TGZ_URL_COLUMN: "",
        }

        if not pmcid:
            output_rows.append(output_row)
            continue

        if pmcid not in oa_links_by_pmcid:
            print(f"[oa] {csv_path.name} row {row_index}/{total_rows} pmcid={pmcid}")
            oa_links, from_cache = fetch_oa_links(pmcid)
            oa_links_by_pmcid[pmcid] = oa_links
            if not from_cache:
                time.sleep(REQUEST_DELAY_SECONDS)

        output_row.update(oa_links_by_pmcid[pmcid])
        output_rows.append(output_row)

    output_path = csv_path.with_name(f"{csv_path.stem}_pmc.csv")
    with output_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"[write] {output_path.name} rows={len(output_rows)}")
    return output_path


def main() -> int:
    args = parse_args()
    configure_cache_dir(args.cache_dir)
    input_dir = args.input_dir.resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"Input folder not found: {input_dir}")

    csv_paths = list(iter_input_csvs(input_dir))
    if not csv_paths:
        raise SystemExit(f"No input CSV files found in {input_dir}")

    total_output_files = 0
    for csv_path in csv_paths:
        print(f"[csv] {csv_path.name}")
        output_path = process_csv(csv_path, email=args.email, tool=args.tool)
        total_output_files += 1
        print(f"[done] {output_path.name}")

    print(f"[summary] csv_outputs={total_output_files}")
    print(f"[summary] added_columns={PMCID_COLUMN}, {PMC_PDF_URL_COLUMN}, {PMC_TGZ_URL_COLUMN}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
