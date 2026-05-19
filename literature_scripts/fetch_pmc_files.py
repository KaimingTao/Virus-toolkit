#!/usr/bin/env python3
"""Write compact PMC identifier and package URL CSVs from PubMed CSV exports."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Dict, Iterator, List, Sequence
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET


IDCONV_URL = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"
OA_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"
IDCONV_BATCH_SIZE = 200
REQUEST_DELAY_SECONDS = 0.34
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
    return parser.parse_args()


def iter_input_csvs(input_dir: Path) -> Iterator[Path]:
    for path in sorted(input_dir.glob("*.csv")):
        if not path.name.endswith(SKIP_SUFFIX):
            yield path


def chunked(values: Sequence[str], size: int) -> Iterator[List[str]]:
    for start in range(0, len(values), size):
        yield list(values[start:start + size])


def request_json(url: str, params: Dict[str, str]) -> Dict[str, object]:
    query = urlencode(params)
    with urlopen(f"{url}?{query}", timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def request_xml(url: str, params: Dict[str, str]) -> ET.Element:
    query = urlencode(params)
    with urlopen(f"{url}?{query}", timeout=60) as response:
        return ET.fromstring(response.read())


def sanitize_pmid(value: str) -> str:
    return value.strip()


def fetch_pmcid_map(pmids: Sequence[str], email: str, tool: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    unique_pmids = [pmid for pmid in dict.fromkeys(pmids) if pmid]

    for batch in chunked(unique_pmids, IDCONV_BATCH_SIZE):
        response = request_json(
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

        time.sleep(REQUEST_DELAY_SECONDS)

    return mapping


def fetch_oa_links(pmcid: str) -> Dict[str, str]:
    root = request_xml(OA_URL, {"id": pmcid})
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
    return links


def process_csv(csv_path: Path, email: str, tool: str) -> Path:
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        if PMID_COLUMN not in fieldnames:
            raise RuntimeError(f"{csv_path.name} is missing the required {PMID_COLUMN!r} column.")
        rows = list(reader)

    pmids = [sanitize_pmid(row.get(PMID_COLUMN, "")) for row in rows]
    pmcid_map = fetch_pmcid_map(pmids, email=email, tool=tool)

    output_fieldnames = [
        PMID_COLUMN,
        PMCID_COLUMN,
        PMC_PDF_URL_COLUMN,
        PMC_TGZ_URL_COLUMN,
    ]

    oa_links_by_pmcid: Dict[str, Dict[str, str]] = {}
    output_rows: List[Dict[str, str]] = []
    for row in rows:
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
            oa_links_by_pmcid[pmcid] = fetch_oa_links(pmcid)
            time.sleep(REQUEST_DELAY_SECONDS)

        output_row.update(oa_links_by_pmcid[pmcid])
        output_rows.append(output_row)

    output_path = csv_path.with_name(f"{csv_path.stem}_pmc.csv")
    with output_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    return output_path


def main() -> int:
    args = parse_args()
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
