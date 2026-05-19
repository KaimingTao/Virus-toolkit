#!/usr/bin/env python3
"""
Search PubMed with a query term and export results to a CSV file.

The output columns are aligned to PubMed's CSV export header as closely as the
E-utilities API allows:
    PMID,Title,Authors,Citation,First Author,Journal/Book,Publication Year,
    Create Date,PMCID,NIHMS ID,DOI

Example:
    uv run pubmed_search.py --term HIV --out hiv_pubmed.csv
"""

from __future__ import annotations

import argparse
import csv
from datetime import date, timedelta
import hashlib
from http.client import IncompleteRead
import json
from pathlib import Path
import re
import sys
import time
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FETCH_BATCH_SIZE = 100
REQUEST_DELAY_SECONDS = 0.34
PUBMED_QUERY_LIMIT = 9999
DATE_FIELD = "PDAT"
EARLIEST_PUBMED_DATE = date(1800, 1, 1)
CACHE_DIR = Path(__file__).with_name("pubmed_search_cache")
OUTPUT_DIR_NAME = "pubmed_search"
CSV_COLUMNS = [
    "PMID",
    "Title",
    "Abstract",
    "Authors",
    "Citation",
    "First Author",
    "Journal/Book",
    "Publication Year",
    "Create Date",
    "PMCID",
    "NIHMS ID",
    "DOI",
]
PMID_INDEX_COLUMNS = ["PMID", "File"]
CHUNK_NAME_RE = re.compile(r"^(?P<stem>.+)_(?P<start>\d{6})_(?P<end>\d{6})\.csv$")


class ChunkedCsvWriter:
    def __init__(
        self,
        base_path: Path,
        rows_per_file: int,
        fieldnames: List[str],
        start_row: int = 1,
        existing_files: int = 0,
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
        self._chunk_start_row = start_row
        self._rows_in_chunk = 0
        self._temp_chunk_path: Optional[Path] = None
        self.files_written = existing_files

    def _build_chunk_path(self, start_row: int, end_row: int) -> Path:
        return self.output_dir / f"{self.base_path.stem}_{start_row:06d}_{end_row:06d}.csv"

    def _open_next_chunk(self) -> None:
        self._temp_chunk_path = self.output_dir / f"{self.base_path.stem}_{self._chunk_start_row:06d}.part"
        self._file_handle = self._temp_chunk_path.open("w", newline="", encoding="utf-8-sig")
        self._writer = csv.DictWriter(self._file_handle, fieldnames=self.fieldnames)
        self._writer.writeheader()
        self.files_written += 1

    def write_row(self, row: Dict[str, str], row_number: int) -> None:
        if self._writer is None or self._rows_in_chunk >= self.rows_per_file:
            self.close()
            self._chunk_start_row = row_number
            self._rows_in_chunk = 0
            self._open_next_chunk()

        assert self._writer is not None
        self._writer.writerow(row)
        self._rows_in_chunk += 1

    def close(self) -> None:
        if self._file_handle is not None:
            self._file_handle.close()
            if self._temp_chunk_path is not None and self._rows_in_chunk > 0:
                end_row = self._chunk_start_row + self._rows_in_chunk - 1
                self._temp_chunk_path.replace(self._build_chunk_path(self._chunk_start_row, end_row))
            self._file_handle = None
            self._writer = None
            self._temp_chunk_path = None


def normalize_cache_params(params: Dict[str, str]) -> Dict[str, str]:
    return {key: params[key] for key in sorted(params) if key not in {"email", "api_key"}}


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

    stored_endpoint = payload.get("endpoint")
    stored_params = payload.get("params")
    expected_params = normalize_cache_params(params)
    if stored_endpoint != endpoint or stored_params != expected_params:
        return None

    return payload


def save_cache_entry(endpoint: str, cache_key: str, entry: Dict[str, object]) -> None:
    cache_path = build_cache_path(endpoint, cache_key)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = cache_path.with_suffix(".tmp")
    temp_path.write_text(json.dumps(entry, ensure_ascii=True, sort_keys=True), encoding="utf-8")
    temp_path.replace(cache_path)


def build_params(email: Optional[str], api_key: Optional[str]) -> Dict[str, str]:
    params: Dict[str, str] = {}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    return params


def request_xml(endpoint: str, params: Dict[str, str], retries: int = 3) -> ET.Element:
    cache_key = build_cache_key(endpoint, params)
    query_summary = summarize_query(params.get("term", ""), limit=80)
    cached_entry = load_cache_entry(endpoint, params, cache_key)
    if isinstance(cached_entry, dict):
        cached_xml = cached_entry.get("xml")
        if isinstance(cached_xml, str):
            if endpoint == "efetch.fcgi":
                batch_start = int(params.get("retstart", "0")) + 1
                batch_size = int(params.get("retmax", "0"))
                batch_end = batch_start + batch_size - 1 if batch_size > 0 else batch_start
                print(
                    f"[cache] {endpoint} rows {batch_start}-{batch_end}"
                )
            elif query_summary:
                print(f"[cache] {endpoint} query: {query_summary}")
            else:
                print(f"[cache] {endpoint}")
            return ET.fromstring(cached_xml)

    query = urlencode(params)
    url = f"{EUTILS_BASE}/{endpoint}?{query}"
    if endpoint == "efetch.fcgi":
        batch_start = int(params.get("retstart", "0")) + 1
        batch_size = int(params.get("retmax", "0"))
        batch_end = batch_start + batch_size - 1 if batch_size > 0 else batch_start
        print(f"[api] {endpoint} rows {batch_start}-{batch_end}")
    elif query_summary:
        print(f"[api] {endpoint} query: {query_summary}")
    else:
        print(f"[api] {endpoint}")

    last_error: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            with urlopen(url, timeout=60) as response:
                xml_bytes = response.read()
                root = ET.fromstring(xml_bytes)
                save_cache_entry(endpoint, cache_key, {
                    "endpoint": endpoint,
                    "params": normalize_cache_params(params),
                    "xml": xml_bytes.decode("utf-8"),
                })
                return root
        except IncompleteRead as exc:
            last_error = exc
            partial = exc.partial
            if partial:
                try:
                    root = ET.fromstring(partial)
                    save_cache_entry(endpoint, cache_key, {
                        "endpoint": endpoint,
                        "params": normalize_cache_params(params),
                        "xml": partial.decode("utf-8"),
                    })
                    return root
                except ET.ParseError:
                    pass
            if attempt == retries:
                break
            time.sleep(attempt)
        except (HTTPError, URLError, ET.ParseError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(attempt)

    raise RuntimeError(f"Request failed for {endpoint}: {last_error}") from last_error


def text_at(node: Optional[ET.Element], path: str, default: str = "") -> str:
    if node is None:
        return default
    found = node.find(path)
    if found is None or found.text is None:
        return default
    return found.text.strip()


def iter_texts(node: Optional[ET.Element], path: str) -> Iterator[str]:
    if node is None:
        return
    for found in node.findall(path):
        if found.text:
            value = found.text.strip()
            if value:
                yield value


def join_nonempty(parts: Iterable[str], sep: str = " ") -> str:
    return sep.join(part for part in parts if part)


def summarize_query(term: str, limit: int = 120) -> str:
    normalized = " ".join(term.split())
    if len(normalized) <= limit:
        return normalized
    return normalized[: limit - 3] + "..."


def chunk_sort_key(path: Path) -> Tuple[int, int, str]:
    match = CHUNK_NAME_RE.match(path.name)
    if match is None:
        return (sys.maxsize, sys.maxsize, path.name)
    return (int(match.group("start")), int(match.group("end")), path.name)


def iter_chunk_paths(output_dir: Path, stem: str) -> List[Path]:
    paths: List[Path] = []
    for path in output_dir.glob(f"{stem}_*.csv"):
        match = CHUNK_NAME_RE.match(path.name)
        if match is None or match.group("stem") != stem:
            continue
        paths.append(path)
    paths.sort(key=chunk_sort_key)
    return paths


def pmid_index_path(output_dir: Path, stem: str) -> Path:
    return output_dir / f"{stem}_pmid_index.csv"


def build_pmid_index(output_dir: Path, stem: str) -> tuple[Set[str], int, int, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    chunk_paths = iter_chunk_paths(output_dir, stem)
    index_path = pmid_index_path(output_dir, stem)
    seen_pmids: Set[str] = set()
    existing_rows = 0

    print(f"Creating PMID index cache from split files in {output_dir}")
    with index_path.open("w", newline="", encoding="utf-8") as index_handle:
        writer = csv.writer(index_handle)
        writer.writerow(PMID_INDEX_COLUMNS)

        for chunk_number, chunk_path in enumerate(chunk_paths, start=1):
            print(f"  Indexing chunk {chunk_number}/{len(chunk_paths)}: {chunk_path.name}")
            with chunk_path.open("r", newline="", encoding="utf-8-sig") as chunk_handle:
                reader = csv.DictReader(chunk_handle)
                if reader.fieldnames is None or "PMID" not in reader.fieldnames:
                    raise RuntimeError(f"Chunk file is missing PMID column: {chunk_path}")
                for row in reader:
                    pmid = (row.get("PMID") or "").strip()
                    if not pmid:
                        continue
                    existing_rows += 1
                    if pmid in seen_pmids:
                        continue
                    seen_pmids.add(pmid)
                    writer.writerow([pmid, chunk_path.name])

    print(
        f"Indexed {len(seen_pmids)} unique PMID(s) across {existing_rows} existing row(s) "
        f"from {len(chunk_paths)} file(s)"
    )
    return seen_pmids, existing_rows, len(chunk_paths), index_path


def load_pmid_index(output_dir: Path, stem: str) -> tuple[Set[str], int, int, Path]:
    index_path = pmid_index_path(output_dir, stem)
    seen_pmids: Set[str] = set()
    referenced_files: Set[str] = set()

    print(f"Loading PMID index cache from {index_path}")
    with index_path.open("r", newline="", encoding="utf-8-sig") as index_handle:
        reader = csv.DictReader(index_handle)
        if reader.fieldnames != PMID_INDEX_COLUMNS:
            raise RuntimeError(
                f"PMID index cache has unexpected columns in {index_path}: {reader.fieldnames}"
            )
        for row in reader:
            pmid = (row.get("PMID") or "").strip()
            file_name = (row.get("File") or "").strip()
            if not pmid:
                continue
            seen_pmids.add(pmid)
            if file_name:
                referenced_files.add(file_name)

    print(
        f"Loaded {len(seen_pmids)} PMID(s) already present in {len(referenced_files)} CSV file(s)"
    )
    return seen_pmids, len(seen_pmids), len(referenced_files), index_path


def load_or_build_pmid_index(output_dir: Path, stem: str) -> tuple[Set[str], int, int, Path]:
    index_path = pmid_index_path(output_dir, stem)
    if index_path.exists():
        return load_pmid_index(output_dir, stem)
    return build_pmid_index(output_dir, stem)


def collect_pubmed_history(
    term: str,
    email: Optional[str],
    api_key: Optional[str],
    max_results: Optional[int],
) -> tuple[int, str, str]:
    print(f"Preparing PubMed history for query: {summarize_query(term)}")
    params = {
        "db": "pubmed",
        "term": term,
        "retmode": "xml",
        "usehistory": "y",
        "retmax": "0",
        **build_params(email, api_key),
    }
    root = request_xml("esearch.fcgi", params)
    count = int(text_at(root, "./Count", "0"))
    webenv = text_at(root, "./WebEnv")
    query_key = text_at(root, "./QueryKey")
    if not webenv or not query_key:
        raise RuntimeError("PubMed search did not return WebEnv/QueryKey.")
    if max_results is not None:
        count = min(count, max_results)
    time.sleep(REQUEST_DELAY_SECONDS)
    return count, webenv, query_key


def count_pubmed_hits(term: str, email: Optional[str], api_key: Optional[str]) -> int:
    print(f"Counting PubMed hits for query: {summarize_query(term)}")
    params = {
        "db": "pubmed",
        "term": term,
        "rettype": "count",
        "retmode": "xml",
        "retmax": "0",
        **build_params(email, api_key),
    }
    root = request_xml("esearch.fcgi", params)
    count = int(text_at(root, "./Count", "0"))
    time.sleep(REQUEST_DELAY_SECONDS)
    return count


def build_dated_term(term: str, start_date: date, end_date: date) -> str:
    date_clause = (
        f'"{start_date:%Y/%m/%d}"[{DATE_FIELD}] : "{end_date:%Y/%m/%d}"[{DATE_FIELD}]'
    )
    return f"({term}) AND ({date_clause})"


def split_date_range(start_date: date, end_date: date) -> Tuple[Tuple[date, date], Tuple[date, date]]:
    total_days = (end_date - start_date).days
    midpoint = start_date + timedelta(days=total_days // 2)
    left = (start_date, midpoint)
    right = (midpoint + timedelta(days=1), end_date)
    return left, right


def build_query_segments_by_halving(
    term: str,
    email: Optional[str],
    api_key: Optional[str],
    start_date: date,
    end_date: date,
    limit: int = PUBMED_QUERY_LIMIT,
) -> List[Tuple[str, int]]:
    dated_term = build_dated_term(term, start_date, end_date)
    print(f"Checking date range {start_date:%Y-%m-%d} to {end_date:%Y-%m-%d}")
    count = count_pubmed_hits(dated_term, email, api_key)
    if count == 0:
        print(f"  Range {start_date:%Y-%m-%d} to {end_date:%Y-%m-%d}: 0 hits")
        return []
    if count <= limit:
        print(
            f"  Range {start_date:%Y-%m-%d} to {end_date:%Y-%m-%d}: "
            f"{count} hits, keeping as one segment"
        )
        return [(dated_term, count)]
    if start_date >= end_date:
        raise RuntimeError(
            f"PubMed returned {count} results for single day {start_date:%Y-%m-%d}, "
            "which still exceeds the supported query limit."
        )

    print(
        f"  Range {start_date:%Y-%m-%d} to {end_date:%Y-%m-%d}: "
        f"{count} hits, splitting further"
    )
    left_range, right_range = split_date_range(start_date, end_date)
    segments = build_query_segments_by_halving(term, email, api_key, *left_range, limit=limit)
    segments.extend(build_query_segments_by_halving(term, email, api_key, *right_range, limit=limit))
    return segments


def iter_half_year_date_ranges(start_date: date, end_date: date) -> Iterator[Tuple[date, date]]:
    for year in range(start_date.year, end_date.year + 1):
        for range_start, range_end in (
            (date(year, 1, 1), date(year, 6, 30)),
            (date(year, 7, 1), date(year, 12, 31)),
        ):
            clipped_start = max(start_date, range_start)
            clipped_end = min(end_date, range_end)
            if clipped_start <= clipped_end:
                yield clipped_start, clipped_end


def iter_quarterly_date_ranges(start_date: date, end_date: date) -> Iterator[Tuple[date, date]]:
    for year in range(start_date.year, end_date.year + 1):
        for range_start, range_end in (
            (date(year, 1, 1), date(year, 3, 31)),
            (date(year, 4, 1), date(year, 6, 30)),
            (date(year, 7, 1), date(year, 9, 30)),
            (date(year, 10, 1), date(year, 12, 31)),
        ):
            clipped_start = max(start_date, range_start)
            clipped_end = min(end_date, range_end)
            if clipped_start <= clipped_end:
                yield clipped_start, clipped_end


def iter_monthly_date_ranges(start_date: date, end_date: date) -> Iterator[Tuple[date, date]]:
    current = date(start_date.year, start_date.month, 1)
    while current <= end_date:
        if current.month == 12:
            next_month = date(current.year + 1, 1, 1)
        else:
            next_month = date(current.year, current.month + 1, 1)
        month_end = next_month - timedelta(days=1)
        clipped_start = max(start_date, current)
        clipped_end = min(end_date, month_end)
        if clipped_start <= clipped_end:
            yield clipped_start, clipped_end
        current = next_month


def build_query_segments_from_ranges(
    term: str,
    email: Optional[str],
    api_key: Optional[str],
    date_ranges: List[Tuple[date, date]],
    label: str,
    limit: int = PUBMED_QUERY_LIMIT,
) -> List[Tuple[str, int]]:
    segments: List[Tuple[str, int]] = []
    print(f"Building {len(date_ranges)} fixed {label} date range(s)")
    for range_start, range_end in date_ranges:
        dated_term = build_dated_term(term, range_start, range_end)
        print(f"Checking date range {range_start:%Y-%m-%d} to {range_end:%Y-%m-%d}")
        count = count_pubmed_hits(dated_term, email, api_key)
        if count == 0:
            print(f"  Range {range_start:%Y-%m-%d} to {range_end:%Y-%m-%d}: 0 hits")
            continue
        if count <= limit:
            print(
                f"  Range {range_start:%Y-%m-%d} to {range_end:%Y-%m-%d}: "
                f"{count} hits, keeping as one segment"
            )
            segments.append((dated_term, count))
            continue
        print(
            f"  Range {range_start:%Y-%m-%d} to {range_end:%Y-%m-%d}: "
            f"{count} hits, splitting to next stage"
        )
        if label == "half-year":
            sub_ranges = list(iter_quarterly_date_ranges(range_start, range_end))
            segments.extend(
                build_query_segments_from_ranges(
                    term, email, api_key, sub_ranges, "quarterly", limit=limit
                )
            )
        elif label == "quarterly":
            sub_ranges = list(iter_monthly_date_ranges(range_start, range_end))
            segments.extend(
                build_query_segments_from_ranges(
                    term, email, api_key, sub_ranges, "monthly", limit=limit
                )
            )
        else:
            segments.extend(
                build_query_segments_by_halving(
                    term, email, api_key, range_start, range_end, limit=limit
                )
            )
    return segments


def build_range_list_query_segments(
    term: str,
    email: Optional[str],
    api_key: Optional[str],
    start_date: date,
    end_date: date,
    limit: int = PUBMED_QUERY_LIMIT,
) -> List[Tuple[str, int]]:
    date_ranges = list(iter_half_year_date_ranges(start_date, end_date))
    return build_query_segments_from_ranges(
        term, email, api_key, date_ranges, "half-year", limit=limit
    )


def fetch_batch(
    webenv: str,
    query_key: str,
    start: int,
    batch_size: int,
    email: Optional[str],
    api_key: Optional[str],
) -> ET.Element:
    params = {
        "db": "pubmed",
        "query_key": query_key,
        "WebEnv": webenv,
        "retstart": str(start),
        "retmax": str(batch_size),
        "retmode": "xml",
        **build_params(email, api_key),
    }
    return request_xml("efetch.fcgi", params)


def extract_author(author: ET.Element) -> str:
    collective = text_at(author, "./CollectiveName")
    if collective:
        return collective

    last_name = text_at(author, "./LastName")
    fore_name = text_at(author, "./ForeName")
    initials = text_at(author, "./Initials")

    if last_name and fore_name:
        return f"{last_name} {fore_name}"
    if last_name and initials:
        return f"{last_name} {initials}"
    return last_name or fore_name or initials


def extract_authors(article: ET.Element) -> List[str]:
    authors: List[str] = []
    for author in article.findall(".//AuthorList/Author"):
        value = extract_author(author)
        if value:
            authors.append(value)
    return authors


def extract_title(article: ET.Element) -> str:
    title_node = article.find("./ArticleTitle")
    if title_node is None:
        return ""
    return "".join(title_node.itertext()).strip()


def extract_journal(article: ET.Element) -> str:
    return (
        text_at(article, "./Journal/Title")
        or text_at(article, "./Book/BookTitle")
        or text_at(article, "./Journal/ISOAbbreviation")
    )


def extract_publication_year(article: ET.Element, pubmed_data: ET.Element) -> str:
    for path in (
        "./Journal/JournalIssue/PubDate/Year",
        "./ArticleDate/Year",
        "./Book/PubDate/Year",
    ):
        value = text_at(article, path)
        if value:
            return value

    medline_date = text_at(article, "./Journal/JournalIssue/PubDate/MedlineDate")
    if medline_date:
        digits = "".join(ch for ch in medline_date if ch.isdigit())
        if len(digits) >= 4:
            return digits[:4]

    for pub_status in pubmed_data.findall("./History/PubMedPubDate"):
        if pub_status.attrib.get("PubStatus") == "pubmed":
            value = text_at(pub_status, "./Year")
            if value:
                return value

    return ""


def extract_create_date(pubmed_data: ET.Element) -> str:
    for pub_status in pubmed_data.findall("./History/PubMedPubDate"):
        if pub_status.attrib.get("PubStatus") == "pubmed":
            year = text_at(pub_status, "./Year")
            month = text_at(pub_status, "./Month").zfill(2)
            day = text_at(pub_status, "./Day").zfill(2)
            if year and month and day:
                return f"{year}/{month}/{day}"
    return ""


def extract_article_ids(pubmed_data: ET.Element) -> Dict[str, str]:
    ids: Dict[str, str] = {}
    for article_id in pubmed_data.findall("./ArticleIdList/ArticleId"):
        id_type = article_id.attrib.get("IdType", "").lower()
        value = (article_id.text or "").strip()
        if id_type and value and id_type not in ids:
            ids[id_type] = value
    return ids


def extract_abstract(article: ET.Element) -> str:
    abstract_texts: List[str] = []
    for abstract_node in article.findall("./Abstract/AbstractText"):
        text = " ".join("".join(abstract_node.itertext()).split())
        if not text:
            continue
        label = (abstract_node.attrib.get("Label") or "").strip()
        if label:
            abstract_texts.append(f"{label}: {text}")
        else:
            abstract_texts.append(text)
    return " ".join(abstract_texts)


def build_citation(article: ET.Element, pubmed_data: ET.Element, doi: str) -> str:
    journal = extract_journal(article)
    year = extract_publication_year(article, pubmed_data)
    volume = text_at(article, "./Journal/JournalIssue/Volume")
    issue = text_at(article, "./Journal/JournalIssue/Issue")
    pages = text_at(article, "./Pagination/MedlinePgn")

    pieces: List[str] = []
    if journal:
        pieces.append(f"{journal}.")

    date_bits = [year]
    month = text_at(article, "./Journal/JournalIssue/PubDate/Month")
    day = text_at(article, "./Journal/JournalIssue/PubDate/Day")
    if month:
        date_bits.append(month)
    if day:
        date_bits.append(day)
    date_text = " ".join(bit for bit in date_bits if bit)
    if date_text:
        pieces.append(date_text + ";")

    volume_issue = volume
    if issue:
        volume_issue = f"{volume}({issue})" if volume else f"({issue})"
    if volume_issue:
        pieces.append(volume_issue)

    if pages:
        if volume_issue:
            pieces[-1] = pieces[-1] + f":{pages}."
        else:
            pieces.append(f"{pages}.")

    citation = " ".join(pieces).strip()
    if doi:
        citation = f"{citation} doi: {doi}.".strip()
    return citation


def article_to_row(article_node: ET.Element) -> Dict[str, str]:
    medline = article_node.find("./MedlineCitation")
    article = medline.find("./Article") if medline is not None else None
    pubmed_data = article_node.find("./PubmedData")

    if medline is None or article is None or pubmed_data is None:
        raise RuntimeError("Unexpected PubMed article structure in efetch response.")

    authors = extract_authors(article)
    author_text = ", ".join(authors)
    first_author = authors[0] if authors else ""
    ids = extract_article_ids(pubmed_data)
    doi = ids.get("doi", "")

    return {
        "PMID": text_at(medline, "./PMID"),
        "Title": extract_title(article),
        "Abstract": extract_abstract(article),
        "Authors": author_text,
        "Citation": build_citation(article, pubmed_data, doi),
        "First Author": first_author,
        "Journal/Book": extract_journal(article),
        "Publication Year": extract_publication_year(article, pubmed_data),
        "Create Date": extract_create_date(pubmed_data),
        "PMCID": ids.get("pmc", ""),
        "NIHMS ID": ids.get("mid", ""),
        "DOI": doi,
    }


def write_csv(
    term: str,
    out_path: str,
    email: Optional[str],
    api_key: Optional[str],
    max_results: Optional[int],
    rows_per_file: int,
) -> tuple[int, int, Path]:
    base_path = Path(out_path)
    output_dir = base_path.with_name(OUTPUT_DIR_NAME)
    existing_pmids, existing_rows, existing_files, _ = load_or_build_pmid_index(
        output_dir, base_path.stem
    )
    print(f"{existing_rows} PMID(s) are already present in the split CSV cache.")

    total_count = count_pubmed_hits(term, email, api_key)
    print(f"Total PubMed hits for {term!r}: {total_count}")
    if total_count == 0:
        with open(out_path, "w", newline="", encoding="utf-8-sig") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=CSV_COLUMNS)
            writer.writeheader()
        return 0, existing_files, base_path

    today = date.today()
    segments = build_range_list_query_segments(term, email, api_key, EARLIEST_PUBMED_DATE, today)
    print(f"Split into {len(segments)} PubMed date segment(s) to stay under the API limit.")

    total_unique_rows = existing_rows
    new_rows_written = 0
    skipped_existing = 0
    hit_max_results = False
    seen_pmids: Set[str] = set(existing_pmids)
    writer = ChunkedCsvWriter(
        base_path,
        rows_per_file,
        CSV_COLUMNS,
        start_row=existing_rows + 1,
        existing_files=existing_files,
    )
    try:
        for segment_index, (segment_term, segment_count) in enumerate(segments, start=1):
            remaining = None if max_results is None else max_results - total_unique_rows
            if remaining is not None and remaining <= 0:
                break

            fetch_count, webenv, query_key = collect_pubmed_history(
                segment_term,
                email,
                api_key,
                remaining if remaining is not None else None,
            )
            fetch_total = min(segment_count, fetch_count)
            print(
                f"Segment {segment_index}/{len(segments)}: fetching {fetch_total} record(s) "
                f"for query {segment_term!r}"
            )

            for start in range(0, fetch_total, FETCH_BATCH_SIZE):
                batch_size = min(FETCH_BATCH_SIZE, fetch_total - start)
                batch_number = start // FETCH_BATCH_SIZE + 1
                batch_end = min(start + batch_size, fetch_total)
                print(
                    f"  Batch {batch_number}: loading segment records "
                    f"{start + 1}-{batch_end} of {fetch_total}"
                )
                root = fetch_batch(webenv, query_key, start, batch_size, email, api_key)
                batch_written = 0
                for article_node in root.findall("./PubmedArticle"):
                    row = article_to_row(article_node)
                    pmid = row["PMID"]
                    if pmid in seen_pmids:
                        skipped_existing += 1
                        continue
                    seen_pmids.add(pmid)
                    total_unique_rows += 1
                    new_rows_written += 1
                    writer.write_row(row, total_unique_rows)
                    batch_written += 1
                    if max_results is not None and total_unique_rows >= max_results:
                        hit_max_results = True
                        break
                segment_progress = min(start + batch_size, fetch_total)
                print(
                    f"  Batch {batch_number}: "
                    f"{segment_progress}/{fetch_total} segment records processed, "
                    f"{total_unique_rows}/{total_count if max_results is None else min(total_count, max_results)} "
                    f"total unique rows, {new_rows_written} new rows written "
                    f"(batch wrote {batch_written}, skipped {skipped_existing})"
                )
                if hit_max_results:
                    break
                time.sleep(REQUEST_DELAY_SECONDS)
            if hit_max_results:
                break
    finally:
        writer.close()

    build_pmid_index(writer.output_dir, base_path.stem)
    return new_rows_written, writer.files_written, writer.output_dir


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Search PubMed and save results to CSV.")
    parser.add_argument("--term", default="HIV", help="Search query for PubMed.")
    parser.add_argument(
        "--out",
        default="hiv_pubmed.csv",
        help="Base CSV file path. Results are written to ./pubmed_search/*.csv.",
    )
    parser.add_argument("--email", help="Your email address for NCBI E-utilities.")
    parser.add_argument("--api-key", help="NCBI API key to increase rate limits.")
    parser.add_argument(
        "--cache-dir",
        default=str(CACHE_DIR),
        help="Directory for API cache files (default: ./pubmed_search_cache next to the script).",
    )
    parser.add_argument("--max-results", type=int, help="Optional cap on number of articles to fetch.")
    parser.add_argument(
        "--rows-per-file",
        type=int,
        default=10000,
        help="Number of data rows per output CSV file (default: 10000).",
    )
    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    configure_cache_dir(args.cache_dir)
    try:
        written, files_written, output_dir = write_csv(
            args.term,
            args.out,
            args.email,
            args.api_key,
            args.max_results,
            args.rows_per_file,
        )
    except Exception as exc:
        print(f"Failed: {exc}", file=sys.stderr)
        return 1

    if written == 0:
        print(f"No new PubMed records were written for term {args.term!r}.")
        return 0

    print(f"Wrote {written} PubMed records across {files_written} file(s) in {output_dir}")
    print("Note: the CSV columns match PubMed's export header, but citation formatting may not be byte-for-byte identical to the website export.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
