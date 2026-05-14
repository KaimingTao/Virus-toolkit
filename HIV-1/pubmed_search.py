#!/usr/bin/env python3
"""
Search PubMed for a query term and export the results to a CSV file.

Example:
    python pubmed_search.py --term "machine learning cancer" --out results.csv --email you@example.com
"""

import argparse
import csv
import sys
import time
from typing import Dict, Iterable, List, Optional

import requests


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DEFAULT_BATCH_SIZE = 200  # E-utilities generally allow up to 200 IDs per summary call.


def _request_json(endpoint: str, params: Dict[str, str]) -> Dict:
    url = f"{EUTILS_BASE}/{endpoint}"
    response = requests.get(url, params=params, timeout=20)
    response.raise_for_status()
    return response.json()


def fetch_pmids(term: str, email: Optional[str], api_key: Optional[str], max_results: Optional[int]) -> List[str]:
    # Initial call to get total count.
    common_params = {"db": "pubmed", "retmode": "json", "term": term}
    if email:
        common_params["email"] = email
    if api_key:
        common_params["api_key"] = api_key

    initial = _request_json("esearch.fcgi", {**common_params, "retmax": 0})
    total_count = int(initial["esearchresult"]["count"])

    if max_results is not None:
        total_count = min(total_count, max_results)

    pmids: List[str] = []
    for start in range(0, total_count, DEFAULT_BATCH_SIZE):
        remaining = total_count - start
        batch_size = min(DEFAULT_BATCH_SIZE, remaining)
        payload = {
            **common_params,
            "retstart": start,
            "retmax": batch_size,
        }
        search = _request_json("esearch.fcgi", payload)
        pmids.extend(search["esearchresult"].get("idlist", []))
        time.sleep(0.34)  # Be gentle with NCBI rate limits (~3 requests/second without key).

    return pmids


def fetch_summaries(pmids: Iterable[str], email: Optional[str], api_key: Optional[str]) -> List[Dict[str, str]]:
    pmid_list = list(pmids)
    summaries: List[Dict[str, str]] = []
    common_params = {"db": "pubmed", "retmode": "json"}
    if email:
        common_params["email"] = email
    if api_key:
        common_params["api_key"] = api_key

    for start in range(0, len(pmid_list), DEFAULT_BATCH_SIZE):
        batch = pmid_list[start : start + DEFAULT_BATCH_SIZE]
        payload = {**common_params, "id": ",".join(batch), "version": "2.0"}
        data = _request_json("esummary.fcgi", payload)
        result = data.get("result", {})
        for uid in result.get("uids", []):
            item = result.get(uid, {})
            authors = item.get("authors") or []
            author_names = "; ".join(a.get("name", "") for a in authors if a.get("name"))
            journal = item.get("fulljournalname") or item.get("source") or ""
            pubdate = item.get("pubdate") or ""
            article_ids = item.get("articleids") or []
            doi = ""
            for aid in article_ids:
                if aid.get("idtype") == "doi":
                    doi = aid.get("value", "")
                    break
            summaries.append(
                {
                    "pmid": uid,
                    "title": item.get("title") or "",
                    "journal": journal,
                    "pubdate": pubdate,
                    "authors": author_names,
                    "doi": doi,
                }
            )
        time.sleep(0.34)

    return summaries


def write_csv(rows: List[Dict[str, str]], path: str) -> None:
    fieldnames = ["pmid", "title", "journal", "pubdate", "authors", "doi"]
    with open(path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Search PubMed and save results to CSV.")
    parser.add_argument("--term", required=True, help="Search query for PubMed.")
    parser.add_argument("--out", required=True, help="Output CSV file path.")
    parser.add_argument("--email", help="Your email (recommended for NCBI E-utilities).")
    parser.add_argument("--api-key", help="NCBI API key to increase rate limits.")
    parser.add_argument("--max-results", type=int, help="Optional cap on number of articles to fetch.")
    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    pmids = fetch_pmids(args.term, args.email, args.api_key, args.max_results)
    if not pmids:
        print("No results found.", file=sys.stderr)
        return 1

    summaries = fetch_summaries(pmids, args.email, args.api_key)
    write_csv(summaries, args.out)
    print(f"Wrote {len(summaries)} records to {args.out}")
    return 0


if __name__ == "__main__":
    main(sys.argv[1:])
