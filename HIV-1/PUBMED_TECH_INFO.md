# PubMed Export Technical Notes

## Script

The PubMed export script is [pubmed_search.py](/Users/kaimingtao/Projects/Virus-toolkit/HIV-1/pubmed_search.py).

Run it with:

```bash
uv run pubmed_search.py --term HIV --out hiv_pubmed.csv
```

## How the query works

The script uses the NCBI PubMed E-utilities API in two phases.

1. `esearch.fcgi`

It sends a search request with:

- `db=pubmed`
- `term=HIV`
- `usehistory=y`
- `retmax=0`
- `retmode=xml`

This step does not download all papers. It asks PubMed to store the full result set on the NCBI side and returns:

- total hit count
- `WebEnv`
- `QueryKey`

2. `efetch.fcgi`

The script then fetches records in batches using:

- `db=pubmed`
- `query_key=<QueryKey>`
- `WebEnv=<WebEnv>`
- `retstart=<offset>`
- `retmax=<batch size>`
- `retmode=xml`

This is the standard approach for large PubMed result sets because it avoids sending every PMID in the URL.

## Why `IncompleteRead(...)` happened

The error:

```text
Failed: IncompleteRead(3715492 bytes read)
```

means the HTTP response from PubMed was cut off before Python finished reading the XML body. This is a transport problem during `efetch`, not a problem with the search term itself.

Likely causes:

- large XML responses for broad terms like `HIV`
- network instability
- server-side connection close during long response bodies

## Hardening added to the script

The script was updated to reduce the chance of partial downloads:

- fetch batch size reduced from `200` to `100`
- explicit retry handling for `IncompleteRead`
- retry handling for HTTP, URL, and XML parse failures
- per-batch delay to stay within PubMed rate limits

## CSV format

The script writes this header to match the official PubMed CSV export as closely as possible:

```text
PMID,Title,Authors,Citation,First Author,Journal/Book,Publication Year,Create Date,PMCID,NIHMS ID,DOI
```

The file is written as `utf-8-sig`, which preserves the BOM seen in the official PubMed CSV download.

## Known limitations

- The output columns match the PubMed website export header, but `Citation` may not be byte-for-byte identical to the website export.
- PubMed API metadata is sufficient to build a close citation string, but the website may apply additional formatting rules that are not exposed directly through E-utilities.
- Very large searches may still need a smaller batch size or resume support if the network is unstable.

## If failures continue

The next practical improvements are:

- add a `--batch-size` CLI argument
- add resume support from the last successful `retstart`
- save intermediate progress periodically during long runs
