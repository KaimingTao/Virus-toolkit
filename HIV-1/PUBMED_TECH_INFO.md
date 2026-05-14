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

## PubMed ESearch limit

The current script uses `esearch.fcgi` on `db=pubmed` and then pages through the saved result set with `WebEnv` and `QueryKey`.

This has an important PubMed limitation:

- for PubMed, `ESearch` can only retrieve the first `10,000` records matching a query

This is documented by NCBI in the E-utilities guide:

- https://www.ncbi.nlm.nih.gov/books/NBK25499/

Practical consequence:

- a very broad query like `HIV` may show hundreds of thousands of results on the PubMed website
- the current script can only access roughly the first `10,000` of those results from one `ESearch` history session

That is why a website search such as:

```text
https://pubmed.ncbi.nlm.nih.gov/?term=hiv&filter=years.1984-2026
```

can show about `456703` results, while the current CSV export run produced only about `9968` records plus the header row.

## Workaround for full exports

To retrieve all PubMed results for a broad term, the query must be split into smaller subqueries so that each subquery returns fewer than `10,000` records.

Typical approaches:

- split by publication year
- split by publication month for years that still exceed `10,000`
- merge all subquery outputs into one CSV

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
