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

## API examples

This script uses two endpoints, with `esearch.fcgi` used in two different ways.

| Endpoint | Purpose | Example input params | Example output fields |
| --- | --- | --- | --- |
| `esearch.fcgi` | Count how many PubMed records match a query | `db=pubmed`, `term=HIV`, `rettype=count`, `retmode=xml`, `retmax=0` | `Count` |
| `esearch.fcgi` | Create a PubMed history session for a query or date-split query | `db=pubmed`, `term=(HIV) AND ("2024/01/01"[PDAT] : "2024/12/31"[PDAT])`, `retmode=xml`, `usehistory=y`, `retmax=0` | `Count`, `QueryKey`, `WebEnv` |
| `efetch.fcgi` | Fetch actual article XML in batches from a history session | `db=pubmed`, `query_key=1`, `WebEnv=MCID_...`, `retstart=0`, `retmax=100`, `retmode=xml` | `PubmedArticle` records with `PMID`, title, authors, journal, dates, IDs |

### `esearch.fcgi` count example

Example request:

```text
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=HIV&rettype=count&retmode=xml&retmax=0
```

Example response shape:

```xml
<eSearchResult>
  <Count>456703</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <IdList />
</eSearchResult>
```

### `esearch.fcgi` history example

Example request:

```text
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%28HIV%29+AND+%28%222024%2F01%2F01%22%5BPDAT%5D+%3A+%222024%2F12%2F31%22%5BPDAT%5D%29&retmode=xml&usehistory=y&retmax=0
```

Example response shape:

```xml
<eSearchResult>
  <Count>8421</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_...</WebEnv>
  <IdList />
</eSearchResult>
```

### `efetch.fcgi` batch example

Example request:

```text
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&query_key=1&WebEnv=MCID_...&retstart=0&retmax=100&retmode=xml
```

Example response shape:

```xml
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>12345678</PMID>
      <Article>
        <ArticleTitle>Example HIV paper title</ArticleTitle>
      </Article>
    </MedlineCitation>
    <PubmedData>
      <ArticleIdList>
        <ArticleId IdType="doi">10.1000/example</ArticleId>
      </ArticleIdList>
    </PubmedData>
  </PubmedArticle>
</PubmedArticleSet>
```

The script converts each `PubmedArticle` into one CSV row with:

- `PMID`
- `Title`
- `Abstract`
- `Authors`
- `Citation`
- `First Author`
- `Journal/Book`
- `Publication Year`
- `Create Date`
- `PMCID`
- `NIHMS ID`
- `DOI`

## Caching

The script caches all API calls it makes.

- `esearch.fcgi` count requests are cached under `.pubmed_search_cache/esearch.fcgi/`
- `esearch.fcgi` history requests are cached under `.pubmed_search_cache/esearch.fcgi/`
- `efetch.fcgi` batch requests are cached under `.pubmed_search_cache/efetch.fcgi/`

Each request gets its own cache file.

- one API request = one JSON cache file
- one `efetch` batch = one cache file
- one `esearch` query = one cache file

The cache key is based on:

- endpoint name
- request parameters

The cache key ignores:

- `email`
- `api_key`

The previous monolithic cache file `.pubmed_search_cache.json` is not used by the new layout.

For `esearch.fcgi` history calls, the script now uses deterministic staged date ranges so cache keys stay reusable across runs.

- each full year starts as two fixed half-year ranges: `Jan-Jun` and `Jul-Dec`
- if a half-year range exceeds the PubMed limit, that half-year is split into fixed quarters
- if a quarter range still exceeds the limit, that quarter is split into fixed months
- if a month range still exceeds the limit, the script recursively halves that date range until it is under the limit
- the whole search period is covered by this range list, with the final current period clipped to today when needed

There is also a separate PMID index cache in `<stem>_split/<stem>_pmid_index.csv`.

- that file is not an API cache
- it is a local index of PMIDs already written to split CSV files
- it is used to skip processing PMIDs that are already present in existing CSV output

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
PMID,Title,Abstract,Authors,Citation,First Author,Journal/Book,Publication Year,Create Date,PMCID,NIHMS ID,DOI
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
