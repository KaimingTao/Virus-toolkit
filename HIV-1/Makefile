TERM := HIV-1
OUT := $(CURDIR)/hiv_pubmed.csv
CACHE_DIR := $(CURDIR)/pubmed_search_cache
SCRIPT := ../pubmed_search.py

.PHONY: pubmed

pubmed:
	uv run $(SCRIPT) --term "$(TERM)" --out "$(OUT)" --cache-dir "$(CACHE_DIR)"
