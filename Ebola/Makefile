TERM := Ebola
BASE_OUT := $(CURDIR)/ebola_pubmed.csv
RESULTS_DIR := $(CURDIR)/pubmed_search
CACHE_DIR := $(CURDIR)/pubmed_search_cache
SCRIPT := ../literature_scripts/pubmed_search.py

.PHONY: pubmed

pubmed:
	@echo "Writing Ebola PubMed CSV chunks to $(RESULTS_DIR)"
	uv run $(SCRIPT) --term "$(TERM)" --out "$(BASE_OUT)" --cache-dir "$(CACHE_DIR)"
