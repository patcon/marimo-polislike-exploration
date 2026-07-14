SHELL := /bin/bash

VALENCY_ANNDATA_VERSION := >=0.4.0
VALENCY_ANNDATA_GIT := https://github.com/patcon/valency-anndata
VALENCY_ANNDATA_LOCAL_PATH := ../valency-anndata

# Picks which notebook*.py to operate on, prompting when more than one
# exists. Override non-interactively with NOTEBOOK=<file>.
define PICK_NOTEBOOK
if [ -n "$(NOTEBOOK)" ]; then nb="$(NOTEBOOK)"; \
else \
	files=(notebook*.py); \
	if [ $${#files[@]} -eq 1 ]; then nb=$${files[0]}; \
	else \
		echo "Select a notebook:" >&2; \
		select nb in "$${files[@]}"; do [ -n "$$nb" ] && break; done; \
	fi; \
fi
endef

install: ## Install dependencies, using the published valency-anndata package
	uv remove valency-anndata --frozen >/dev/null 2>&1 || true
	uv add "valency-anndata$(VALENCY_ANNDATA_VERSION)"

install-local: ## Install valency-anndata as an editable checkout at ../valency-anndata
	uv remove valency-anndata --frozen >/dev/null 2>&1 || true
	uv add --editable $(VALENCY_ANNDATA_LOCAL_PATH)

install-remote-patch: ## Install valency-anndata from a git branch (usage: make install-remote-patch BRANCH=my-branch)
	@test -n "$(BRANCH)" || (echo "Usage: make install-remote-patch BRANCH=<branch>" >&2; exit 1)
	uv remove valency-anndata --frozen >/dev/null 2>&1 || true
	uv add "valency-anndata @ git+$(VALENCY_ANNDATA_GIT)" --branch $(BRANCH)

dev: ## Open a notebook in the marimo editor (interactive, editable; prompts if multiple notebooks exist)
	@$(PICK_NOTEBOOK); \
	uv run marimo edit "$$nb" --watch --no-token

run-app: ## Serve a notebook as a read-only app (prompts if multiple notebooks exist)
	@$(PICK_NOTEBOOK); \
	uv run marimo run "$$nb"

html: ## Export a notebook to a static HTML file (prompts if multiple notebooks exist)
	@$(PICK_NOTEBOOK); \
	uv run marimo export html "$$nb" -o "$${nb%.py}.html"

script: ## Export a notebook to a plain Python script (prompts if multiple notebooks exist)
	@$(PICK_NOTEBOOK); \
	uv run marimo export script "$$nb" -o "$${nb%.py}_script.py"

clean: ## Remove marimo session/autosave state
	rm -rf __marimo__

# These make tasks allow the default help text to work properly.
%:
	@true

.PHONY: install install-local install-remote-patch dev run-app html script clean help

help:
	@echo 'Usage: make <command>'
	@echo
	@echo 'where <command> is one of the following:'
	@echo
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
