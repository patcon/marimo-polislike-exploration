install: ## Install dependencies
	uv sync

dev: ## Open the notebook in the marimo editor (interactive, editable)
	uv run marimo edit notebook.py --watch --no-token

run-app: ## Serve the notebook as a read-only app
	uv run marimo run notebook.py

html: ## Export the notebook to a static HTML file (notebook.html)
	uv run marimo export html notebook.py -o notebook.html

script: ## Export the notebook to a plain Python script (notebook_script.py)
	uv run marimo export script notebook.py -o notebook_script.py

clean: ## Remove marimo session/autosave state
	rm -rf __marimo__

# These make tasks allow the default help text to work properly.
%:
	@true

.PHONY: install dev run-app html script clean help

help:
	@echo 'Usage: make <command>'
	@echo
	@echo 'where <command> is one of the following:'
	@echo
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
