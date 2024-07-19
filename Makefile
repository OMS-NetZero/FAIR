# Gratitude to IIASA Climate Assessment, which this is based upon

.DEFAULT_GOAL := help

VENV_DIR ?= ./venv

FILES_TO_FORMAT_PYTHON=src tests setup.py

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([0-9a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

.PHONY: help
help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

.PHONY: test
test: $(VENV_DIR)  ## run all the tests and produce a coverage report
	$(VENV_DIR)/bin/pytest -rfsxEX tests -r a --cov=fair --cov-report term-missing:skip-covered -o log_cli=true

.PHONY: checks
checks: $(VENV_DIR)  ## run all the checks
	@echo "=== bandit ==="; $(VENV_DIR)/bin/bandit -c .bandit.yml -r src/fair || echo "--- bandit failed ---" >&2; \
		echo "\n\n=== black ==="; $(VENV_DIR)/bin/black --check src tests setup.py --exclude fair/_version.py || echo "--- black failed ---" >&2; \
		echo "\n\n=== isort ==="; $(VENV_DIR)/bin/isort --check-only --quiet src tests setup.py || echo "--- isort failed ---" >&2; \
		echo "\n\n=== flake8 ==="; $(VENV_DIR)/bin/flake8 src tests setup.py || echo "--- flake8 failed ---" >&2; \
		echo

.PHONY: format
format:  ## re-format files
	make isort
	make black

.PHONY: black
black: $(VENV_DIR)  ## use black to autoformat code
	$(VENV_DIR)/bin/black --target-version py311 $(FILES_TO_FORMAT_PYTHON)

isort: $(VENV_DIR)  ## format the code
	$(VENV_DIR)/bin/isort $(FILES_TO_FORMAT_PYTHON)

virtual-environment: $(VENV_DIR)  ## update venv, create a new venv if it doesn't exist
$(VENV_DIR): setup.py
	[ -d $(VENV_DIR) ] || python3.11 -m venv $(VENV_DIR)

	$(VENV_DIR)/bin/pip install --upgrade pip
	$(VENV_DIR)/bin/pip install wheel
	$(VENV_DIR)/bin/pip install -e .[dev]

	touch $(VENV_DIR)

.PHONY: test_notebooks
test_notebooks: $(VENV_DIR)
	$(VENV_DIR)/bin/pytest  --nbmake ./examples

.PHONY: docs
docs:
	$(VENV_DIR)/bin/nbstripout examples/*.ipynb
	$(VENV_DIR)/bin/jupyter-nbconvert --to rst examples/*.ipynb --output-dir="./docs/examples"
	$(VENV_DIR)/bin/sphinx-build -M html docs docs/_build
