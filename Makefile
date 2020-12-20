.DEFAULT_GOAL := help

VENV_DIR ?= venv
TESTS_DIR=$(PWD)/tests

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

.PHONY: help
help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

checks: $(VENV_DIR)  ## run all the checks
	@echo "\n\n=== black ==="; $(VENV_DIR)/bin/black --check fair tests setup.py --exclude fair/_version.py || echo "--- black failed ---" >&2; \
		echo "\n\n=== flake8 ==="; $(VENV_DIR)/bin/flake8 fair tests setup.py || echo "--- flake8 failed ---" >&2; \
		echo "\n\n=== isort ==="; $(VENV_DIR)/bin/isort --check-only --quiet fair tests setup.py || echo "--- isort failed ---" >&2; \
		echo "\n\n=== tests ==="; $(VENV_DIR)/bin/pytest tests --cov -rfsxEX --cov-report term-missing || echo "--- tests failed ---" >&2; \
		echo "\n\n=== docs ==="; $(VENV_DIR)/bin/sphinx-build -M html docs/source docs/build -qW || echo "--- docs failed ---" >&2; \
		echo

.PHONY: format
format:  ## re-format files
	make isort
	make black

black: $(VENV_DIR)  ## apply black formatter to source and tests
	@status=$$(git status --porcelain setup.py fair tests docs); \
	if test ${FORCE} || test "x$${status}" = x; then \
		$(VENV_DIR)/bin/black --exclude _version.py setup.py fair tests docs; \
	else \
		echo Not trying any formatting. Working directory is dirty ... >&2; \
	fi;

isort: $(VENV_DIR)  ## format the code
	@status=$$(git status --porcelain fair tests setup.py); \
	if test ${FORCE} || test "x$${status}" = x; then \
		$(VENV_DIR)/bin/isort fair tests setup.py; \
	else \
		echo Not trying any formatting. Working directory is dirty ... >&2; \
	fi;

flake: $(VENV_DIR)  ## PEP8 compliance
	$(VENV_DIR)/bin/flake8 fair tests setup.py;


.PHONY: docs
docs: $(VENV_DIR)  ## build the docs
	$(VENV_DIR)/bin/sphinx-build -M html docs/source docs/build

.PHONY: test
test:  $(VENV_DIR) ## run the full testsuite
	$(VENV_DIR)/bin/pytest --cov -r a --cov-report term-missing tests

virtual-environment:  ## update venv, create a new venv if it doesn't exist
	make $(VENV_DIR)

$(VENV_DIR): setup.py
	[ -d $(VENV_DIR) ] || python3 -m venv $(VENV_DIR)

	$(VENV_DIR)/bin/pip install --upgrade pip wheel
	$(VENV_DIR)/bin/pip install -e .[dev]

	touch $(VENV_DIR)
