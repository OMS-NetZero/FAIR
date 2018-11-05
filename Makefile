venv: setup.py
	[ -d ./venv ] || python3 -m venv venv
	./venv/bin/pip install --upgrade pip
	./venv/bin/pip install -e .[docs,dev,test]
	touch venv

.PHONY: test
test: venv
	./venv/bin/pytest -rfsxEX tests

.PHONY: test_notebooks
test_notebooks: venv
	./venv/bin/pytest -rfsxEX --nbval ./notebooks --sanitize ./notebooks/tests_sanitize.cfg

.PHONY: docs
docs:
	./venv/bin/sphinx-build -M html docs docs/_build
