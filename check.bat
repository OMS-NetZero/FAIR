REM runs formatting and tests on Windows.
REM to run this script, you should first create a new environment (e.g. with conda)
REM and run
REM     pip install -e .[dev]

flake8 fair tests setup.py
black --exclude _version.py setup.py fair tests docs
isort fair tests setup.py
pytest --cov -r a --cov-report term-missing tests