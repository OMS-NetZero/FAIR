venv: setup.py
	[ -d ./venv ] || python3 -m venv venv
	./venv/bin/pip install --upgrade pip
	./venv/bin/pip install -e .[docs,dev]
	touch venv
