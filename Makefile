venv:
	[ -d ./venv ] || python3 -m venv venv
	./venv/bin/pip install --upgrade pip
