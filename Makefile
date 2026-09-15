.PHONY: test docs

test:
	python -m pytest -q

docs:
	python -m mkdocs build
