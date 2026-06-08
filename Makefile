.PHONY: run check clean

CONFIG ?= configs/default.ini

run:
	./scripts/run_case.sh $(CONFIG)

check:
	python3 -m py_compile scripts/postprocess.py
	bash -n scripts/run_case.sh
	python3 scripts/postprocess.py tests/fixtures/diagnostics.csv

clean:
	rm -rf build tmp
