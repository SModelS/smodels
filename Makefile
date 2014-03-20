all: externaltools

externaltools:
	cd tools/external && make

regression: .PHONY
	cd regression && ./run.py
	cd test && ./run_all.py

.PHONY:
