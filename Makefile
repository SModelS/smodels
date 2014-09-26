all: externaltools

externaltools:
	cd lib && make

regression: .PHONY
	cd regression && ./run.py
	cd test && ./run_all.py

.PHONY:
