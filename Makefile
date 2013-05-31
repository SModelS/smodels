all: pythia nllfast
nllfast: nllfast-1.2 nllfast-2.1

pythia:
	cd pythia_src && make
nllfast-1.2:
	cd nllfast/nllfast-1.2 && make
nllfast-2.1:
	cd nllfast/nllfast-2.1 && make
regression: .PHONY
	cd regression && ./run.py
	cd test && ./run_all.py

.PHONY:
