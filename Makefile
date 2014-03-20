all: pythia nllfast
nllfast: nllfast-1.2 nllfast-2.1 nllfast-3.0 nllfast-4.01dcpl nllfast-5.01dcpl

pythia:
	cd tools/external/pythia6 && make
nllfast-1.2:
	cd tools/external/nllfast/nllfast-1.2 && make
nllfast-2.1:
	cd tools/external/nllfast/nllfast-2.1 && make
nllfast-3.0:
	cd tools/external/nllfast/nllfast-3.0 && make
nllfast-4.01dcpl:
	cd tools/external/nllfast/nllfast-4.01dcpl && make
nllfast-5.01dcpl:
	cd tools/external/nllfast/nllfast-5.01dcpl && make
regression: .PHONY
	cd regression && ./run.py
	cd test && ./run_all.py

.PHONY:
