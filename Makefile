all: pythia nllfast

pythia:
	cd pythia_src && make

nllfast:
	cd nllfast && make
