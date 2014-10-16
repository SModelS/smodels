all: pythia nllfast
clean: clean_pythia clean_nllfast

nllfast: nllfast-1.2 nllfast-2.1 nllfast-3.0 nllfast-4.01dcpl nllfast-5.01dcpl
clean_nllfast: clean_nllfast-1.2 clean_nllfast-2.1 clean_nllfast-3.0 clean_nllfast-4.01dcpl clean_nllfast-5.01dcpl
pythia: pythia6
clean_pythia: clean_pythia6

pythia6:
	cd pythia6 && make
clean_pythia6:
	cd pythia6 && make clean

nllfast-1.2:
	cd nllfast/nllfast-1.2 && make
clean_nllfast-1.2:
	cd nllfast/nllfast-1.2 && make clean
nllfast-2.1:
	cd nllfast/nllfast-2.1 && make
clean_nllfast-2.1:
	cd nllfast/nllfast-2.1 && make clean
nllfast-3.0:
	cd nllfast/nllfast-3.0 && make
clean_nllfast-3.0:
	cd nllfast/nllfast-3.0 && make clean
nllfast-4.01dcpl:
	cd nllfast/nllfast-4.01dcpl && make
clean_nllfast-4.01dcpl:
	cd nllfast/nllfast-4.01dcpl && make clean
nllfast-5.01dcpl:
	cd nllfast/nllfast-5.01dcpl && make
clean_nllfast-5.01dcpl:
	cd nllfast/nllfast-5.01dcpl && make clean

.PHONY:
