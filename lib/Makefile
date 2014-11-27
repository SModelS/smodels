all: pythia nllfast
clean: clean_pythia clean_nllfast

nllfast: nllfast-1.2 nllfast-2.1
clean_nllfast: clean_nllfast-1.2 clean_nllfast-2.1
pythia: pythia6
clean_pythia: clean_pythia6

pythia6: .PHONY
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

.PHONY:
