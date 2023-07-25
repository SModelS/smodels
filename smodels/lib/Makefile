all: pythia nllfast resummino
clean: clean_pythia clean_nllfast clean_resummino

nllfast: nllfast-1.2 nllfast-2.1 nllfast-3.1
clean_nllfast: clean_nllfast-1.2 clean_nllfast-2.1 clean_nllfast-3.1

pythia: pythia8 pythia6
clean_pythia: clean_pythia8 clean_pythia6 ## we clean both

resummino: resummino-3.1.2
clean_resummino: clean_resummino-3.1.2

pythia6: .PHONY
	cd pythia6 && make
clean_pythia6:
	cd pythia6 && make clean
pythia8: .PHONY
	cd pythia8 && make
clean_pythia8:
	cd pythia8 && make clean
nllfast-1.2:
	cd nllfast/nllfast-1.2 && make
clean_nllfast-1.2:
	cd nllfast/nllfast-1.2 && make clean
nllfast-2.1:
	cd nllfast/nllfast-2.1 && make
clean_nllfast-2.1:
	cd nllfast/nllfast-2.1 && make clean
nllfast-3.1:
	cd nllfast/nllfast-3.1 && make
clean_nllfast-3.1:
	cd nllfast/nllfast-3.1 && make clean
resummino-3.1.2:
	cd resummino && make
clean_resummino-3.1.2:
	cd resummino && make clean
.PHONY:
