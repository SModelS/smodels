FCC=$(shell ../check_fortran_compiler.sh)

all: pythia_lhe

tarball:
	tar czvf pythia.tar.gz pythia_lhe.f pythia-6.4.27.f Makefile

pythia_lhe:
	$(FCC) -o pythia_lhe pythia_lhe.f pythia-6.4.27.f 

test:
	cp ../../inputFiles/slha/compressedSpec.slha fort.61
	pythia_lhe < ../../etc/pythia_test.card | tee out.log

clean:
	rm -f pythia_lhe
