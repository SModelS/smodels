FCC=$(shell ../check_fortran_compiler.sh)

all: pythia_lhe

tarball:
	tar czvf pythia.tar.gz pythia_lhe.f pythia-6.4.27.f Makefile

pythia_lhe:
	$(FCC) -o pythia_lhe pythia_lhe.f pythia-6.4.27.f 

test:
	cp ../../../inputFiles/slha/complicated.slha fort.61
	cat ../../etc/pythia.card | sed -e 's/NEVENTS/100/' | sed -e 's/SQRTSD0/13000/' > pythia_test.card
	pythia_lhe < pythia_test.card | tee out.log
	rm -f pythia_test.card

clean:
	rm -f pythia_lhe
