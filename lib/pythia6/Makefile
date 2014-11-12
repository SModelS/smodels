all: pythia_lhe

tarball:
	tar czvf pythia.tar.gz pythia_lhe.f pythia-6.4.27.f Makefile

pythia_lhe:
	gfortran -o pythia_lhe pythia_lhe.f pythia-6.4.27.f 

clean:
	rm -f pythia_lhe
