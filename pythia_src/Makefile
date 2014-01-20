all: compile

tarball:
	tar czvf pythia.tar.gz pythia_lhe.f pythia-6.4.27.f Makefile

compile:
	gfortran -o ../pythia_lhe pythia_lhe.f pythia-6.4.27.f 
