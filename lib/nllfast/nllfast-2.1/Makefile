all: compile

tarball:
	tar czvf nllfast8.tar.gz nllfast-2.1.f *.grid Makefile

compile:
	gfortran -o nllfast_8TeV nllfast-2.1.f

clean:
	rm -f nllfast_8TeV
