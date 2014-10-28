all: compile

tarball:
	tar czvf nllfast14.tar.gz nllfast-4.01dcpl.f *.grid Makefile

compile:
	gfortran -o nllfast_14TeV nllfast-4.01dcpl-14TeV.f

clean:
	rm -rf nllfast_14TeV
