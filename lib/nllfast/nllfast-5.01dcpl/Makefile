all: compile

tarball:
	tar czvf nllfast33.tar.gz nllfast-5.01dcpl.f *.grid Makefile

compile:
	gfortran -o nllfast_33TeV nllfast-5.01dcpl-33TeV.f

clean:
	rm -f nllfast_33TeV
