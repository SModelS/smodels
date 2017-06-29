VER=$(shell cat smodels/version)

all: externaltools

version:
	@echo $(VER)

externaltools:
	cd smodels/lib && make

pythia6:
	cd smodels/lib && make pythia6

pythia8:
	cd smodels/lib && make pythia8

clean:
	cd smodels/lib && make clean

buildrpm:
	$(PYTHON) setup.py bdist_rpm --force-arch x86_64

builddeb: buildrpm
	cd dist && fakeroot alien smodels-$(VER)-1.x86_64.rpm

tarballs:
	cd distribution && make tarballs

.PHONY:
