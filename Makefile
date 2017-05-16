VER=$(shell cat smodels/version)

version:
	@echo $(VER)

all: externaltools

externaltools:
	cd lib && make

pythia6:
	cd lib && make pythia6

pythia8:
	cd lib && make pythia8

clean:
	cd lib && make clean

buildrpm:
	$(PYTHON) setup.py bdist_rpm --force-arch x86_64

builddeb: buildrpm
	cd dist && fakeroot alien smodels-$(VER)-1.x86_64.rpm

tarballs:
	cd distribution && make tarballs

.PHONY:
