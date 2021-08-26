## Makefile, static linking. Pythia8 is fetched and built
## automatically.
## on a mac, for now please use the DYLD_LIBRARY_PATH variable  
## for pythia8.exe to find the shared object.

## PYTHIAVER := `cat pythiaversion`
PYTHIAVER := $(shell ./installer.py -v)

PYTHIA8HOME := pythia$(PYTHIAVER)/
CXXFLAGS := -O3 -I$(PYTHIA8HOME)/include -I$(PYTHIA8HOME)/include/Pythia8/
LDFLAGS  := -L$(PYTHIA8HOME)/lib/ -L$(PYTHIA8HOME)/lib -Wl,-rpath,$(PWD)/$(PYTHIA8HOME)/lib
XMLDOC   := $(PYTHIA8HOME)/share/Pythia8/xmldoc

all: install_pythia8 pythia8.exe

version:
	@echo $(PYTHIAVER)

install_pythia8:
	./installer.py

clean:
	rm -rf $(EXES) test.slha pythia$(PYTHIAVER) .pythia*.tgz* pythia*.tgz*  pythia8.exe events.lhe

pythia8.exe: pythia8.cc
	echo $(XMLDOC) > xml.doc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ pythia8.cc -lpythia8 -ldl

run:
	./pythia8.exe -f ../../../inputFiles/slha/simplyGluino.slha -n 10 -c ../../etc/pythia8.cfg -x $(XMLDOC)

help:
	./pythia8.exe -h

particles:
	## used to create mappings of pythia8 particle names to pdg ids
	## the user shouldnt have to do this
	./parseXml.py

upload_particles: particles
	## only needed if the particle names have changed
	cp pythia8particles.py ../../smodels/tools/
	

test: all run

.PHONY:
