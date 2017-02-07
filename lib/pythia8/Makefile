CXX      := g++
CXXFLAGS := -O2 -I/usr/include/Pythia8
LDFLAGS  := -L/usr/lib/
SRC := pythia8.cc

EXES := $(SRC:.cc=.exe)

all: $(EXES)

clean:
	rm -rf $(EXES) test.slha

%.exe: %.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< -lpythia8 -lLHAPDF

run:
	pythia8.exe -f test.slha -n 10 -c ./pythia8.cfg

help:
	pythia8.exe -h

copy_slha:
	cp ../../inputFiles/slha/simplyGluino.slha test.slha

particles:
	./parseXml.py

upload_particles: particles
	cp pythia8particles.py ../../smodels/tools/
	

test: all copy_slha run

.PHONY:
