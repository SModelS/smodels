CXX      := g++
CXXFLAGS := -O2 -I/usr/include/Pythia8
SRC := pythia8.cc

EXES := $(SRC:.cc=.exe)

all: $(EXES)

clean:
	rm -rf $(EXES) test.slha

%.exe: %.cc
	$(CXX) $(CXXFLAGS) -o $@ $< -lpythia8 -lLHAPDF

run:
	pythia8.exe -h

copy_slha:
	cp ../../inputFiles/slha/simplyGluino.slha test.slha

test: all copy_slha run

.PHONY:
