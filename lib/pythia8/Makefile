CXX      := g++
CXXFLAGS := -O2 -I/usr/include/Pythia8
SRC := pythia8.cc

EXES := $(SRC:.cc=.exe)

all: $(EXES)

%.exe: %.cc
	$(CXX) $(CXXFLAGS) -o $@ $< -lpythia8 -lLHAPDF

run:
	pythia8.exe
