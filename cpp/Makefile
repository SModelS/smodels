python_version_full := $(wordlist 2,4,$(subst ., ,$(shell python --version 2>&1)))
python_version_major := $(word 1,${python_version_full})
python_version_minor := $(word 2,${python_version_full})
libpython := -lpython$(python_version_major).$(python_version_minor)

GCC=g++

all: library exe

library:
	$(GCC) SModelS.cc -I. $(libpython) -shared -fPIC -o libSModelS.so

exe:
	$(GCC) run.cpp -I. -L. -o run -l SModelS

run: .PHONY
	LD_LIBRARY_PATH=. ./run

.PHONY:
