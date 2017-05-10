python_version_full := $(wordlist 2,4,$(subst ., ,$(shell python --version 2>&1)))
python_version_major := $(word 1,${python_version_full})
python_version_minor := $(word 2,${python_version_full})
libpython := -lpython$(python_version_major).$(python_version_minor)
#LDFLAGS := $(shell python-config --ldflags)
#LIBDIR := $(shell python-config --prefix)
#CFLAGS := $(shell python-config --cflags)
#EXECPATH := $(shell python-config --exec-prefix)/lib/python$(python_version_major).$(python_version_minor)

GCC=g++

all: library exe

library:
	$(GCC) SModelS.cc -I. $(libpython) -shared -fPIC -o libSModelS.so # $(LDFLAGS) -L$(LIBDIR)/lib

exe:
	$(GCC) run.cpp -I. -Wl,-rpath,`pwd` -L. -o run -l SModelS # $(LDFLAGS) -L$(LIBDIR)/lib

run: .PHONY 
  # consider setting also PYTHONPATH or even PYTHONHOME
	LD_LIBRARY_PATH=. ./run

.PHONY:
