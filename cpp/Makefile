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
	$(GCC) SModelS.cc -I. $(libpython) -shared -fPIC -o libSModelS.so

exe:
	## run.cpp is meant as an example, tweak the code and integrate it in your framework!
	$(GCC) run.cpp -I. -Wl,-rpath,`pwd` -L. -o run -ldl -lpthread -lutil -l SModelS
	## done. now adapt parameters.ini and execute ./run

.PHONY:
