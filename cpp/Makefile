python_version_full := $(wordlist 2,4,$(subst ., ,$(shell python3 --version 2>&1)))
python_version_major := $(word 1,${python_version_full})
python_version_minor := $(word 2,${python_version_full})
# libpython := -lpython$(python_version_major).$(python_version_minor)
libpython := $(shell ./find_libpython.py)
LDFLAGS := $(shell python3-config --ldflags)
#LIBDIR := $(shell python3-config --prefix)
CFLAGS := $(shell python3-config --cflags)
#EXECPATH := $(shell python3-config --exec-prefix)/lib/python$(python_version_major).$(python_version_minor)

# GCC=g++

all: library executable

clean:
	rm -rf run libSModelS.so

library:
	$(CXX) SModelS.cc -I. $(CFLAGS) $(LDFLAGS) $(libpython) -shared -fPIC -o libSModelS.so

executable:
	@echo "run.cpp is meant as an example, tweak the code and integrate it in your framework!"
	$(CXX) run.cpp -I. -Wl,-rpath,`pwd` -L. -o run -ldl -lpthread -lutil -l SModelS
	@echo "done. now adapt parameters.ini and execute ./run"

run: .PHONY
	./run

.PHONY:
