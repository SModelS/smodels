# --- Check for python3-config --------------------------------------------------

define PYCONFIG_ERROR
python3-config not found. Please install the Python development package:
  • Debian/Ubuntu: sudo apt install python3-dev
  • Fedora/RHEL:  sudo dnf install python3-devel
  • Arch Linux:   sudo pacman -S python
  • macOS:        install Python from python.org or brew install python
endef

PYTHON3_CONFIG := $(shell command -v python3-config 2>/dev/null)

ifeq ($(PYTHON3_CONFIG),)
$(error $(PYCONFIG_ERROR))
endif

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
