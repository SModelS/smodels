VER=$(shell cat smodels/version)

HAS_FC := $(shell smodels/lib/check_fortran_compiler.sh 2> /dev/null)
HAS_CXX := $(shell command -v $(CXX) 2> /dev/null)

all: resolve_deps externaltools

check_compilers: .PHONY
ifndef HAS_FC
	$(error "Fortran compiler not found. Cannot compile external tools. You may try to give the explicit path to the compiler via the FC variable (make FC=... smodels). Alternatively, you can still build smodels proper, via 'make smodels_noexternaltools'" )
endif
ifndef HAS_CXX
	$(error "C++ compiler not found. Cannot compile external tools. You may try to give the explicit path to the compiler via the CXX variable (make CXX=... smodels). Alternatively, you can still build smodels proper, via 'make smodels_noexternaltools'" )
endif

resolve_deps: ## resolve the deps via pip
	@echo "try to resolve the python dependencies via pip"
	smodels/installation.py -R

smodels: all tidy
	@echo
	@echo "done. you can now run the software directly from this source directory.\n"
	@echo "Try e.g. \n\n ./runSModelS.py --help\n"
	@echo "The latest SModelS documentation can be found at: http://smodels.readthedocs.io/en/latest/"
	@echo "For this version documentation go to: https://smodels.readthedocs.io/en/v$(VER)"

smodels_noexternaltools: resolve_deps tidy
	@echo
	@echo "done. you can now run the software directly from this source directory.\n"
	@echo "Try e.g. \n\n ./runSModelS.py --help\n"
	@echo "The latest SModelS documentation can be found at: http://smodels.readthedocs.io/en/latest/"
	@echo "For this version documentation go to: https://smodels.readthedocs.io/en/v$(VER)"

tidy:
	# tidy up the directory, remove files not needed for users
	yes | rm -rf build dist test docs apt.txt environment.yml .git .gitattributes .gitignore smodels.egg*

version:
	@echo $(VER)

externaltools: check_compilers
	cd smodels/lib && make -j 4

pythia6:
	cd smodels/lib && make pythia6

pythia8:
	cd smodels/lib && make pythia8

cpp: .PHONY
	cd cpp && make

clean:
	yes | rm -rf build dist
	cd smodels/lib && make clean

buildrpm:
	$(PYTHON) setup.py bdist_rpm --force-arch x86_64

builddeb: buildrpm
	cd dist && fakeroot alien smodels-$(VER)-1.x86_64.rpm

pypi:
	## pypi user is walten, repository is https://upload.pypi.org/legacy/
	rm -rf dist
	python3 setup.py sdist bdist_wheel
	twine upload dist/smodels-*.tar.gz

testpypi: 
	## testpypi user is smodels, repository is https://test.pypi.org/legacy/
	# to install from testpypi: 
	# pip3 install --user --upgrade --index-url https://test.pypi.org/simple/ smodels
	rm -rf dist
	python3 setup.py sdist bdist_wheel
	twine upload -r pypitest dist/smodels-*.tar.gz

tarballs:
	cd distribution && make tarballs

.PHONY:
