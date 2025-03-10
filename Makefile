VER=$(shell cat smodels/version)

HAS_FC := $(shell smodels/lib/check_fortran_compiler.sh 2> /dev/null)
HAS_CXX := $(shell command -v $(CXX) 2> /dev/null)

all: resolve_deps # make all just resolves dependencies

check_compilers: .PHONY
ifndef HAS_FC
	$(error "Fortran compiler not found. Cannot compile external tools. You may try to give the explicit path to the compiler via the FC variable (make FC=... smodels). Alternatively, you can still build smodels proper, via 'make smodels_noexternaltools'" )
endif
ifndef HAS_CXX
	$(error "C++ compiler not found. Cannot compile external tools. You may try to give the explicit path to the compiler via the CXX variable (make CXX=... smodels). Alternatively, you can still build smodels proper, via 'make smodels_noexternaltools'" )
endif

resolve_deps: ## resolve the deps via pip
	@echo "trying to resolve the python dependencies via pip:"
	@smodels/installation.py -R

smodels_externaltools: resolve_deps externaltools
	@echo
	@echo "done. you can now run the software directly from this source directory.\n"
	@echo "Try e.g. \n\n ./runSModelS.py --help\n"
	@echo "The latest SModelS documentation can be found at: http://smodels.readthedocs.io/en/latest/"
	@echo "For this version documentation go to: https://smodels.readthedocs.io/en/v$(VER)"

smodels: resolve_deps

smodels_noexternaltools: resolve_deps
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

resummino:
	cd smodels/lib && make resummino

nllfast:
	cd smodels/lib && make nllfast

cpp: .PHONY
	cd cpp && make

clean:
	yes | rm -rf build dist
	cd smodels/lib && make clean

buildrpm:
	$(PYTHON) setup.py bdist_rpm --force-arch x86_64

builddeb: buildrpm
	cd dist && fakeroot alien smodels-$(VER)-1.x86_64.rpm

pypi: clean
	## pypi user is walten, repository is https://upload.pypi.org/legacy/
	python -m pip install build
	rm -rf dist
#	python3 setup.py sdist bdist_wheel
	python -m build .
	twine check --strict dist/*
	twine upload dist/smodels-*.tar.gz

testpypi: clean
	## testpypi user is smodels, repository is https://test.pypi.org/legacy/
	# to install from testpypi:
	# pip3 install --user --upgrade --index-url https://test.pypi.org/simple/ smodels
	python -m pip install build
	rm -rf dist
#python3 setup.py sdist bdist_wheel
	python -m build .
	twine check --strict dist/*
	twine upload -r pypitest dist/smodels-*.tar.gz

testpypiold: clean
	## testpypi user is smodels, repository is https://test.pypi.org/legacy/
	# to install from testpypi:
	# pip3 install --user --upgrade --index-url https://test.pypi.org/simple/ smodels
	rm -rf dist
	python3 setup.py sdist bdist_wheel
	twine upload -r pypitest dist/smodels-*.tar.gz


tarballs:
	cd distribution && make tarballs

.PHONY:
