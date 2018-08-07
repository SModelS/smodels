VER=$(shell cat smodels/version)

all: resolve_deps externaltools

resolve_deps: ## resolve the deps via pip
	@echo "try to resolve the python dependencies via pip"
	smodels/installation.py -R

install: all tidy
	@echo "done. you can now run the software directly from this source directory."
	@echo "Try e.g. ./runSModelS.py --help"

tidy:
	# tidy up the directory, remove files not needed for users
	yes | rm -rf build dist test

version:
	@echo $(VER)

externaltools:
	cd smodels/lib && make

pythia6:
	cd smodels/lib && make pythia6

pythia8:
	cd smodels/lib && make pythia8

clean:
	yes | rm -rf build dist
	cd smodels/lib && make clean

buildrpm:
	$(PYTHON) setup.py bdist_rpm --force-arch x86_64

builddeb: buildrpm
	cd dist && fakeroot alien smodels-$(VER)-1.x86_64.rpm

pypi:
	rm -rf dist
	python setup.py sdist bdist_wheel
	twine upload dist/smodels-*.tar.gz

testpypi: 
	# to install from testpypi: 
	# pip3 install --user --upgrade --index-url https://test.pypi.org/simple/ smodels
	rm -rf dist
	python setup.py sdist bdist_wheel
	twine upload -r pypitest dist/smodels-*.tar.gz

tarballs:
	cd distribution && make tarballs

.PHONY:
