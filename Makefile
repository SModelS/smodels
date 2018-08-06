VER=$(shell cat smodels/version)

all: externaltools

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
	python setup.py sdist bdist_wheel
	twine upload dist/*

testpypi: 
	# to install from testpypi: 
	# pip3 install --user --upgrade --index-url https://test.pypi.org/simple/ smodels
	python setup.py sdist bdist_wheel
	twine upload -r pypitest dist/*

tarballs:
	cd distribution && make tarballs

.PHONY:
