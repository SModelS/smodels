all: externaltools

externaltools:
	cd lib && make

clean:
	cd lib && make clean

buildrpm:
	$(PYTHON) setup.py bdist_rpm --force-arch x86_64

builddeb: buildrpm
	cd dist && fakeroot alien smodels-1.0.93-1.x86_64.rpm

tarballs:
	cd distribution && make tarballs

.PHONY:
