all: externaltools

externaltools:
	cd lib && make

builddeb:
	# build the source package in the parent directory
	# then rename it to project_version.orig.tar.gz
	$(PYTHON) setup.py sdist $(COMPILE) --dist-dir=../ ## --prune
	rename -f 's/$(PROJECT)-(.*)\.tar\.gz/$(PROJECT)_$$1\.orig\.tar\.gz/' ../*
	# build the package
	dpkg-buildpackage -i -I -rfakeroot

buildrpm:
	$(PYTHON) setup.py bdist_rpm --force-arch x86_64  ## --post-install=rpm/postinstall --pre-uninstall=rpm/preuninstall
##	$(PYTHON) setup.py bdist --formats rpm ## --post-install=rpm/postinstall --pre-uninstall=rpm/preuninstall


regression: .PHONY
	cd regression && ./run.py
	cd test && ./run_all.py

.PHONY:
