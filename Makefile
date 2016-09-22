all: externaltools

externaltools:
	cd lib && make

buildrpm:
	$(PYTHON) setup.py bdist_rpm --force-arch x86_64  ## --post-install=rpm/postinstall --pre-uninstall=rpm/preuninstall
##	$(PYTHON) setup.py bdist --formats rpm ## --post-install=rpm/postinstall --pre-uninstall=rpm/preuninstall

tarballs:
	cd distribution && make tarballs

.PHONY:
