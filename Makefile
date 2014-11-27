all: externaltools

externaltools:
	cd lib && make

buildrpm:
	$(PYTHON) setup.py bdist_rpm --force-arch x86_64  ## --post-install=rpm/postinstall --pre-uninstall=rpm/preuninstall
##	$(PYTHON) setup.py bdist --formats rpm ## --post-install=rpm/postinstall --pre-uninstall=rpm/preuninstall

regression: .PHONY
	cd regression && ./run.py
	cd test && ./run_all.py

.PHONY:
