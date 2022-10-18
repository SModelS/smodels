all: test

test:
	./runCompleteTestSuite.py -v

test_norecipes:
	./runCompleteTestSuite.py -S -v

reduced:
	./runCompleteTestSuite.py -S -r -v

obs:
	py.test -v

run_coverage:
	coverage run runCompleteTestSuite.py

coverage:
	coverage report --omit="/usr/lib/*,testFiles/scanExample/*,*/.local/*,unitTestOutput/*,test*py,test*py" -i
	rm -rf htmlcov/
	coverage html --omit="/usr/lib/*,testFiles/scanExample/*,*/.local/*,unitTestOutput/*,test*py,test*py" -i
	cp htmlcov/* ../../smodels.github.io/coverage/

notebook:
	ipython notebook
