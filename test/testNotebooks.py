#!/usr/bin/env python3

"""
.. module:: testNotebooks
   :synopsis: Tests that the recipes runs

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import unittest
import os
testdir = os.getcwd()

class NotebookTest(unittest.TestCase):
    nbdir = "../docs/manual/source/recipes"

    def restViaPyTest ( self ):
        """ test the notebooks by calling pytest --nbmake """
        import sys, time
        os.chdir( self.nbdir )
        cmd = "pytest --nbmake"
        cmd = "pytest --nbmake ascii_graph_from_lhe.ipynb"
        import subprocess
        ps = subprocess.Popen ( cmd, shell=True, stdout = subprocess.PIPE, \
                stderr = subprocess.STDOUT )
        os.chdir ( testdir )
        #time.sleep(.1)
        #for c in iter(lambda: ps.stdout.read(1), b''):
        #    sys.stdout.buffer.write(c)
        #    sys.stdout.buffer.flush()

    def restIpynb ( self ):
        """ test the notebooks via the testipynb module """
        try:
            import testipynb
        except ImportError as e:
            print ( "testipynb is not installed. You may want to do:" )
            print ( "pip install testipynb"  )
        ignores = []
        ignores += [ "compareUL" ]
        ignores += [ "compute_likelihood" ]
        ignores += [ "runWithParameterFile" ]
        ignores += [ "runAsLibrary" ]
        ignores += [ "browserExample3" ]
        ignores += [ "missingTopologies" ]
        ignores += [ "marginalize" ]
        ignores += [ "print_theoryPrediction" ]
        ignores += [ "load_database" ]
        ignores += [ "lheLLPExample" ]
        Test = testipynb.TestNotebooks ( directory = self.nbdir, ignore=ignores, 
                                         timeout=3 )
        self.assertTrue ( Test.run_tests() )
        os.chdir ( testdir )

    def listOfNotebooks ( self ):
        import glob
        notebooks = glob.glob ( f"{self.nbdir}/*.ipynb" )
        notebooks = [ x[:-6].replace(self.nbdir+"/", "") for x in notebooks ]
        return notebooks

    def testRun( self ):
        """ test via the ipynb module, which overrides importlib
            to enable import of ipynb files """
        import importlib
        try:
            import ipynb
        except ImportError as e:
            print ( "ipynb is not installed. You may wish to do:" )
            print ( "pip install ipynb"  )
        to_unlink = [] # keep track of files that need to be unlinked
        try:
            import redirector
            import os
            notebooks = self.listOfNotebooks()
            for x in [ "smodels_paths.py", "inputFiles" ]:
                if os.path.exists ( x ):
                    os.unlink ( x )
                os.symlink ( f"{self.nbdir}/{x}", x )
                to_unlink.append ( x )
            pre = "ipynb.fs.full"
            for nb in notebooks: 
                # for now check only the ascii graph notebook!
                if not "ascii_" in nb:
                    continue
                nbfile = nb + ".ipynb"
                # print ( "Testing", nb )
                module = f"{pre}.{nb}"
                if os.path.exists ( nbfile ):
                    os.unlink ( nbfile )
                os.symlink ( f"{self.nbdir}/{nb}.ipynb", nbfile )
                to_unlink.append ( nbfile )
                to = os.devnull
                with redirector.stdout_redirected ( to = to ):
                    a=importlib.import_module ( module )
                os.unlink ( nbfile )
                to_unlink.remove ( nbfile )
            for x in to_unlink:
                if os.path.exists ( x ):
                    os.unlink ( x )
        except Exception as e:
            for x in to_unlink:
                if os.path.exists ( x ):
                    os.unlink ( x )
            raise(e)

if __name__ == "__main__":
    unittest.main()
