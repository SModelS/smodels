#!/usr/bin/env python3

"""
.. module:: testNotebooks
   :synopsis: Tests that the recipes runs

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import unittest
import os
import sys
import subprocess
sys.path.append('../')
from smodels.tools.smodelsLogging import logger
from smodels.tools.smodelsLogging import setLogLevel
setLogLevel( "debug" )



class NotebookTest(unittest.TestCase):
    nbdir = "../docs/manual/source/recipes"

    def listOfNotebooks(self):
        import glob
        notebooks = glob.glob(f"{self.nbdir}/*.ipynb")
        notebooks = list(notebooks)
        return notebooks

    def cleanUp(self, to_unlink):
        """ clean up files in to_unlink, heeding wildcards """
        import glob
        import shutil
        for x in to_unlink:
            files = glob.glob(x)
            for f in files:
                if os.path.isdir(f) and not os.path.islink(f):
                    shutil.rmtree(f)
                else:
                    os.unlink(f)
        to_unlink = []
        return to_unlink

    def testRunNotebooks(self):

        try:
            import pytest
        except:
            logger.error("pip install pytest and nbmake for testing notebooks")
            return

        notebooks = self.listOfNotebooks()
        for notebookFile in notebooks:
            try:
                r = subprocess.check_output(["pytest --nbmake %s" %notebookFile], shell=True)
            except subprocess.CalledProcessError as e:
                logger.debug(e.output.decode())
                self.assertEqual(e.returncode,0)


    # def testRun(self):
    #     """ test via the ipynb module, which overrides importlib
    #         to enable import of ipynb files """
    #     import importlib
    #     try:
    #         import ipynb
    #     except ImportError as e:
    #         print("ipynb is not installed. You may wish to do:")
    #         print("pip install ipynb")
    #     to_unlink = []  # keep track of files that need to be unlinked
    #     try:
    #         import redirector
    #         import os
    #         cwd = os.getcwd()
    #         if cwd.endswith ( "/test" ):
    #             os.chdir ( "../docs/manual/source/recipes/" )
    #         notebooks = self.listOfNotebooks()
    #         #for x in [ "smodels_paths.py" ]:
    #         # for x in ["smodels_paths.py", "inputFiles"]:
    #         #    if os.path.exists(x):
    #         #        os.unlink(x)
    #         #    os.symlink(f"{self.nbdir}/{x}", x)
    #         #    to_unlink.append(x)
    #         pre = "ipynb.fs.full"
    #         for nb in notebooks:
    #             # for now check only the ascii graph notebook!
    #             good = ["lookup_efficiency", "interactivePlotsExample",
    #                     "browserExample2", "print_decomposition",
    #                     "ascii_graph_from_lhe",
    #                     "nll_xsecs_from_slha", "lookup_upper_limit"]
    #             # good = [ "lo_xsecs_from_slha" ]
    #             # bad because much memory FIXME should eventually
    #             # also be run
    #             bad = ["load_database", "missingTopologies", "lheLLPExample",
    #                    "print_theoryPrediction", "browserExample3",
    #                    "runAsLibrary", "runWithParameterFile",
    #                    "compareUL", "compute_likelihood", "lo_xsecs_from_slha"]
    #             # if not nb in good:
    #             #    continue
    #             if nb in bad:
    #                 continue
    #             #if os.getcwd().endswith ( "/test" ):
    #             #    f = open ( "smodels_paths.py", "wt" )
    #             #    f.close()
    #                 #os.chdir ( "../docs/manual/source/recipes/" )
    #             nbfile = nb + ".ipynb"
    #             module = f"{pre}.{nb}"
    #             if os.path.exists(nbfile):
    #                 os.unlink(nbfile)
    #             os.symlink(f"{self.nbdir}/{nb}.ipynb", nbfile)
    #             to_unlink.append(nbfile)
    #             to = os.devnull
    #             with redirector.stdout_redirected(to=to):
    #                 a = importlib.import_module(module)
    #             os.unlink(nbfile)
    #             to_unlink.remove(nbfile)
    #             if nb == "interactivePlotsExample":
    #                 self.cleanUp(["*.html", "data_frame.txt", "iplots/"])
    #         os.chdir ( cwd )
    #         to_unlink = self.cleanUp(to_unlink)
    #     except Exception as e:
    #         self.cleanUp(to_unlink)
    #         raise(e)


if __name__ == "__main__":
    unittest.main()
