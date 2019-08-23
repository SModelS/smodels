#!/usr/bin/env python3
 
"""
.. module:: testCombined
   :synopsis: Tests the combination code
 
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
 
"""
 
import sys,os
import importlib
sys.path.insert(0,"../")
import unittest
from unitTestHelpers import equalObjs, runMain, importModule
from smodels.tools.smodelsLogging import logger, setLogLevel
 
class CombinedTest(unittest.TestCase):

    def defineTest(self):
        """ define the current output as the default output.
        Use with care! """
        filename = "./testFiles/slha/gluino_squarks.slha"
        outputfile = runMain(filename,inifile="testParameters_agg.ini", suppressStdout=True )
        smodelsOutput = importModule ( outputfile )
        f=open("gluino_squarks_default_agg.py","w")
        f.write ( "smodelsOutputDefault = %s\n" % smodelsOutput )
        f.close()

    def testCombinedResult(self):
        filename = "./testFiles/slha/gluino_squarks.slha"
        outputfile = runMain(filename,inifile="testParameters_agg.ini", suppressStdout=True )
        smodelsOutput = importModule ( outputfile )
        from gluino_squarks_default_agg import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'Element', 'database version', 'Total missed xsec', 
                        'Missed xsec long-lived', 'Missed xsec displaced', 'Missed xsec MET', 'Total outside grid xsec',
                        'Total xsec for missing topologies (fb)','Total xsec for missing topologies with displaced decays (fb)',
                        'Total xsec for missing topologies with prompt decays (fb)', 
                        'Total xsec for topologies outside the grid (fb)']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields, fname = outputfile )
        if equals != True:
            logger.error ( "%s differs from %s!" % ( "gluino_squarks_default_agg.py", outputfile) ) 
        self.assertTrue(equals)
        for i in [ outputfile, outputfile.replace(".py",".pyc") ]:
            if os.path.exists ( i ):
                os.remove ( i )
 
if __name__ == "__main__":
    setLogLevel ( "debug" )
    unittest.main()
