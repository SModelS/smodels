#!/usr/bin/env python

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys,shutil,os
sys.path.insert(0,"../")
import unittest
import glob
from os.path import join, basename
from smodels.installation import installDirectory as iDir
from smodels.tools import crashReport
from smodels.tools.timeOut import NoTime
from runSModelS import main
import logging as logger
import redirector
import unum
import time

def equalObjs(obj1,obj2,allowedDiff,ignore=[]):
    """
    Compare two objects.
    The numerical values are compared up to the precision defined by allowedDiff.

    :param obj1: First python object to be compared
    :param obj2: Second python object to be compared
    :param allowedDiff: Allowed % difference between two numerical values
    :param ignore: List of keys to be ignored
    :return: True/False
    """

    if type(obj1) != type(obj2):
        logger.info("Data types differ (%s,%s)" %(type(obj1),type(obj2)))
        return False

    if isinstance(obj1,unum.Unum):
        if obj1 == obj2:
            return True
        diff = 2.*abs(obj1-obj2)/abs(obj1+obj2)
        return diff.asNumber() < allowedDiff
    elif isinstance(obj1,float):
        if obj1 == obj2:
            return True
        diff = 2.*abs(obj1-obj2)/abs(obj1+obj2)
        return diff < allowedDiff
    elif isinstance(obj1,str):
        return obj1 == obj2
    elif isinstance(obj1,dict):
        for key in obj1:
            if key in ignore: continue
            if not key in obj2:
                logger.warning("Key %s missing" %key)
                return False
            if not equalObjs(obj1[key],obj2[key],allowedDiff, ignore=ignore ):
                logger.warning( "Dictionaries differ in key ``%s''" % key )
                s1,s2 = str(obj1[key]),str(obj2[key]) 
                if len(s1) + len(s2) > 200:
                    logger.warning ( "The values are too long to print." )
                else:
                    logger.warning( 'The values are: %s (this run) versus %s (default)'%\
                                ( s1,s2 ) )
                return False
    elif isinstance(obj1,list):
        for ival,val in enumerate(obj1):
            if not equalObjs(val,obj2[ival],allowedDiff):
                logger.warning('Lists differ:\n   %s (this run)\n and\n   %s (default)' %\
                                (str(val),str(obj2[ival])))
                return False
    else:
        return obj1 == obj2

    return True


class RunSModelSTest(unittest.TestCase):
    def runMain(self, filename, timeout = 0, suppressStdout=True, development=False,
                 inifile = "testParameters.ini" ):
        to = None
        verbosity='debug'
        if suppressStdout:
            verbosity='error'
            to = os.devnull
        with redirector.stdout_redirected ( to = to ):
            out = join( iDir(), "test/unitTestOutput" )
            main(filename, parameterFile=join ( iDir(), "test/%s" % inifile ),
                 outputDir= out, verbosity = verbosity, db= None, timeout = timeout,
                 development = development)
            sfile = join(iDir(),"test/unitTestOutput/%s.py" % basename(filename))
            return sfile

    def testMultipleFiles( self ):
        out = join( iDir(), "test/unitTestOutput")
        for i in os.listdir( out ):
            if i[-8:]==".smodels":
                os.unlink( os.path.join ( out, i ))
        dirname = join( iDir(), "inputFiles/slha/" )
        self.runMain(dirname)
        nout = len([i for i in glob.iglob("unitTestOutput/*smodels") if not "~" in i])
        nin = len([i for i in glob.iglob("%s/*slha" % dirname) if not "~" in i])
        if nout != nin:
            logger.error("Number of output file (%d) differ from number of input files (%d)" % (nout, nin))
        self.assertTrue( nout == nin )

    def timeoutRun(self):
        filename = join ( iDir(), "inputFiles/slha/complicated.slha" )
        outputfile = self.runMain(filename, timeout=1, suppressStdout=True,
                     development=True, inifile = "timeout.ini" )

    def testTimeout(self):
        self.assertRaises(NoTime, self.timeoutRun)

    def testGoodFile(self):
        filename = join ( iDir(), "inputFiles/slha/gluino_squarks.slha" )
        outputfile = self.runMain(filename)
        shutil.copyfile(outputfile,'./output.py')
        from gluino_squarks_default import smodelsOutputDefault
        from output import smodelsOutput
        ignoreFields = ['input file','smodels version', 'ncpus']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: [res['theory prediction (fb)'],res['TxNames'],
                    res['AnalysisID'],res['DataSetID']])
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields)
        try:
            os.remove('./output.py')
            os.remove('./output.pyc')
        except OSError: pass
        self.assertTrue(equals)

    def testBadFile(self):
        filename = join (iDir(), "inputFiles/slha/I_dont_exist.slha" )
        outputfile = self.runMain(filename)
        shutil.copyfile(outputfile,'./bad_output.py')
        from bad_default import smodelsOutputDefault
        from bad_output import smodelsOutput
        ignoreFields = ['input file','smodels version', 'ncpus']
        equals = equalObjs( smodelsOutput,smodelsOutputDefault,allowedDiff=0.,
                            ignore=ignoreFields)
        os.remove('./bad_output.py')
        os.remove('./bad_output.pyc')
        self.assertTrue(equals )

    def cleanUp ( self ):
        for f in os.listdir("."):
            if ".crash" in f: os.remove(f)
        for i in [ "crash_report_parameter", "crash_report_input" ]:
            if os.path.exists ( i ): os.remove ( i )

    def testCrash(self):
        filename = join ( iDir(), "inputFiles/slha/gluino_squarks.slha" )
        ctr=0
        crash_file = None
        self.cleanUp()
        outputfile = self.runMain(filename, timeout=1, suppressStdout=True,
                                   inifile= "timeout.ini" )
        try:
            ## trying to sync!!
            import ctypes
            libc = ctypes.CDLL("libc.so.6")
            libc.sync()
        except (OSError,AttributeError,ImportError),e:
            pass
        time.sleep(.1)
        for f in os.listdir("."):
            if ".crash" in f:
                crash_file = f
                ctr+=1
        self.assertEquals ( ctr, 1 )
        inp, par = crashReport.readCrashReportFile(crash_file)

        self.assertEquals(open(filename).readlines(), open(inp).readlines())
        self.assertEquals( open("timeout.ini").readlines(),
                           open(par).readlines())
        self.cleanUp()

if __name__ == "__main__":
    unittest.main()
