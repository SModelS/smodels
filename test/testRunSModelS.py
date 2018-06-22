#!/usr/bin/env python3
 
"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS
 
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
 
"""
 
import sys,os,imp
sys.path.insert(0,"../")
import unittest
import glob
from os.path import join, basename
from smodels.installation import installDirectory as iDir
from smodels.tools import crashReport
from smodels.tools.timeOut import NoTime
from databaseLoader import database ## to make sure the db exists
from smodels.tools.runSModelS import run
import redirector
import unum
import time
 
from smodels.tools.smodelsLogging import logger, setLogLevel
 
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
    if type(obj1) in [ float, int ] and type ( obj2) in [ float, int ]:
        obj1,obj2=float(obj1),float(obj2)
 
    if type(obj1) != type(obj2):
        logger.warning("Data types differ (%s,%s)" %(type(obj1),type(obj2)))
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
                if False: # len(s1) + len(s2) > 200:
                    logger.warning ( "The values are too long to print." )
                else:
                    logger.warning( 'The values are: >>%s<< (this run) versus >>%s<< (default)'%\
                                ( s1[:20],s2[:20] ) )
                return False
    elif isinstance(obj1,list):
        if len(obj1) != len(obj2):
            logger.warning('Lists differ in length:\n   %i (this run)\n and\n   %i (default)' %\
                                (len(obj1),len(obj2)))
            return False
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
        level = 'debug'
        if suppressStdout:
            level = 'error'
            to = os.devnull
        with redirector.stdout_redirected ( to = to ):
            out = join( iDir(), "test/unitTestOutput" )
            setLogLevel ( level )
            run(filename, parameterFile=join ( iDir(), "test/%s" % inifile ),
                 outputDir= out, db= database, timeout = timeout,
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

    def removeOutputs ( self, f ):
        """ remove cruft outputfiles """
        for i in [ f, f.replace(".py",".pyc") ]:
            if os.path.exists ( i ): os.remove ( i )
  
    def testGoodFile(self):
        filename = join ( iDir(), "inputFiles/slha/gluino_squarks.slha" )
        outputfile = self.runMain(filename)
        with open( outputfile, 'rb') as fp: ## imports file with dots in name
            output_module = imp.load_module("output",fp,outputfile, ('.py', 'rb', imp.PY_SOURCE) )
            smodelsOutput = output_module.smodelsOutput
        from gluino_squarks_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'database version']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: [res['theory prediction (fb)'],res['TxNames'],
                    res['AnalysisID'],res['DataSetID']])
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields)
        self.assertTrue(equals)
        self.removeOutputs ( outputfile )
      
    def testGoodFile13(self):
            
        filename = join ( iDir(), "inputFiles/slha/simplyGluino.slha" )
        outputfile = self.runMain(filename,suppressStdout = True )
        with open( outputfile, 'rb') as fp: ## imports file with dots in name
            output_module = imp.load_module("output",fp,outputfile, ('.py', 'rb', imp.PY_SOURCE) )
            smodelsOutput = output_module.smodelsOutput
        from simplyGluino_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'Element', 'database version' ]
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: [res['theory prediction (fb)'],res['TxNames'],
                    res['AnalysisID'],res['DataSetID']])
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields)
        if not equals:
            e =  "output13.py and simplyGluino_default.py differ!" 
            logger.error ( e )
            # raise AssertionError ( e )

        self.assertTrue(equals)        

        ## test went through, so remove the output files
        self.removeOutputs ( outputfile )
    
    def testBadFile(self):
        # since 112 we skip non-existing slha files!
        filename = join (iDir(), "inputFiles/slha/I_dont_exist.slha" )
        of="unitTestOutput/I_dont_exist.slha.py"
        self.removeOutputs ( of )
        outputfile = self.runMain(filename  )
        self.assertTrue ( of in outputfile )
        self.assertTrue ( not os.path.exists ( outputfile ) )
    
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
        except (OSError,AttributeError,ImportError) as e:
            pass
        time.sleep(.1)
        for f in os.listdir("."):
            if ".crash" in f:
                crash_file = f
                ctr+=1
        self.assertEqual ( ctr, 1 )
        inp, par = crashReport.readCrashReportFile(crash_file)
     
        with open(filename) as f:
            with open(inp) as g:
                self.assertEqual(f.readlines(), g.readlines())
        with open("timeout.ini") as f:
            with open(par) as g:
               self.assertEqual( f.readlines(), g.readlines())
        self.cleanUp()
 
if __name__ == "__main__":
    unittest.main()
