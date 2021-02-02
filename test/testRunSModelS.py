#!/usr/bin/env python3

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys,os
sys.path.insert(0,"../")
import unittest
import glob
from smodels.tools import crashReport
from smodels.tools.timeOut import NoTime
from unitTestHelpers import equalObjs, runMain, importModule
import time
import subprocess

from smodels.tools.smodelsLogging import logger

class RunSModelSTest(unittest.TestCase):
    definingRun = False ## meant only to adapt to changes in output format
    ## use with super great care!!

    def testMultipleFiles( self ):
        out = "./unitTestOutput"
        for i in os.listdir( out ):
            if i[-8:]==".smodels":
                os.unlink(os.path.join(out, i))
        dirname = "./testFiles/slha/"
        runMain(dirname)
        nout = len([i for i in glob.iglob("unitTestOutput/*smodels") if not "~" in i])
        nin = len([i for i in glob.iglob("%s/*slha" % dirname) if not "~" in i])
        if nout != nin:
            logger.error("Number of output file(%d) differ from number of input files(%d)" %
                          (nout, nin))
        self.assertEqual(nout,nin)

    def testTimeout(self):
        try:
            filename = "./testFiles/slha/complicated.slha"
            t0=time.time()
            runMain(filename, timeout=1, suppressStdout=True,
                                 development=True, inifile = "timeout.ini" )
            print ( "should never get here. time spent:%.1fs " % ( time.time()-t0 ) )
            self.assertTrue ( False )
        except NoTime:
            self.assertTrue  ( True )
        except Exception as e:
            print ( "wrong exception %s %s" % ( type(e), e ) )
            self.assertTrue ( False )

    def removeOutputs( self, f ):
        """ remove cruft outputfiles """
        for i in [ f, f.replace(".py",".pyc") ]:
            if os.path.exists( i ): os.remove( i )

    def testGoodFile(self):
        filename = "./testFiles/slha/gluino_squarks.slha"
        outputfile = runMain(filename)
        if self.definingRun:
            logger.error ( "This is a definition run! Know what youre doing!" )
            default = "gluino_squarks_default.py"
            cmd = "cat %s | sed -e 's/smodelsOutput/smodelsOutputDefault/' > %s" % ( outputfile, default )
            a = subprocess.getoutput ( cmd )
        smodelsOutput = importModule ( outputfile )
        from gluino_squarks_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'Element', 'database version', 'Total missed xsec',
                            'Missed xsec long-lived', 'Missed xsec displaced', 'Missed xsec MET', 'Total outside grid xsec',
                            'Total xsec for missing topologies (fb)','Total xsec for missing topologies with displaced decays (fb)',
                            'Total xsec for missing topologies with prompt decays (fb)',
                            'Total xsec for topologies outside the grid (fb)']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields, fname = outputfile )
        for i in [ './output.py', './output.pyc' ]:
            if os.path.exists( i ): os.remove( i )
        if not equals:
            p = outputfile.find ( "unitTestOutput" )
            fname = outputfile
            if p > 0:
                fname = fname[p:]
            print ( "[testRunSModelS] %s != %s" % \
                    ( fname, "gluino_squarks_default.py" ) )
        self.assertTrue(equals)
        self.removeOutputs( outputfile )


    def testGoodFileWithModelFromSLHA(self):
        filename = "./testFiles/slha/gluino_squarks.slha"
        outputfile = runMain(filename,inifile='testParametersB.ini')
        smodelsOutput = importModule ( outputfile )
        from gluino_squarks_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'Element', 'database version', 'Total missed xsec',
                            'Missed xsec long-lived', 'Missed xsec displaced', 'Missed xsec MET', 'Total outside grid xsec',
                            'Total xsec for missing topologies (fb)','Total xsec for missing topologies with displaced decays (fb)',
                            'Total xsec for missing topologies with prompt decays (fb)',
                            'Total xsec for topologies outside the grid (fb)']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields, fname = outputfile )
        for i in [ './output.py', './output.pyc' ]:
            if os.path.exists( i ): os.remove( i )
        if not equals:
            p = outputfile.find ( "unitTestOutput" )
            fname = outputfile
            if p > 0:
                fname = fname[p:]
            print ( "[testRunSModelS] %s != %s" % \
                    ( fname, "gluino_squarks_default.py" ) )
        self.assertTrue(equals)
        self.removeOutputs( outputfile )

    def testPyhfCombination(self):
        filename = "./testFiles/slha/T6bbHH_pyhf.slha"
        inifile = "./testParameters_pyhf.ini"
        outputfile = runMain(filename, inifile=inifile, suppressStdout = True)
        smodelsOutput = importModule ( outputfile )
        from T6bbHH_pyhf_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'database version']
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields, fname=outputfile,
	                         fname2 = "T6bbHH_pyhf_default.py" )
        self.assertTrue(equals)

    def testGoodFile13(self):

        filename = "./testFiles/slha/simplyGluino.slha"
        outputfile = runMain(filename,suppressStdout = True )
        if self.definingRun:
            logger.error ( "This is a definition run! Know what youre doing!" )
            default = "simplyGluino_default.py"
            cmd = "cat %s | sed -e 's/smodelsOutput/smodelsOutputDefault/' > %s" % ( outputfile, default )
            a = subprocess.getoutput ( cmd )
        smodelsOutput = importModule ( outputfile )
        from simplyGluino_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'Element', 'database version', 'Total missed xsec',
                            'Missed xsec long-lived', 'Missed xsec displaced', 'Missed xsec MET', 'Total outside grid xsec',
                            'Total xsec for missing topologies (fb)','Total xsec for missing topologies with displaced decays (fb)',
                            'Total xsec for missing topologies with prompt decays (fb)',
                            'Total xsec for topologies outside the grid (fb)']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.08,
                           ignore=ignoreFields, fname = outputfile )
        if not equals:
            e =  "output13.py != simplyGluino_default.py"
            logger.error( e )
            # raise AssertionError( e )

        self.assertTrue(equals)

        ## test went through, so remove the output files
        self.removeOutputs( outputfile )

    def testGoodFileHSCP(self):
        filename = "./testFiles/slha/longLived.slha"
        outputfile = runMain(filename)
        if self.definingRun:
            logger.error ( "This is a definition run! Know what youre doing!" )
            default = "longLived_default.py"
            cmd = "cat %s | sed -e 's/smodelsOutput/smodelsOutputDefault/' > %s" % ( outputfile, default )
            a = subprocess.getoutput ( cmd )
        smodelsOutput = importModule ( outputfile )
        from longLived_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'Element', 'database version', 'Total missed xsec',
                        'Missed xsec long-lived', 'Missed xsec displaced', 'Missed xsec MET', 'Total outside grid xsec',
                        'Total xsec for missing topologies (fb)','Total xsec for missing topologies with displaced decays (fb)',
                        'Total xsec for missing topologies with prompt decays (fb)',
                        'Total xsec for topologies outside the grid (fb)']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields, fname = outputfile )
        for i in [ './outputHSCP.py', './outputHSCP.pyc' ]:
            if os.path.exists( i ): os.remove( i )
        self.assertTrue(equals)

    def testLifeTimeDependent(self):
        filename = "./testFiles/slha/lifetime.slha"
        outputfile = runMain(filename)
        if self.definingRun:
            logger.error ( "This is a definition run! Know what youre doing!" )
            default = "lifetime_default.py"
            cmd = "cat %s | sed -e 's/smodelsOutput/smodelsOutputDefault/' > %s" % ( outputfile, default )
            a = subprocess.getoutput ( cmd )
        smodelsOutput = importModule ( outputfile )
        from lifetime_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'database version']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields, fname=outputfile )
        for i in [ './outputHSCP.py', './outputHSCP.pyc' ]:
            if os.path.exists( i ): os.remove( i )
        self.assertTrue(equals)


    def testBadFile(self):
        # since 112 we skip non-existing slha files!
        filename = "./testFiles/slha/I_dont_exist.slha"
        of="unitTestOutput/I_dont_exist.slha.py"
        self.removeOutputs(of)
        outputfile = runMain(filename)
        self.assertTrue( of in outputfile )
        self.assertTrue( not os.path.exists( outputfile ) )

    def cleanUp( self ):
        for f in os.listdir("."):
            if ".crash" in f: os.remove(f)
        for i in [ "crash_report_parameter", "crash_report_input" ]:
            if os.path.exists( i ): os.remove( i )

    def testCrash(self):
        filename = "./testFiles/slha/gluino_squarks.slha"
        ctr=0
        crash_file = None
        self.cleanUp()
        runMain(filename, timeout=1, suppressStdout=True,
                                   inifile= "timeout.ini" )
        """
        try:
            ## trying to sync!!
            import ctypes
            libc = ctypes.CDLL("libc.so.6")
            libc.sync()
        except(OSError,AttributeError,ImportError) as e:
            print ( "This shouldnt throw %s" % e )
            # pass
        """
        time.sleep(.2)
        for f in os.listdir("."):
            if ".crash" in f:
                crash_file = f
                ctr+=1
        self.assertEqual( ctr, 1 )
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
