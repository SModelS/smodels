#!/usr/bin/env python3

"""
.. module:: testPrinter
   :synopsis: Tests printer facilities with runSModelS

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
import os
import importlib
import subprocess
sys.path.insert(0, "../")
import unittest
from smodels.theory import decomposer
from smodels.tools.physicsUnits import fb, GeV
from smodels.theory.theoryPrediction import theoryPredictionsFor, TheoryPredictionList
from smodels.tools import printer, ioObjects
from smodels.tools import coverage
from xml.etree import ElementTree
from databaseLoader import database
from unitTestHelpers import equalObjs, Summary, runMain
from imp import reload
from smodels.tools import runtime
from smodels import particlesLoader
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
import pyslha
from smodels.tools.smodelsLogging import logger


def sortXML(xmltree):
    for el in xmltree:
        sortXML(el)
    xmltree[:] = sorted(xmltree, key=lambda el: [el.tag, ElementTree.tostring(el)])


def log(strng):
    with open("debug.log", "at") as f:
        f.write(strng + "\n")
        f.close()


def compareXML(xmldefault, xmlnew, allowedDiff, ignore=[]):
    if len(xmldefault) != len(xmlnew):
        log("lengths of document %d != %d" % (len(xmldefault), len(xmlnew)))
        return False
    for i, el in enumerate(xmldefault):
        newel = xmlnew[i]
        if len(el) != len(newel):
            log("lengths of elements %d != %d" % (len(el), len(newel)))
            return False
        if len(el) == 0:
            if el.tag in ignore:
                continue
            if type(el.text) == str and not "[" in el.text:
                try:
                    el.text = eval(el.text)
                    newel.text = eval(newel.text)
                except (TypeError, NameError, SyntaxError):
                    pass
            if isinstance(el.text, float) and isinstance(newel.text, float) \
                    and newel.text != el.text:
                diff = 2.*abs(el.text-newel.text)/abs(el.text+newel.text)
                if diff > allowedDiff:
                    log("values %s and %s differ" % (el.text, newel.text))
                    return False
            else:
                if el.text != newel.text:
                    log("texts %s and %s differ" % (el.text, newel.text))
                    return False
            if el.tag != newel.tag:
                log("tags %s and %s differ" % (el.tag, newel.tag))
                return False
        else:
            compareXML(el, newel, allowedDiff, ignore)

    return True


def compareSLHA(slhadefault, slhanew, allowedDiff):

    newData = pyslha.read(slhanew, ignorenomass=True, ignorenobr=True, ignoreblocks=["SMODELS_SETTINGS"])
    defaultData = pyslha.read(slhadefault, ignorenomass=True, ignorenobr=True, ignoreblocks=["SMODELS_SETTINGS"])
    defaultBlocks = sorted([defaultData.blocks[b].name for b in defaultData.blocks])
    newBlocks = sorted([newData.blocks[b].name for b in newData.blocks])
    if defaultBlocks != newBlocks:
        logger.error('Block structure differs!')
        return False

    for b in defaultData.blocks:
        if len(defaultData.blocks[b].entries) != len(newData.blocks[b].entries):
            logger.error('Numbers of entries in block %s differ' % (defaultData.blocks[b].name))
            return False
        keys = defaultData.blocks[b].keys()
        bkeys = newData.blocks[b].keys()
        if keys != bkeys:
            logger.error(f"Keys of blocks {b} differ!")
            return False
        for k in keys:
            ei = defaultData.blocks[b].entries[k]
            ej = newData.blocks[b].entries[k]
            if type(ei) == float and type(ej) == float:
                denom = ei + ej
                if denom == 0.:
                    denom = 1e-6
                de = 2. * abs(ei - ej) / denom
                if de > allowedDiff:
                    logger.error(f'Entries in block differ: {ei}!={ej} {type(ei)}')
                    return False
            elif ei != ej:
                logger.error(f'Entries in block differ: {ei}!={ej} {type(ei)}')
                return False
    return True


class RunPrinterTest(unittest.TestCase):
    definingRun = False  # meant only to adapt to changes in output format
    # use with super great care!!

    def runPrinterMain(self, slhafile, mprinter, addTopList=False):
        """
        Main program. Displays basic use case.

        """
        runtime.modelFile = 'mssm'
        reload(particlesLoader)

        # Set main options for decomposition:
        sigmacut = 0.03 * fb
        mingap = 5. * GeV

        """ Decompose model  """
        model = Model(BSMList, SMList)
        model.updateParticles(slhafile)
        smstoplist = decomposer.decompose(model, sigmacut,
                         doCompress=True, doInvisible=True, minmassgap=mingap)

        # Add the decomposition result to the printers
        if addTopList:
            mprinter.addObj(smstoplist)
        listOfExpRes = database.getExpResults(analysisIDs=['*:8*TeV', 'CMS-PAS-SUS-15-002', 'CMS-PAS-SUS-16-024'])
        # Compute the theory predictions for each analysis
        allPredictions = []
        for expResult in listOfExpRes:
            predictions = theoryPredictionsFor(expResult, smstoplist)
            if not predictions:
                continue
            allPredictions += predictions._theoryPredictions

        for theoPred in allPredictions:
            if theoPred.dataType() == 'efficiencyMap':  # and hasattr(theoPred,'expectedUL') and not theoPred.expectedUL is None:
                theoPred.computeStatistics()

        maxcond = 0.2
        theoryPredictions = TheoryPredictionList(allPredictions, maxcond)
        mprinter.addObj(theoryPredictions)

        # Add coverage information:
        coverageInfo = coverage.Uncovered(smstoplist)
        mprinter.addObj(coverageInfo)

        # Add additional information:
        databaseVersion = database.databaseVersion
        outputStatus = ioObjects.OutputStatus([1, 'Input file ok'], slhafile,
                                              {'sigmacut': sigmacut.asNumber(fb),
                                               'minmassgap': mingap.asNumber(GeV),
                                               'maxcond': maxcond},
                                              databaseVersion)
        outputStatus.status = 1
        mprinter.addObj(outputStatus)
        mprinter.flush()

    def testTextPrinter(self):
        outputfile = "./unitTestOutput/printer_output.smodels"
        self.removeOutputs(outputfile)

        slhafile = "./testFiles/slha/gluino_squarks.slha"
        outputfile = runMain(slhafile)
        outputfile = outputfile.replace('.py', '.smodels')

        defaultfile = "gluino_squarks_default.smodels"
        # Test summary output
        output = Summary(outputfile, allowedDiff=0.05)
        default = Summary(defaultfile, allowedDiff=0.05)
        try:
            self.assertEqual(default, output)
            # self.removeOutputs(outputfile)
        except AssertionError:
            msg = "%s != %s" % (default, output)
            raise AssertionError(msg)

    def removeOutputs(self, f):
        """ remove cruft outputfiles """
        for i in [f, f.replace(".py", ".pyc")]:
            if os.path.exists(i):
                os.remove(i)

    def testPythonPrinter(self):

        self.removeOutputs('./unitTestOutput/printer_output.py')
        mprinter = printer.MPrinter()
        mprinter.Printers['python'] = printer.PyPrinter(output='file')
        # Define options:
        mprinter.Printers['python'].addElementList = False

        slhafile = "./testFiles/slha/gluino_squarks.slha"
        mprinter.setOutPutFiles('./unitTestOutput/printer_output', silent=True)
        self.runPrinterMain(slhafile, mprinter)

        try:
            smodelsOutput = importlib.import_module("unitTestOutput.printer_output").smodelsOutput
        except ImportError:  # Python2 fallback
            import imp
            pM = imp.load_source("smodels", "./unitTestOutput/printer_output.py")
            smodelsOutput = pM.smodelsOutput
        # Test python output
        from gluino_squarks_default import smodelsOutputDefault
        ignoreFields = ['input file', 'smodels version', 'ncpus', 'database version',
                      'Total missed xsec', 'Missed xsec long-lived', 'Missed xsec displaced',
                      'Missed xsec MET', 'Total outside grid xsec',
                      'Total xsec for missing topologies (fb)', 'Total xsec for missing topologies with displaced decays (fb)',
                      'Total xsec for missing topologies with prompt decays (fb)',
                      'Total xsec for topologies outside the grid (fb)']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                               key=lambda res: res['r'], reverse=True)
        smodelsOutput['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                        key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedDiff=0.05,
                         ignore=ignoreFields, where="top",
                         fname="./unitTestOutput/printer_output.py")
        try:
            self.assertTrue(equals)
        except AssertionError as e:
            print("Error: %s, when comparing %s \nwith %s." % (e, "output.py", "gluino_squarks_default.py"))
            raise AssertionError(e)
        self.removeOutputs('./unitTestOutput/printer_output.py')
        self.removeOutputs('./debug.log')

    def testPythonPrinterSimple(self):
        outfile = './unitTestOutput/printer_output_simple.py'
        self.removeOutputs(outfile)

        mprinter = printer.MPrinter()
        mprinter.Printers['python'] = printer.PyPrinter(output='file')
        # Define options:
        mprinter.Printers['python'].addelementlist = True

        slhafile = "./testFiles/slha/simplyGluino.slha"
        mprinter.setOutPutFiles(outfile.replace(".py", ""), silent=True)
        self.runPrinterMain(slhafile, mprinter, addTopList=True)

        if self.definingRun:
            from smodels.tools.smodelsLogging import logger
            logger.error("This is a definition run! Know what youre doing!")
            default = "simplyGluino_default.py"
            outputfile = './unitTestOutput/printer_output_simple.py'
            cmd = "cat %s | sed -e 's/smodelsOutput/smodelsOutputDefault/' > %s" % (outputfile, default)
            subprocess.getoutput(cmd)

        impfile = outfile.replace(".py", "").replace("/", ".").replace("..", "")
        try:
            smodelsOutput = importlib.import_module(impfile).smodelsOutput
        except ImportError:  # Python2 fallback
            import imp
            pM = imp.load_source("smodels", outfile)
            smodelsOutput = pM.smodelsOutput

        from simplyGluino_default_extended import smodelsOutputDefault

        ignoreFields = ['input file', 'smodels version', 'ncpus', 'database version', 'Total missed xsec',
                      'Missed xsec long-lived', 'Missed xsec displaced', 'Missed xsec MET', 'Total outside grid xsec',
                      'Total xsec for missing topologies (fb)', 'Total xsec for missing topologies with displaced decays (fb)',
                      'Total xsec for missing topologies with prompt decays (fb)',
                      'Total xsec for topologies outside the grid (fb)']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                               key=lambda res: res['r'], reverse=True)
        smodelsOutput['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                        key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedDiff=0.05,
                         ignore=ignoreFields, fname="./unitTestOutput/printer_output_simple.py")
        try:
            self.assertTrue(equals)
            self.removeOutputs(outfile)
        except AssertionError as e:
            print("Error: %s, when comparing %s \nwith %s." % (e, "outputSimple.py", "simplyGluino_default.py"))
            raise AssertionError(e)

    def testXmlPrinter(self):
        self.removeOutputs('./unitTestOutput/printer_output.xml')

        mprinter = printer.MPrinter()
        mprinter.Printers['xml'] = printer.XmlPrinter(output='file')
        #Define options:
        mprinter.Printers['xml'].addelementlist = False

        slhafile = "./testFiles/slha/gluino_squarks.slha"
        mprinter.setOutPutFiles('./unitTestOutput/printer_output', silent=True)
        self.runPrinterMain(slhafile, mprinter)

        defFile = "default_output.xml"
        outFile = "./unitTestOutput/printer_output.xml"

        #Test xml output
        xmlDefault = ElementTree.parse(defFile).getroot()
        xmlNew = ElementTree.parse(outFile).getroot()
        sortXML(xmlDefault)
        sortXML(xmlNew)
        try:
            self.assertTrue(compareXML(xmlDefault, xmlNew,
                                       allowedDiff=0.05,
                                       ignore=['input_file', 'smodels_version', 'ncpus',
                                              'Total missed xsec', 'Missed xsec long-lived',
                                              'Missed xsec displaced', 'Missed xsec MET',
                                              'Total outside grid xsec',
                                              'Total xsec for missing topologies (fb)',
                                              'Total xsec for missing topologies with displaced decays (fb)',
                                              'Total xsec for missing topologies with prompt decays (fb)',
                                              'Total xsec for topologies outside the grid (fb)']))
        except AssertionError as e:
            msg = "%s != %s" % (defFile, outFile) + "\n" + str(e)
            msg += ". Try and consult debug.log for more info."
            raise AssertionError(msg)
        self.removeOutputs('./unitTestOutput/printer_output.xml')
        self.removeOutputs('./debug.log')

    def testXmlPrinterSimple(self):
        self.removeOutputs('./unitTestOutput/printer_output_simple.xml')

        mprinter = printer.MPrinter()
        mprinter.Printers['xml'] = printer.XmlPrinter(output='file')
        #Define options:
        mprinter.Printers['xml'].addelementlist = True

        slhafile = "./testFiles/slha/simplyGluino.slha"
        mprinter.setOutPutFiles('./unitTestOutput/printer_output_simple', silent=True)
        self.runPrinterMain(slhafile, mprinter, addTopList=True)

        defFile = "default_outputSimplyGluino.xml"
        outFile = "./unitTestOutput/printer_output_simple.xml"

        #Test xml output
        xmlDefault = ElementTree.parse(defFile).getroot()
        xmlNew = ElementTree.parse(outFile).getroot()
        sortXML(xmlDefault)
        sortXML(xmlNew)
        try:
            self.assertTrue(compareXML(xmlDefault, xmlNew, allowedDiff=0.05,
                                       ignore=['input_file', 'smodels_version', 'ncpus', 'Total missed xsec',
                                               'Missed xsec long-lived', 'Missed xsec displaced',
                                               'Missed xsec MET', 'Total outside grid xsec',
                                               'Total xsec for missing topologies (fb)',
                                               'Total xsec for missing topologies with displaced decays (fb)',
                                               'Total xsec for missing topologies with prompt decays (fb)',
                                               'Total xsec for topologies outside the grid (fb)']))
        except AssertionError as e:
            msg = "%s != %s" % (defFile, outFile) + "\n" + str(e)
            raise AssertionError(msg)
        self.removeOutputs('./unitTestOutput/printer_output_simple.xml')
        self.removeOutputs('./debug.log')

    def testSLHAPrinter(self):
        self.removeOutputs('./unitTestOutput/printer_output.smodelsslha')

        mprinter = printer.MPrinter()
        mprinter.Printers['slha'] = printer.SLHAPrinter(output='file')
        # Define options:
        mprinter.Printers['slha'].addelementlist = True
        mprinter.Printers['slha'].docompress = 1

        slhafile = "./testFiles/slha/gluino_squarks.slha"
        mprinter.setOutPutFiles('./unitTestOutput/printer_output', silent=False)
        self.runPrinterMain(slhafile, mprinter, addTopList=True)

        slhaDefaultFile = "./gluino_squarks_default.slha.smodelsslha"
        slhaNewFile = './unitTestOutput/printer_output.smodelsslha'
        try:
            self.assertTrue(compareSLHA(slhaDefaultFile, slhaNewFile,
                                        allowedDiff=0.05))
        except AssertionError:
            msg = "%s != %s" % (slhaDefaultFile, slhaNewFile)
            raise AssertionError(msg)
        self.removeOutputs('./unitTestOutput/printer_output.smodelsslha')


if __name__ == "__main__":
    unittest.main()
