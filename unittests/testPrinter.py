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
from xml.etree import ElementTree
from unitTestHelpers import equalObjs, Summary, runMain, importModule
import pyslha
from smodels.base.smodelsLogging import logger


def sortXML(xmltree):
    for el in xmltree:
        sortXML(el)
    xmltree[:] = sorted(xmltree, key=lambda el: [el.tag, ElementTree.tostring(el)])


def compareXML(xmldefault, xmlnew, allowedRelDiff, ignore=[]):

    if len(xmldefault) != len(xmlnew):
        logger.warning("lengths of document %d != %d" % (len(xmldefault), len(xmlnew)))
        return False
    for i, el in enumerate(xmldefault):
        newel = xmlnew[i]
        if len(el) != len(newel):
            logger.warning("lengths of elements %s and %s differ (%d != %d)" % (el.tag,newel.tag,len(el), len(newel)))
            return False
        if len(el) == 0:
            if el.tag in ignore:
                continue
            if el.tag != newel.tag:
                logger.warning("tags %s and %s differ" % (el.tag, newel.tag))
                return False

            if el.text == newel.text:
                continue

            if type(el.text) == str and "[" not in el.text:
                try:
                    el.text = eval(el.text)
                    newel.text = eval(newel.text)
                except (TypeError, NameError, SyntaxError):
                    el.text = el.text.replace(" ","")
                    newel.text = newel.text.replace(" ","")

            if isinstance(el.text, float) and isinstance(newel.text, float):
                diff = 2.*abs(el.text-newel.text)/abs(el.text+newel.text)
                if diff > allowedRelDiff:
                    logger.warning("values %s and %s differ" % (el.text, newel.text))
                    return False
            elif newel.text != el.text:
                logger.warning("texts %s and %s differ" % (el.text, newel.text))
                return False
        else:
            if not compareXML(el, newel, allowedRelDiff, ignore):
                return False

    return True


def compareSLHA(slhadefault, slhanew, allowedRelDiff):

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
                if de > allowedRelDiff:
                    logger.error(f'Entries in block differ: {ei}!={ej} {type(ei)}')
                    return False
            elif ei != ej:
                logger.error(f'Entries in block differ: {ei}!={ej} {type(ei)}')
                return False
    return True


class RunPrinterTest(unittest.TestCase):
    definingRun = False  # meant only to adapt to changes in output format
    # use with super great care!!

    def removeOutputs(self, f):
        """ remove cruft outputfiles """
        for i in [f, f.replace(".py", ".pyc")]:
            if os.path.exists(i):
                os.remove(i)

    def testTextPrinter(self):

        slhafile = "./testFiles/slha/gluino_squarks.slha"
        outputfile = runMain(slhafile)
        outputfile = outputfile.replace('.py', '.smodels')

        defaultfile = "gluino_squarks_default.smodels"
        # Test summary output
        output = Summary(outputfile, allowedRelDiff=0.05)
        default = Summary(defaultfile, allowedRelDiff=0.05)

        self.assertEqual(default, output)
        self.removeOutputs(outputfile)

    def testPythonPrinter(self):

        slhafile = "./testFiles/slha/gluino_squarks.slha"
        outputfile = runMain(slhafile)
        smodelsOutput = importModule(outputfile)

        # Test python output
        from gluino_squarks_default import smodelsOutputDefault
        ignoreFields = ['input file', 'smodels version', 'ncpus', 'database version']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                                 key=lambda res: res['r'], reverse=True)
        smodelsOutput['ExptRes'] = sorted(smodelsOutput['ExptRes'],
                                          key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedRelDiff=0.05,
                           ignore=ignoreFields, where="top",
                           fname="./unitTestOutput/printer_output.py")

        self.assertTrue(equals)
        self.removeOutputs(outputfile)
        self.removeOutputs('./debug.log')

    def testPythonPrinterSimple(self):

        slhafile = "./testFiles/slha/simplyGluino.slha"
        outputfile = runMain(slhafile, inifile='testParameters_exp.ini')
        smodelsOutput = importModule(outputfile)

        if self.definingRun:
            from smodels.base.smodelsLogging import logger
            logger.error("This is a definition run! Know what youre doing!")
            default = "simplyGluino_default.py"
            outputfile = './unitTestOutput/printer_output_simple.py'
            cmd = "cat %s | sed -e 's/smodelsOutput/smodelsOutputDefault/' > %s" % (outputfile, default)
            subprocess.getoutput(cmd)
        from simplyGluino_default_extended import smodelsOutputDefault

        ignoreFields = ['input file', 'smodels version', 'ncpus', 'database version', 'model', 'checkinput',
                        'doinvisible', 'docompress', 'computestatistics', 'testcoverage']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                                 key=lambda res: res['r'], reverse=True)
        smodelsOutput['ExptRes'] = sorted(smodelsOutput['ExptRes'],
                                          key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedRelDiff=0.05,
                           ignore=ignoreFields)

        self.assertTrue(equals)
        self.removeOutputs(outputfile)

    def testXmlPrinter(self):

        slhafile = "./testFiles/slha/gluino_squarks.slha"
        outputfile = runMain(slhafile)
        outputfile = outputfile.replace('.py', '.xml')

        defFile = "default_output.xml"

        # Test xml output
        xmlDefault = ElementTree.parse(defFile).getroot()
        xmlNew = ElementTree.parse(outputfile).getroot()
        sortXML(xmlDefault)
        sortXML(xmlNew)
        self.assertTrue(compareXML(xmlDefault, xmlNew,
                                   allowedRelDiff=0.05,
                                   ignore=['input_file', 'smodels_version', 'ncpus']))
        # self.removeOutputs(outputfile)
        # self.removeOutputs('./debug.log')

    def testXmlPrinterSimple(self):

        slhafile = "./testFiles/slha/simplyGluino.slha"
        outputfile = runMain(slhafile, inifile='testParameters_exp.ini')
        outputfile = outputfile.replace('.py', '.xml')

        defFile = "default_outputSimplyGluino.xml"

        # Test xml output
        xmlDefault = ElementTree.parse(defFile).getroot()
        xmlNew = ElementTree.parse(outputfile).getroot()
        sortXML(xmlDefault)
        sortXML(xmlNew)

        self.assertTrue(compareXML(xmlDefault, xmlNew, allowedRelDiff=0.05,
                                   ignore=['input_file', 'smodels_version', 'ncpus']))

        self.removeOutputs(outputfile)
        self.removeOutputs('./debug.log')

    def testSLHAPrinter(self):

        slhafile = "./testFiles/slha/gluino_squarks.slha"
        outputfile = runMain(slhafile)
        outputfile = outputfile.replace('.py', '.smodelsslha')

        slhaDefaultFile = "./gluino_squarks_default.slha.smodelsslha"
        self.assertTrue(compareSLHA(slhaDefaultFile, outputfile,
                                    allowedRelDiff=0.05))
        self.removeOutputs(outputfile)


if __name__ == "__main__":
    unittest.main()
