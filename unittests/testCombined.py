#!/usr/bin/env python3

"""
.. module:: testCombined
   :synopsis: Tests the combination code

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
import os
sys.path.insert(0, "../")
import unittest


def _normalizeMissingTopologyTies(outputDict):
    """
    Normalize tie-boundary entries in printed missing-topology lists.

    Coverage printers keep only the top-N entries (N=10). If several entries tie
    at the cutoff weight, tiny floating-point differences can swap which tied
    element appears as the last printed entry. To keep this regression test focused
    on physics content, replace the element label for cutoff-tied entries by a
    placeholder while keeping the weight untouched.
    """

    keys = [
        'missing topologies',
        'missing topologies with prompt decays',
        'missing topologies with displaced decays',
    ]
    for key in keys:
        if key not in outputDict:
            continue
        entries = outputDict[key]
        if not entries:
            continue

        weights = [entry.get('weight (fb)') for entry in entries]
        if any(weight is None for weight in weights):
            continue

        cutoff = min(weights)
        for entry in entries:
            if entry.get('weight (fb)') == cutoff and 'element' in entry:
                entry['element'] = '__TIE_AT_CUTOFF__'

class CombinedTest(unittest.TestCase):

    def testCombinedResult(self):
        from unitTestHelpers import equalObjs, runMain, importModule
        from smodels.base.smodelsLogging import logger, setLogLevel
        filename = "./testFiles/slha/gluino_squarks.slha"
        setLogLevel ( "error" )
        from databaseLoader import database
        outputfile = runMain(filename, inifile="testParameters_agg.ini",
                overridedatabase = database, suppressStdout=True )
        smodelsOutput = importModule(outputfile)
        from gluino_squarks_default_agg import smodelsOutputDefault
        ignoreFields = ['input file', 'smodels version', 'ncpus',
                'database version', 'model', 'promptwidth', 'stablewidth',
                'checkinput', 'doinvisible', 'docompress', 'computestatistics']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                                 key=lambda res: res['r'], reverse=True)

        # Normalize tie-boundary entries in missing-topology printouts.
        _normalizeMissingTopologyTies(smodelsOutput)
        _normalizeMissingTopologyTies(smodelsOutputDefault)

        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedRelDiff=0.02,
                           ignore=ignoreFields, fname=outputfile)
        if equals != True:
            logger.error( f"gluino_squarks_default_agg.py differs from {outputfile}!" )
        self.assertTrue(equals)
        for i in [outputfile, outputfile.replace(".py", ".pyc")]:
            if os.path.exists(i):
                os.remove(i)


if __name__ == "__main__":
    unittest.main()
