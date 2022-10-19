
"""
.. module:: summaryPrinter
   :synopsis: Class for decribing a summary printer in text format.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import sys
import os
import copy
from smodels.decomposition.topologyDict import TopologyDict
from smodels.matching.theoryPrediction import TheoryPredictionList
from smodels.matching.theoryPredictionsCombiner import TheoryPredictionsCombiner
from smodels.experiment.databaseObj import ExpResultList
from smodels.tools.ioObjects import OutputStatus
from smodels.tools.coverage import Uncovered
from smodels.base.physicsUnits import GeV, fb, TeV
from smodels.base.smodelsLogging import logger
from smodels.tools.printers.txtPrinter import TxTPrinter
import numpy as np
from collections import OrderedDict
from xml.dom import minidom
from xml.etree import ElementTree
import unum
import time


class SummaryPrinter(TxTPrinter):
    """
    Printer class to handle the printing of one single summary output.
    It uses the facilities of the TxTPrinter.
    """

    def __init__(self, output='stdout', filename=None, outputFormat='current'):
        TxTPrinter.__init__(self, output, filename, outputFormat)
        self.name = "summary"
        self.printingOrder = [
            OutputStatus, TheoryPredictionList, TheoryPredictionsCombiner, Uncovered]
        self.toPrint = [None]*len(self.printingOrder)

    def setOutPutFile(self, filename, overwrite=True, silent=False):
        """
        Set the basename for the text printer. The output filename will be
        filename.smodels.
        :param filename: Base filename
        :param overwrite: If True and the file already exists, it will be removed.
        :param silent: dont comment removing old files
        """

        self.filename = filename + '.smodels'
        if overwrite and os.path.isfile(self.filename):
            if not silent:
                logger.warning("Removing old output file " + self.filename)
            os.remove(self.filename)

    def _formatTheoryPredictionList(self, obj):
        """
        Format data of the TheoryPredictionList object.

        :param obj: A TheoryPredictionList object to be printed.
        """
        obj.sortTheoryPredictions()
        if hasattr(self, "expandedsummary") and not self.expandedsummary:
            theoPredictions = [obj._theoryPredictions[0]]
        else:
            theoPredictions = obj._theoryPredictions

        output = ""

        maxr = {"obs": -1., "exp": -1, "anaid": "?"}
        maxcoll = {"CMS": {"obs": -1., "exp": -1, "anaid": "?"},
                   "ATLAS": {"obs": -1., "exp": -1, "anaid": "?"}}
        for theoPred in obj._theoryPredictions:
            r = theoPred.getRValue(expected=False)
            r_expected = theoPred.getRValue(expected=self.getTypeOfExpected())
            expResult = theoPred.expResult
            coll = "ATLAS" if "ATLAS" in expResult.globalInfo.id else "CMS"
            if (r_expected is not None) and (r_expected > maxcoll[coll]["exp"]):
                maxcoll[coll] = {"obs": r, "exp": r_expected,
                                 "anaid": expResult.globalInfo.id}
            if (r is not None) and (r > maxr["obs"]):
                maxr = {"obs": r, "exp": r_expected,
                        "anaid": expResult.globalInfo.id}

        output += "#Analysis  Sqrts  Cond_Violation  Theory_Value(fb)  Exp_limit(fb)  r  r_expected"
        output += "\n\n"
        for theoPred in theoPredictions:
            expResult = theoPred.expResult
            txnames = theoPred.txnames
            ul = theoPred.getUpperLimit(expected=False)
            uls = str(ul)
            if isinstance(ul, unum.Unum):
                uls = "%10.3E" % ul.asNumber(fb)
            signalRegion = theoPred.dataset.getID()
            if signalRegion is None:
                signalRegion = '(UL)'
            value = theoPred.xsection
            r = theoPred.getRValue(expected=False)
            r_expected = theoPred.getRValue(expected=self.getTypeOfExpected())
            rs = str(r)
            rs_expected = str(r_expected)
            if type(r) in [int, float, np.float64]:
                rs = "%10.3E" % r
            if type(r_expected) in [int, float, np.float64]:
                rs_expected = "%10.3E" % r_expected

            output += "%19s  " % (expResult.globalInfo.id)  # ana
            # output += "%4s " % (expResult.globalInfo.sqrts/ TeV)  # sqrts
            # sqrts
            output += "%2.2E  " % (expResult.globalInfo.sqrts.asNumber(TeV))
            output += "%5s " % theoPred.getmaxCondition()  # condition violation
            # theory cross section , expt upper limit
            output += "%10.3E %s " % (value.asNumber(fb), uls)
            if r_expected:
                output += "%s %s" % (rs, rs_expected)
            elif r is None:
                output += "N/A  N/A"
            else:
                output += "%10.3E  N/A" % r
            output += "\n"
            output += " Signal Region:  "+signalRegion+"\n"
            txnameStr = str(sorted(list(set([str(tx) for tx in txnames]))))
            txnameStr = txnameStr.replace(
                "'", "").replace("[", "").replace("]", "")
            output += " Txnames:  " + txnameStr + "\n"
            # print L, L_max and L_SM instead of chi2 and llhd; SK 2021-05-14
            llhd = theoPred.likelihood()
            if llhd is not None:
                lmax = theoPred.lmax()
                lsm = theoPred.lsm()
                lvals = [llhd, lmax, lsm]
                for i, lv in enumerate(lvals):
                    if isinstance(lv, (float, np.float64)):
                        lv = "%10.3E" % lv
                    else:
                        lv = str(lv)
                    lvals[i] = lv
                llhd, lmax, lsm = lvals[:]
                if llhd == lmax == lsm == "None":
                    output += " Likelihoods: L, L_max, L_SM = N/A\n"
                else:
                    output += " Likelihoods: L, L_max, L_SM = %s, %s, %s\n" % (
                        llhd, lmax, lsm)

            if not (theoPred is obj[-1]):
                output += 80 * "-" + "\n"

        output += "\n \n"
        output += 80 * "=" + "\n"
        output += "The highest r value is = %.5f from %s" % \
            (maxr["obs"], maxr["anaid"])
        if maxr["exp"] is not None and maxr["exp"] > -.5:
            output += " (r_expected=%.5f)" % maxr["exp"]
        else:
            output += " (r_expected not available)"
        output += "\n"
        for coll, values in maxcoll.items():
            if values["obs"] < -.5:
                continue
            output += "%s analysis with highest available r_expected: %s, r_expected=%.5f, r_obs=%.5f\n" % \
                      (coll, values["anaid"], values["exp"], values["obs"])

        return output

    def _formatTheoryPredictionsCombiner(self, obj):
        """
        Format data of the TheoryPredictionsCombiner object.

        :param obj: A TheoryPredictionsCombiner object to be printed.
        """

        output = "===================================================== \n"

        # Get list of analyses used in combination:
        expIDs = obj.analysisId()
        # Get r-value:
        r = obj.getRValue()
        r_expected = obj.getRValue(expected=True)
        # Get likelihoods:
        lsm = obj.lsm()
        llhd = obj.likelihood()
        lmax = obj.lmax()
        output += "Combined Analyses: %s\n" % (expIDs)
        output += "Likelihoods: L, L_max, L_SM = %10.3E, %10.3E, %10.3E\n" % (
            llhd, lmax, lsm)
        output += "combined r-value: %10.3E\n" % r
        output += "combined r-value (expected): %10.3E" % r_expected
        output += "\n===================================================== \n"
        output += "\n"

        return output
