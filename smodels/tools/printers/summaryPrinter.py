"""
.. module:: summaryPrinter
   :synopsis: Class for describing a summary printer in text format.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import os
from smodels.matching.theoryPrediction import TheoryPredictionList,TheoryPrediction,TheoryPredictionsCombiner
from smodels.tools.ioObjects import OutputStatus
from smodels.tools.coverage import Uncovered
from smodels.base.physicsUnits import fb, TeV
from smodels.base.smodelsLogging import logger
from smodels.tools.printers.txtPrinter import TxTPrinter
import numpy as np
import unum


class SummaryPrinter(TxTPrinter):
    """
    Printer class to handle the printing of one single summary output.
    It uses the facilities of the TxTPrinter.
    """

    def __init__(self, output='stdout', filename=None, outputFormat='current'):
        TxTPrinter.__init__(self, output, filename, outputFormat)
        self.name = "summary"
        self.printingOrder = [
            OutputStatus, TheoryPredictionList, 
            TheoryPredictionsCombiner, 
            TheoryPrediction, Uncovered]
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
            r = theoPred.getRValue(evaluationType=False)
            r_expected = theoPred.getRValue(evaluationType=self.getTypeOfExpected())
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
            ul = theoPred.getUpperLimit(evaluationType=False)
            uls = str(ul)
            if isinstance(ul, unum.Unum):
                uls = f"{ul.asNumber(fb):10.3E}"
            signalRegion = theoPred.dataset.getID()
            if signalRegion is None:
                signalRegion = '(UL)'
            value = theoPred.xsection
            r = theoPred.getRValue(evaluationType=False)
            r_expected = theoPred.getRValue(evaluationType=self.getTypeOfExpected())
            if r is not None:
                rs = f"{r:10.3E}"
            else:
                rs = "NaN" # r = None means the calculation failed
            if r_expected is not None:
                rs_expected = f"{r_expected:10.3E}"
            else:
                rs_expected = "N/A" # r_exp could not be available

            output += "%19s  " % (expResult.globalInfo.id)  # ana
            # output += "%4s " % (expResult.globalInfo.sqrts/ TeV)  # sqrts
            # sqrts
            output += f"{expResult.globalInfo.sqrts.asNumber(TeV):2.2E}  "
            output += "%5s " % theoPred.getmaxCondition()  # condition violation
            # theory cross section , expt upper limit
            output += f"{value.asNumber(fb):10.3E} {uls} "
            output += f"{rs} {rs_expected}"
            
            output += "\n"
            output += " Signal Region:  "+signalRegion+"\n"
            txnameStr = str(sorted(list(set([str(tx) for tx in txnames]))))
            txnameStr = txnameStr.replace(
                "'", "").replace("[", "").replace("]", "")
            output += " Txnames:  " + txnameStr + "\n"
            nll = theoPred.likelihood( return_nll = True )
            if nll is not None:
                nllmin = theoPred.lmax( return_nll = True )
                nllsm = theoPred.lsm( return_nll = True )
                lvals = [nll, nllmin, nllsm]
                for i, lv in enumerate(lvals):
                    if isinstance(lv, (float, np.float64)):
                        lv = f"{lv:10.3E}"
                    else:
                        lv = str(lv)
                    lvals[i] = lv
                nll, nllmin, nllsm = lvals[:]
                if nll == nllmin == nllsm == "None":
                    output += " Likelihoods: nll, nll_min, nll_SM = N/A\n"
                else:
                    output += f" Likelihoods: nll, nll_min, nll_SM = {nll}, {nllmin}, {nllsm}\n"

            if not (theoPred is obj[-1]):
                output += 80 * "-" + "\n"

        output += "\n \n"
        output += 80 * "=" + "\n"
        output += f"The highest r value is = {maxr['obs']:.5f} from {maxr['anaid']}"
        if maxr["exp"] is not None and maxr["exp"] >= 0.0:
            output += f" (r_expected={maxr['exp']:.5f})"
        else:
            output += " (r_expected not available)"
        output += "\n"
        for coll, values in maxcoll.items():
            if values["obs"] == None or values["obs"] < 0.0:
                continue
            output += "%s analysis with highest available r_expected: %s, r_expected=%.5f, r_obs=%.5f\n" % \
                      (coll, values["anaid"], values["exp"], values["obs"])

        return output

    def _formatTheoryPrediction(self,obj):
        return self._formatTheoryPredictionsCombiner(obj)

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
        r_expected = obj.getRValue(evaluationType=self.getTypeOfExpected())
        # Get likelihoods:
        nllsm = obj.lsm( return_nll = True )
        nll = obj.likelihood( return_nll = True )
        nllmin = obj.lmax( return_nll = True )
        output += f"Combined Analyses: {expIDs}\n"
        output += f"Likelihoods: nll, nll_min, nll_SM = {nll:.3f}, {nllmin:.3f}, {nllsm:.3f}\n"
        if r is not None:
            output += f"combined r-value: {r:10.3E}\n"
        else:
            output += f"combined r-value: NaN (failed to compute r-value)\n"
        if r_expected is not None:
            output += f"combined r-value (expected): {r_expected:10.3E}\n"
        else:
            output += f"combined r-value (expected): NaN (failed to compute r-value)\n"
        output += "\n===================================================== \n"
        output += "\n"

        return output
