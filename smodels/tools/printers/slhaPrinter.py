"""
.. module:: slhaPrinter
   :synopsis: Class for describing a printer in SLHA format

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import os
from smodels.matching.theoryPrediction import TheoryPredictionList,TheoryPrediction,TheoryPredictionsCombiner
from smodels.tools.ioObjects import OutputStatus
from smodels.tools.coverage import Uncovered
from smodels.base.physicsUnits import GeV, fb, TeV
from smodels.base.smodelsLogging import logger
from smodels.tools.printers.txtPrinter import TxTPrinter
from smodels.statistics.basicStats import observed, apriori, aposteriori
import numpy as np
import unum


class SLHAPrinter(TxTPrinter):
    """
    Printer class to handle the printing of slha format summary output.
    It uses the facilities of the TxTPrinter.
    """

    def __init__(self, output='file', filename=None, outputFormat='current'):
        TxTPrinter.__init__(self, output, filename, outputFormat)
        self.name = "slha"
        self.docompress = 0
        self.combinesr = 0
        self.combineanas = 0
        self.printingOrder = [OutputStatus, TheoryPredictionList,
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

        self.filename = filename + '.smodelsslha'
        if overwrite and os.path.isfile(self.filename):
            if not silent:
                logger.warning("Removing old output file " + self.filename)
            os.remove(self.filename)

    def _formatOutputStatus(self, obj):

        smodelsversion = obj.smodelsVersion
        if not smodelsversion.startswith("v"):
            smodelsversion = "v" + smodelsversion

        keysDict = {0: f"{smodelsversion:25} #SModelS version\n",
                    1: f"{obj.databaseVersion.replace(' ', ''):25} #database version\n",
                    2: f"{obj.parameters['maxcond']:25} #maximum condition violation\n",
                    3: f"{self.docompress:25} #compression (0 off, 1 on)\n",
                    4: f"{obj.parameters['minmassgap']:25} #minimum mass gap for mass compression [GeV]\n",
                    5: f"{obj.parameters['sigmacut']:25} #sigmacut [fb]\n",
                    6: f"{self.combinesr:25} #signal region combination (0 off, 1 on)\n",
                    7: f"{self.combineanas:25} #analyses combination (0 off, 1 on)\n"}

        if 'promptwidth' in obj.parameters:
            keysDict[8] = f"{obj.parameters['promptwidth']:25} #prompt width [GeV] \n"
        if 'stablewidth' in obj.parameters:
            keysDict[9] = f"{obj.parameters['stablewidth']:25} #stable width [GeV] \n"
        if 'minmassgapISR' in obj.parameters:
            keysDict[10] = f"{obj.parameters['minmassgapISR']:25} #minimum mass gap for ISR mass compression [GeV]\n"

        output = "BLOCK SModelS_Settings\n"
        for key in sorted(list(keysDict.keys())):
            output += " %i %s" % (key, keysDict[key])
        output += '\n'

        # for SLHA output we always want to have SModelS_Exclusion block, if no results we write it here
        if obj.status <= 0:
            output += "BLOCK SModelS_Exclusion\n"
            output += " 0 0 %-30s #output status (-1 not tested, 0 not excluded, 1 excluded)\n\n" % (-1)

        return output

    def _formatTheoryPredictionList(self, obj):

        printAll = True  # Controls which theory predictions are printed
        if hasattr(self, "expandedoutput") and not self.expandedoutput:
            printAll = False

        output = "BLOCK SModelS_Exclusion\n"
        if not obj._theoryPredictions[0]:
            excluded = -1
        else:
            obj.sortTheoryPredictions()
            firstResult = obj._theoryPredictions[0]
            r = firstResult.getRValue()
            excluded = 0
            if r is not None and r > 1:
                excluded = 1
        output += " 0 0 %-30s #output status (-1 not tested, 0 not excluded, 1 excluded)\n" % (
            excluded)
        if excluded == -1:
            rList = []
        elif not printAll:
            rList = [firstResult] + [res for res in obj._theoryPredictions[1:]
                                   if (res.getRValue() is not None 
                                       and res.getRValue() >= 1.0)]
        else:
            rList = obj._theoryPredictions[:]

        for iTP, theoPred in enumerate(rList):
            cter = iTP + 1
            expResult = theoPred.expResult
            txnames = theoPred.txnames
            signalRegion = theoPred.dataId()
            if signalRegion is None:
                signalRegion = '(UL)'
            r = theoPred.getRValue()
            r_expected = theoPred.getRValue(evaluationType=self.getTypeOfExpected())
            txnameStr = str(sorted(list(set([str(tx) for tx in txnames]))))
            txnameStr = txnameStr.replace(
                "'", "").replace("[", "").replace("]", "")

            output += " %d 0 %-30s #txname \n" % (cter, txnameStr)
            if r is not None:
                output += " %d 1 %-30.3E #r value\n" % (cter, r)
            else:
                output += " %d 1 NaN                            #r value (failed to compute r-value)\n" % (cter)
            if not r_expected:
                output += " %d 2 N/A                            #expected r value\n" % (cter) # r_expected could fail or simply not be available
            else:
                output += " %d 2 %-30.3E #expected r value\n" % (
                    cter, r_expected)
            output += " %d 3 %-30.2f #condition violation\n" % (
                cter, theoPred.getmaxCondition())
            output += " %d 4 %-30s #analysis\n" % (cter,
                                                   expResult.globalInfo.id)
            output += " %d 5 %-30s #signal region \n" % (
                cter, signalRegion.replace(" ", "_"))
            nll = theoPred.likelihood( return_nll  = True )
            if nll is not None:
                nllmin = theoPred.lmax(return_nll=True)
                nllsm = theoPred.lsm( return_nll=True )
                lvals = [nll, nllmin, nllsm]
                for i, lv in enumerate(lvals):
                    if isinstance(lv, (float, np.float64)):
                        lv = f"{lv:-30.2E}"
                    else:
                        lv = str(lv)
                    lvals[i] = lv
                nll, nllmin, nllsm = lvals[:]
                output += " %d 6 %s #nll\n" % (cter, nll)
                output += " %d 7 %s #nll_min\n" % (cter, nllmin)
                output += " %d 8 %s #nll_SM\n" % (cter, nllsm)
            else:
                output += " %d 6 N/A                            #nll\n" % (
                    cter)
                output += " %d 7 N/A                            #nll_min\n" % (
                    cter)
                output += " %d 8 N/A                            #nll_SM\n" % (
                    cter)
            output += "\n"

        return output

    def _formatUncovered(self, obj):

        # First sort groups by label
        groups = sorted(obj.groups[:], key=lambda g: g.label)
        # Get summary of groups:
        output = "\nBLOCK SModelS_Coverage"
        for i, group in enumerate(sorted(groups, key=lambda g: g.label)):
            output += "\n %d 0 %-30s      # %s" % (
                i, group.label, group.description)
            output += "\n %d 1 %-30.3E      # %s" % (
                i, group.getTotalXSec(), "Total cross-section (fb)")
        output += "\n"
        return output

    def _formatTheoryPrediction(self,obj):
        return self._formatTheoryPredictionsCombiner(obj)

    def _formatTheoryPredictionsCombiner(self, obj):
        """
        Format data of the TheoryPredictionsCombiner object.

        :param obj: A TheoryPredictionsCombiner object to be printed.
        """

        output = "BLOCK SModelS_CombinedAnas\n"

        combRes = [obj]  # For now use a dummy list (only a single combined result is expected)
        for icomb, cRes in enumerate(combRes):
            cter = icomb + 1
            # Get list of analyses IDs used in combination:
            expIDs = cRes.analysisId()
            ul = cRes.getUpperLimit()
            ulExpected = cRes.getUpperLimit(evaluationType=self.getTypeOfExpected())
            if isinstance(ul, unum.Unum):
                ul = ul.asNumber(fb)
            if isinstance(ulExpected, unum.Unum):
                ulExpected = ulExpected.asNumber(fb)

            r = self._round(cRes.getRValue(evaluationType=observed))
            r_expected = self._round(cRes.getRValue(evaluationType=self.getTypeOfExpected()))

            nll = cRes.likelihood(return_nll=True)
            nllmin = cRes.lmax(return_nll=True)
            nllsm = cRes.lsm(return_nll=True)
            lvals = [nll, nllmin, nllsm]
            for i, lv in enumerate(lvals):
                if isinstance(lv, (float, np.float64)):
                    lv = f"{lv:-30.2E}"
                else:
                    lv = str(lv)
                lvals[i] = lv
            nll, nllmin, nllsm = lvals[:]

            if r is not None:
                output += " %d 1 %-30.3E #r value\n" % (cter, r)
            else:
                output += " %d 1 NaN                            #r value (failed to compute r-value)\n" % (cter)
            if r_expected is not None:
                output += " %d 2 %-30.3E #expected r value\n" % (cter, r_expected)
            else:
                output += " %d 2 NaN                            #expected r value (failed to compute evaluationType r-value)\n" % (cter)
            output += " %d 3 %s #nll\n" % (cter, nll)
            output += " %d 4 %s #nll_min\n" % (cter, nllmin)
            output += " %d 5 %s #nll_SM\n" % (cter, nllsm)
            output += " %d 6 %s #IDs of combined analyses\n" % (cter, expIDs)
            output += "\n"

        return output
