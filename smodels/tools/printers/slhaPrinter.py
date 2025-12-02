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
from smodels.base.physicsUnits import fb
from smodels.base.smodelsLogging import logger
from smodels.tools.printers.txtPrinter import TxTPrinter
from smodels.statistics.basicStats import observed
import numpy as np
import unum


class SLHAPrinter(TxTPrinter):
    """
    Printer class to handle the printing of slha format summary output.
    It uses the facilities of the TxTPrinter.
    """

    def __init__(self, output='file', filename=None, outputFormat='version3'):
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

        keysDict = {0: f"{smodelsversion:<25} #SModelS version\n",
                    1: f"{obj.databaseVersion.replace(' ', ''):<25} #database version\n",
                    2: f"{obj.parameters['maxcond']:<25} #maximum condition violation\n",
                    3: f"{self.docompress:<25} #compression (0 off, 1 on)\n",
                    4: f"{obj.parameters['minmassgap']:<25} #minimum mass gap for mass compression [GeV]\n",
                    5: f"{obj.parameters['sigmacut']:<25} #sigmacut [fb]\n",
                    6: f"{self.combinesr:<25} #signal region combination (0 off, 1 on)\n",
                    7: f"{self.combineanas:<25} #analyses combination (0 off, 1 on)\n"}

        if 'promptwidth' in obj.parameters:
            keysDict[8] = f"{obj.parameters['promptwidth']:<25} #prompt width [GeV] \n"
        if 'stablewidth' in obj.parameters:
            keysDict[9] = f"{obj.parameters['stablewidth']:<25} #stable width [GeV] \n"
        if 'minmassgapisr' in obj.parameters:
            keysDict[10] = f"{obj.parameters['minmassgapisr']:<25} #minimum mass gap for ISR mass compression [GeV]\n"

        output = "BLOCK SModelS_Settings\n"
        for key in sorted(list(keysDict.keys())):
            output += f" {key:<2d} {keysDict[key]}"
        output += '\n'

        # for SLHA output we always want to have SModelS_Exclusion block, if no results we write it here
        if obj.status <= 0:
            output += "BLOCK SModelS_Exclusion\n"
            output += f" 0 0 {-1:<30d} #output status (-1 not tested, 0 not excluded, 1 excluded)\n\n"

        return output

    def _formatTheoryPredictionList(self, obj):

        printAll = True  # Controls which theory predictions are printed
        if hasattr(self, "expandedoutput") and not self.expandedoutput:
            printAll = False

        output = "BLOCK SModelS_Exclusion\n"
        if not obj._theoryPredictions[0]:
            excluded = -1
            firstResult = None
        else:
            obj.sortTheoryPredictions()
            firstResult = obj._theoryPredictions[0]
            r = firstResult.getRValue()
            excluded = 0
            if r is not None and r > 1:
                excluded = 1
        output += f" 0 0 {excluded:<30d} #output status (-1 not tested, 0 not excluded, 1 excluded)\n"
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
            # Get sorted txnames
            txWeightsDict = theoPred.getTxNamesWeights(sort=True)
            txnames = []
            for tx in txWeightsDict:
                if tx.txName not in txnames:
                    txnames.append(tx.txName)
            txnameStr = ', '.join(txnames)
            txnameStr = txnameStr.replace("'", "").replace("[", "").replace("]", "")

            signalRegion = theoPred.dataId()
            if signalRegion is None:
                signalRegion = '(UL)'
            r = theoPred.getRValue()
            r_expected = theoPred.getRValue(evaluationType=self.getTypeOfExpected())
            
            
            # Get TxNames final states:
            fStates = []
            for txname in txWeightsDict:
                for sms in txname.smsMap:
                    fs = sms.getFinalStateStr()
                    if fs not in fStates:
                        fStates.append(fs)

            max_length = 3
            fStates_str = ', '.join(fStates[:max_length])
            if len(fStates) > max_length:
                fStates_str += f',...{len(fStates)-max_length:d} more'



            output += f" {cter} 0 {txnameStr:<30} #txname (final states = {fStates_str})\n"
            if r is not None:
                output += f" {cter} 1 {r:<30.3E} #r value\n"
            else:
                output += f" {cter} 1 {"NaN":<30} #r value (failed to compute r-value)\n"
            if not r_expected:
                output += f" {cter} 2 {"N/A":<30} #expected r value\n"  # r_expected could fail or simply not be available
            else:
                output += f" {cter} 2 {r_expected:<30.3E} #expected r value\n"
            output += f" {cter} 3 {theoPred.getmaxCondition():<30.2f} #condition violation\n"
            output += f" {cter} 4 {expResult.globalInfo.id:<30} #analysis\n"
            output += f" {cter} 5 {signalRegion.replace(" ", "_"):<30} #signal region \n"
            nll = theoPred.likelihood( return_nll  = True )
            if nll is not None:
                nllmin = theoPred.lmax(return_nll=True)
                nllsm = theoPred.lsm( return_nll=True )
                lvals = [nll, nllmin, nllsm]
                for i, lv in enumerate(lvals):
                    if isinstance(lv, (float, np.floating)):
                        lv = f"{lv:1.4E}"
                    else:
                        lv = str(lv)
                    lvals[i] = lv
                nll, nllmin, nllsm = lvals[:]
            else:
                nll, nllmin, nllsm = "N/A","N/A","N/A"
    
            output += f" {cter} 6 {nll:<30} #nll\n"
            output += f" {cter} 7 {nllmin:<30} #nll_min\n"
            output += f" {cter} 8 {nllsm:<30} #nll_SM\n"
            output += "\n"

        return output

    def _formatUncovered(self, obj):

        # First sort groups by label
        groups = sorted(obj.groups[:], key=lambda g: g.label)
        # Get summary of groups:
        output = "\nBLOCK SModelS_Coverage"
        for i, group in enumerate(sorted(groups, key=lambda g: g.label)):
            output += f"\n {i} 0 {group.label:<30}      # {group.description}"
            output += f"\n {i} 1 {group.getTotalXSec():<30.3E}      # {"Total cross-section (fb)"}"
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
                if isinstance(lv, (float, np.floating)):
                    lv = f"{lv:1.4E}"
                else:
                    lv = str(lv)
                lvals[i] = lv
            nll, nllmin, nllsm = lvals[:]

            # Get sorted txnames
            txnames = []
            for tx in obj.getTxNamesWeights(sort=True):
                if tx.txName not in txnames:
                    txnames.append(tx.txName)

            if r is not None:
                output += f" {cter} 1 {r:<30.3E} #r value\n"
            else:
                output += f" {cter} 1 {"NaN":<30} #r value (failed to compute r-value)\n"
            if r_expected is not None:
                output += f" {cter} 2 {r_expected:<30.3E} #expected r value\n"
            else:
                output += f" {cter} 2 {"NaN":<30} #expected r value (failed to compute evaluationType r-value)\n"
            output += f" {cter} 3 {nll:<30} #nll\n"
            output += f" {cter} 4 {nllmin:<30} #nll_min\n"
            output += f" {cter} 5 {nllsm:<30} #nll_SM\n"
            output += f" {cter} 6 {expIDs:<30} #IDs of combined analyses\n"
            output += f" {cter} 7 {','.join(txnames):<30} #txnames\n"
            output += "\n"

        return output
