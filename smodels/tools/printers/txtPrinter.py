"""
.. module:: txtPrinter
   :synopsis: Class for defining a log (stdout) printer

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import os
from smodels.tools.printers.basicPrinter import BasicPrinter
from smodels.decomposition.topologyDict import TopologyDict
from smodels.matching.theoryPrediction import TheoryPredictionList,TheoryPrediction,TheoryPredictionsCombiner
from smodels.experiment.databaseObj import Database
from smodels.tools.ioObjects import OutputStatus
from smodels.tools.coverage import Uncovered
from smodels.base.physicsUnits import GeV, fb, TeV
from smodels.base.smodelsLogging import logger
from smodels.statistics.basicStats import observed, apriori, aposteriori
import numpy as np
from itertools import groupby
import time

class TxTPrinter(BasicPrinter):
    """
    Printer class to handle the printing of one single text output
    """

    def __init__(self, output='stdout', filename=None, outputFormat='current'):
        BasicPrinter.__init__(self, output, filename, outputFormat)
        self.name = "log"
        self.printtimespent = False
        self.printingOrder = [OutputStatus, Database, TopologyDict,
                              TheoryPredictionList, TheoryPredictionsCombiner,
                              TheoryPrediction, Uncovered]
        self.toPrint = [None] * len(self.printingOrder)

    def setOutPutFile(self, filename, overwrite=True, silent=False):
        """
        Set the basename for the text printer. The output filename will be
        filename.log.

        :param filename: Base filename
        :param overwrite: If True and the file already exists, it will be removed.
        :param silent: dont comment removing old files
        """

        self.filename = filename + '.' + self.name
        if overwrite and os.path.isfile(self.filename):
            if not silent:
                logger.warning("Removing old output file " + self.filename)
            os.remove(self.filename)

    def _formatDoc(self, obj):

        return False

    def _formatOutputStatus(self, obj):
        """
        Format data for a OutputStatus object.

        :param obj: A OutputStatus object to be printed.
        """

        output = ""
        output += "Input status: " + str(obj.filestatus) + "\n"
        # hidden feature, printtimespent, turn on in ini file, e.g.
        # [summary-printer] printtimespent = True
        if self.printtimespent:
            output += f"Time spent: {time.time() - self.time:.2f}s\n"
        output += "Decomposition output status: " + str(obj.status) + " "
        st = "unknown status"
        if obj.status in obj.statusStrings:
            st = obj.statusStrings[obj.status]
        output += st + "\n"
        if obj.filestatus < 0:
            output += str(obj.warnings) + "\n"
        output += "# Input File: " + obj.inputfile + "\n"
        labels = list(obj.parameters.keys())
        labels.sort()
        # for label, par in obj.parameters.items():
        for label in labels:
            par = obj.parameters[label]
            output += "# " + label + " = " + str(par) + '\n'
        if obj.smodelsVersion:
            output += f"# SModelS version: {obj.smodelsVersion}\n"
        if obj.databaseVersion:
            output += f"# Database version: {obj.databaseVersion}\n"
        output += "=" * 80 + "\n"
        return output

    def _formatTopologyDict(self, obj):
        """
        Format data for a TopologyDict object.

        :param obj: A TopologyDict object to be printed.
        """

        if not hasattr(self, 'printdecomp') or not self.printdecomp:
            return None
        
        slabel = "Topologies Table"
        output = ""
        output += "  " + "="*56 + "  \n"
        output += "||" + " "*56 + "||\n"
        xspace = int((56-len(slabel))/2.)
        output += "||" + " "*xspace+slabel+" "*(56-xspace-len(slabel))+"||\n"
        output += "||" + " "*56 + "||\n"
        output += "  " + "="*56 + "  \n"

        if self.outputFormat == "version2":
            baseLabel = 'Element'
        else:
            baseLabel = 'SMS'

        # Get topology names:
        topoNames = {canonName : canonName for canonName in obj}
        topoNames_v2 = {}
        if self.outputFormat == "version2":            
            for canonName in obj:
                smsList = obj[canonName]       
                try:
                    sms = smsList[0]
                    evenParticles = sms.treeToBrackets()[0]
                    vertnumb = str([len(v) for v in evenParticles[0]])
                    vertnumb += str([len(v) for v in evenParticles[1]])
                    vertnumb = vertnumb.replace(' ','')
                    topoName = vertnumb
                    topoNames_v2[canonName] = topoName
                except:
                    logger.info("Could not format SMS using version2, switching to current format.")
                    self.outputFormat = 'current'
        
        if self.outputFormat == 'version2':
            topoNames = topoNames_v2

        for canonName in obj:
            smsList = obj[canonName]
            topoName = topoNames[canonName]

            output += "===================================================== \n"
            output += f"Topology: {topoName} \n"
            totxsec = smsList[0].weightList
            for sms in smsList[1:]:
                totxsec += sms.weightList

            output += "Total Global topology weight :\n" + totxsec.niceStr() + '\n'
            output += f"Total Number of {baseLabel}: " + \
                str(len(smsList)) + '\n'
            if not hasattr(self, 'addsmsinfo') or not self.addsmsinfo:
                continue
            for sms in smsList:
                output += "\t\t " + 73 * "." + "\n"
                output += self._formatSMS(sms) + "\n"

        return output

    def _formatSMS(self, obj):
        """
        Format data for an SMS object.

        :param obj: SMS object to be printed.
        """

        output = ""

        if self.outputFormat == 'version2':
            output += "\t\t Element: \n"
            branchList, finalState, _ = obj.treeToBrackets()
            masses = []
            pidlist = []
            for bIndex in obj.daughterIndices(obj.rootIndex):
                branch = obj.indexToNode(bIndex)
                if branch.isSM:
                    continue
                bMasses = [branch.mass]
                pids = [branch.pdg]
                for n in obj.dfsIndexIterator(bIndex):
                    node = obj.indexToNode(n)
                    if node.isSM:
                        continue
                    bMasses.append(node.mass)
                    pids.append(node.pdg)
                masses.append(bMasses)
                pidlist.append(pids)
            # Sort pidlist for consistent output
            pidlist = sorted(pidlist, key = lambda x: str(x))

            output += "\t\t Element ID: " + str(obj.smsID)
            output += "\n"
            output += "\t\t Particles in element: " + str(branchList).replace("'","").replace(" ","")
            output += "\n"
            output += "\t\t Final states in element: " + str(finalState).replace("'","").replace(" ","")
            output += "\n"
            output += "\t\t The element masses are \n"
            for i, mass in enumerate(masses):
                output += "\t\t Branch %i: " % i + str(mass) + "\n"
            output += "\n"
            output += "\t\t The element PIDs are \n"
            pidstr = sorted([str(p) for p in pidlist])
            output += "\t\t PIDs: " + str(pidlist) + "\n"
            output += "\t\t The element weights are: \n \t\t " + \
                obj.weightList.niceStr().replace("\n", "\n \t\t ")

        else:
            output += "\t\t SMS ID: %i\n" %obj.smsID
            output += f"\t\t SMS: {str(obj)}\n"
            output += "\t\t Masses: %s\n" %str([(node,mass) for node,mass in zip(obj.nodes,obj.mass)
                                                if not node.isSM and not node is obj.root])
            output += "\t\t Cross-Sections:\n \t\t "+obj.weightList.niceStr().replace("\n", "\n \t\t ")

        return output

    def _formatDatabase(self, obj):
        """
        Format data for a Database object.

        :param obj: A Database object to be printed.
        """

        if not hasattr(self, "printdatabase") or not self.printdatabase:
            return None

        expResults = obj.expResultList
        slabel = "Selected Experimental Results"
        output = ""
        output += "  " + "="*56 + "  \n"
        output += "||" + " "*56 + "||\n"
        xspace = int((56-len(slabel))/2.)
        output += "||" + " "*xspace+slabel+" "*(56-xspace-len(slabel))+"||\n"
        output += "||" + " "*56 + "||\n"
        output += "  " + "="*56 + "  \n"

        for expRes in expResults:
            output += self._formatExpResult(expRes)

        return output+"\n"

    def _formatExpResult(self, obj):
        """
        Format data for a ExpResult object.

        :param obj: A ExpResult object to be printed.
        """

        if self.outputFormat == "version2":
            baseLabel = 'Elements'
        else:
            baseLabel = 'SMS'


        txnames = []
        for dataset in obj.datasets:
            for txname in dataset.txnameList:
                tx = txname.txName
                if tx not in txnames:
                    txnames.append(tx)

        txnames = sorted(txnames)
        output = ""
        output += "========================================================\n"
        output += "Experimental Result ID: " + obj.globalInfo.id + '\n'
        output += "Tx Labels: " + str(txnames) + '\n'
        output += f"Sqrts: {obj.globalInfo.sqrts.asNumber(TeV):2.2E}\n"
        if hasattr(self, "addanainfo") and self.addanainfo:
            output += "\t -----------------------------\n"
            output += f"\t {baseLabel} tested by analysis:\n"
            listOfSMS = []
            for dataset in obj.datasets:
                for txname in dataset.txnameList:
                    for sms in txname.smsMap:
                        if self.outputFormat == 'version2':
                            smsStr = str(sms.treeToBrackets()[0])
                            smsStr = smsStr.replace("'","").replace(" ","")
                        else:
                            smsStr = str(sms)
                        if smsStr not in listOfSMS:
                            listOfSMS.append(smsStr)
            for sms in listOfSMS:
                output += f"\t    {sms} \n"

        return output

    def _formatNumber(self, number, n=4):
        """ format a number <number> to have n digits,
            but allow also for None, strings, etc """
        if type(number) not in [float, np.float64]:
            return str(number)
        fmt = ".%dg" % n
        return ("%"+fmt) % number

    def _formatTheoryPredictionList(self, obj):
        """
        Format data for a TheoryPredictionList object.

        :param obj: A TheoryPredictionList object to be printed.
        """

        if self.outputFormat == "version2":
            baseLabel = 'Elements'
        else:
            baseLabel = 'SMS'


        slabel = "Theory Predictions and"
        output = ""
        output += "  " + "="*56 + "  \n"
        output += "||" + " "*56 + "||\n"
        xspace = int((56-len(slabel))/2.)
        output += "||" + " "*xspace+slabel+" "*(56-xspace-len(slabel))+"||\n"
        slabel = "Experimental Constraints"
        xspace = int((56-len(slabel))/2.)
        output += "||" + " "*xspace+slabel+" "*(56-xspace-len(slabel))+"||\n"
        output += "||" + " "*56 + "||\n"
        output += "  " + "="*56 + "  \n"

        for theoryPrediction in obj._theoryPredictions:
            expRes = theoryPrediction.expResult
            dataId = theoryPrediction.dataId()
            txnames = [str(txname) for txname in theoryPrediction.txnames]
            txnames = sorted(list(set(txnames)))
            output += "\n"
            output += "---------------Analysis Label = " + expRes.globalInfo.id + "\n"
            output += "-------------------Dataset Label = " + \
                str(dataId).replace("None", "(UL)") + "\n"
            output += "-------------------Txname Labels = " + \
                str(txnames) + "\n"
            output += "Analysis sqrts: " + str(expRes.globalInfo.sqrts) + \
                "\n"

            output += "Theory prediction: " + \
                str(theoryPrediction.xsection) + "\n"
            output += "Theory conditions:"
            output += "  " + str(theoryPrediction.conditions) + "\n"

            # Get upper limit for the respective prediction:
            upperLimit = theoryPrediction.getUpperLimit(evaluationType=observed)
            upperLimitExp = theoryPrediction.getUpperLimit(
                evaluationType=self.getTypeOfExpected())

            output += "Observed experimental limit: " + str(upperLimit) + "\n"
            if upperLimitExp is not None:
                output += "Expected experimental limit: " + \
                    str(upperLimitExp) + "\n"
            srv = self._formatNumber(
                theoryPrediction.getRValue(evaluationType=observed), 4)
            output += f"Observed r-value: {srv}\n"
            if upperLimitExp is not None:
                serv = self._formatNumber(theoryPrediction.getRValue(
                    evaluationType=self.getTypeOfExpected()), 4)
                output += f"Expected r-value: {serv}\n"
            nll = theoryPrediction.likelihood( return_nll = True )
            if nll is not None:
                chi2, chi2sm = None, None
                nllsm = theoryPrediction.lsm( return_nll = True )
                nllmin = theoryPrediction.lmax( return_nll = True )
                try:
                    # chi2sm = -2*np.log(llhd/theoryPrediction.lsm())
                    chi2sm = 2*(nll - nllsm )
                except TypeError:
                    pass
                try:
                    # chi2 = -2*np.log(llhd/theoryPrediction.lmax())
                    chi2 = 2*(nll - nllmin )
                except TypeError:
                    pass
                output += "nll: " + self._formatNumber(nll, 4) + "\n"
                output += "nll_min: " + self._formatNumber(nllmin, 4) + \
                          "   -2log(L/L_max): " + self._formatNumber(chi2, 4) + "\n"
                output += "nll_SM: " + self._formatNumber(nllsm, 4) + \
                          "   -2log(L/L_SM): " + \
                    self._formatNumber(chi2sm, 4) + "\n"

            if hasattr(self, "printextendedresults") and self.printextendedresults:
                if theoryPrediction.mass:
                    for ibr, br in enumerate(theoryPrediction.mass):
                        output += "Masses in branch %i: " % ibr + \
                            str(br) + "\n"
                IDList = list(
                    set([sms.smsID for sms in theoryPrediction.smsList]))
                if IDList:
                    output += f"Contributing {baseLabel}: " + str(IDList) + "\n"

                if self.outputFormat == 'version2':
                    smsPIDs = []
                    for sms in theoryPrediction.smsList:
                        pidList = []
                        for bIndex in sms.daughterIndices(sms.rootIndex):
                            pids = [sms.indexToNode(bIndex).pdg]
                            for nodeIndex in sms.dfsIndexIterator(bIndex):
                                node = sms.indexToNode(nodeIndex)
                                if node.isSM:
                                    continue
                                pids.append(node.pdg)
                            pidList.append(pids)
                        if pidList not in smsPIDs:
                            smsPIDs.append(pidList)
                    for pidList in smsPIDs:
                        # Sort pidlist for consistent output
                        pidList = sorted(pidList, key = lambda x: str(x))
                        output += "PIDs:" + str(pidList) + "\n"

        return output

    def _formatUncovered(self, obj):
        """
        Format all uncovered data.

        :param obj: Uncovered object to be printed.
        """

        if self.outputFormat == "version2":
            baseLabel = 'Element'
        else:
            baseLabel = 'SMS'


        # Number of missing topologies to be printed (ordered by cross sections)
        nprint = 10

        # First sort groups by label
        groups = sorted(obj.groups[:], key=lambda g: g.label)
        # Get summary of groups:
        output = "\n"
        for group in groups:
            output += "Total cross-section for %s (fb): %10.3E\n" % (
                group.description, group.getTotalXSec())

        output += "\n#Full information on unconstrained cross sections\n"
        output += "================================================================================\n"

        # Get detailed information:
        for group in groups:
            description = group.description
            sqrts = group.sqrts.asNumber(TeV)
            if not group.finalStateSMS:
                output += f"No {description} found\n"
                output += "================================================================================\n"
                continue
            output += "%s with the highest cross sections (up to %i):\n" % (
                description, nprint)
            output += f"Sqrts (TeV)   Weight (fb)                  {baseLabel} description\n"
            for fsSMS in group.finalStateSMS[:nprint]:
                if self.outputFormat == 'version2':
                    smsStr = fsSMS.oldStr()
                else:
                    smsStr = str(fsSMS)

                output += "%5s         %10.3E    # %53s \n" % (
                    str(sqrts), fsSMS.missingX, smsStr)
                if hasattr(self, "addcoverageid") and self.addcoverageid:
                    contributing = []
                    for sms in fsSMS._contributingSMS:
                        contributing.append(sms.smsID)
                    output += f"Contributing {baseLabel} {str(contributing)}\n"
            output += "================================================================================\n"
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
        r = self._formatNumber(obj.getRValue(),4)
        r_expected = self._formatNumber(obj.getRValue(evaluationType=self.getTypeOfExpected()),4)
        # Get likelihoods:
        nllsm = obj.lsm( return_nll = True )
        nll = obj.likelihood( return_nll = True )
        nllmin = obj.lmax( return_nll = True )
        output += f"Combined Analyses: {expIDs}\n"
        output += f"Likelihoods: nll, nll_min, nll_SM = {nll:.3f}, {nllmin:.3f}, {nllsm:.3f}\n" 
        output += f"combined r-value: {r:s}\n"
        output += f"combined r-value (expected): {r_expected:s}"
        output += "\n===================================================== \n"
        output += "\n"

        return output
