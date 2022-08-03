
"""
.. module:: pythonPrinter
   :synopsis: Class for decribing a printer in python format

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import sys
import os
import copy
from smodels.theory.topologyDict import TopologyDict
from smodels.matcher.theoryPrediction import TheoryPredictionList
from smodels.tools.theoryPredictionsCombiner import TheoryPredictionsCombiner
from smodels.experiment.databaseObj import ExpResultList
from smodels.tools.ioObjects import OutputStatus
from smodels.tools.coverage import Uncovered
from smodels.tools.physicsUnits import GeV, fb, TeV
from smodels.tools.smodelsLogging import logger
from smodels.tools.printers.basicPrinter import BasicPrinter
import numpy as np
from collections import OrderedDict
from xml.dom import minidom
from xml.etree import ElementTree
import unum
import time


class PyPrinter(BasicPrinter):
    """
    Printer class to handle the printing of one single pythonic output
    """

    def __init__(self, output='stdout', filename=None, outputFormat='current'):
        BasicPrinter.__init__(self, output, filename, outputFormat)
        self.name = "py"
        self.printtimespent = False
        self.printingOrder = [OutputStatus, TopologyDict,
                              TheoryPredictionList, TheoryPredictionsCombiner, Uncovered]
        self.toPrint = [None]*len(self.printingOrder)

    def setOutPutFile(self, filename, overwrite=True, silent=False):
        """
        Set the basename for the text printer. The output filename will be
        filename.py.
        :param filename: Base filename
        :param overwrite: If True and the file already exists, it will be removed.
        :param silent: dont comment removing old files
        """

        self.filename = filename + '.py'
        if overwrite and os.path.isfile(self.filename):
            if not silent:
                logger.warning("Removing old output file " + self.filename)
            os.remove(self.filename)

    def flush(self):
        """
        Write the python dictionaries generated by the object formatting
        to the defined output
        """

        outputDict = {}
        for obj in self.toPrint:
            if obj is None:
                continue
            output = self._formatObj(obj)
            if not output:
                continue  # Skip empty output
            outputDict.update(output)

        output = 'smodelsOutput = '+str(outputDict)
        if self.output == 'stdout':
            sys.stdout.write(output)
        elif self.output == 'file':
            if not self.filename:
                logger.error('Filename not defined for printer')
                return False
            with open(self.filename, "a") as outfile:
                outfile.write(output)
                outfile.close()

        self.toPrint = [None]*len(self.printingOrder)
        # it is a special feature of the python printer
        # that we also return the output dictionary
        return outputDict

    def _formatTopologyDict(self, obj):
        """
        Format data for a TopologyDict object.

        :param obj: A TopologyDict object to be printed.
        """

        if not hasattr(self, 'addsmslist') or not self.addsmslist:
            return None


        smsList = []
        for sms in obj.getSMSList():
            smsList.append(self._formatSMS(sms))

        if self.outputFormat == 'version2':
            return {'Element' : smsList}
        else:
            return {"SMS Decomposition": smsList}

    def _formatSMS(self, obj):
        """
        Format data for a SMS object.

        :param obj: A SMS object to be printed.
        """

        smsDict = {}
        if self.outputFormat == 'version2':
            branchList, finalState, intermediateState = obj.treeToBrackets()
            masses = []
            pidlist = []
            for bIndex in obj.daughterIndices(obj.rootIndex):
                branch = obj.indexToNode(bIndex)
                if branch.isSM:
                    continue
                mass = float('%1.3e' %branch.mass.asNumber(GeV))
                bMasses = [mass]
                pids = [branch.pdg]
                for n in obj.dfsIndexIterator(bIndex):
                    node = obj.indexToNode(n)
                    if node.isSM:
                        continue
                    mass = float('%1.3e' %node.mass.asNumber(GeV))
                    bMasses.append(mass)
                    pids.append(node.pdg)
                masses.append(bMasses)
                pidlist.append(pids)

            smsDict["ID"] = obj.smsID
            smsDict["Particles"] =  str(branchList).replace("'","").replace(" ","")
            smsDict["final states"] =  finalState
            smsDict["Masses (GeV)"] = masses
            smsDict["PIDs"] = pidlist
        else:
            smsDict["ID"] = obj.smsID
            smsDict["SMS"] = str(obj)
            smsDict["Masses (GeV)"] = [(str(node),float('%1.3e' %node.mass.asNumber(GeV)))
                                      for node in obj.nodes if not node.isSM]
            smsDict["PIDs"] = [(str(node),node.pdg)
                              for node in obj.nodes if not node.isSM]


        smsDict["Weights (fb)"] = {}
        sqrts = [info.sqrts.asNumber(TeV) for info in obj.weightList.getInfo()]
        allsqrts = sorted(list(set(sqrts)))
        for ssqrts in allsqrts:
            sqrts = ssqrts * TeV
            xsecs = [xsec.value.asNumber(fb)
                     for xsec in obj.weightList.getXsecsFor(sqrts)]
            if len(xsecs) != 1:
                logger.warning("Cross section lists contain multiple values for %s .\
                Only the first cross section will be printed" % str(sqrts))
            xsecs = float('%1.2e' %xsecs[0])
            sqrtsStr = 'xsec '+str(sqrts.asNumber(TeV))+' TeV'
            smsDict["Weights (fb)"][sqrtsStr] = xsecs

        return smsDict

    def _formatOutputStatus(self, obj):
        """
        Format data for a OutputStatus object.

        :param obj: A OutputStatus object to be printed.
        """

        infoDict = {}
        for key, val in obj.parameters.items():
            try:
                infoDict[key] = eval(val)
            except (NameError, TypeError, SyntaxError):
                infoDict[key] = val
        infoDict['file status'] = obj.filestatus
        infoDict['decomposition status'] = obj.status
        infoDict['warnings'] = obj.warnings
        infoDict['input file'] = obj.inputfile
        infoDict['database version'] = obj.databaseVersion
        infoDict['smodels version'] = obj.smodelsVersion
        # hidden feature, printtimespent, turn on in ini file, e.g.
        # [summary-printer] printtimespent = True
        if self.printtimespent:
            infoDict['time spent'] = "%.2fs" % (time.time() - self.time)
        return {'OutputStatus': infoDict}

    def _formatTheoryPredictionList(self, obj):
        """
        Format data of the TheoryPredictionList object.

        :param obj: A TheoryPredictionList object to be printed.
        """
        obj.sortTheoryPredictions()
        ExptRes = []
        for theoryPrediction in obj._theoryPredictions:
            expResult = theoryPrediction.expResult
            expID = expResult.globalInfo.id
            datasetID = theoryPrediction.dataId()
            dataType = theoryPrediction.dataType()
            ul = theoryPrediction.getUpperLimit()
            ulExpected = theoryPrediction.getUpperLimit(
                expected=self.getTypeOfExpected())
            if isinstance(ul, unum.Unum):
                ul = ul.asNumber(fb)
            if isinstance(ulExpected, unum.Unum):
                ulExpected = ulExpected.asNumber(fb)

            value = theoryPrediction.xsection.asNumber(fb)
            txnamesDict = {}
            for sms in theoryPrediction.smsList:
                if sms.txname.txName not in txnamesDict:
                    txnamesDict[sms.txname.txName] = sms.weight.asNumber(fb)
                else:
                    txnamesDict[sms.txname.txName] += sms.weight.asNumber(fb)
            maxconds = theoryPrediction.getmaxCondition()

            def _convWidth(x):
                if type(x) == type(GeV):
                    x = float(x.asNumber(GeV))
                if x == float("inf"):
                    x = "prompt"
                if x == 0.:
                    x = "stable"
                return x

            def roundme(x):
                if type(x) == tuple:
                    return (round(x[0].asNumber(GeV), 2), x[1].asNumber(GeV))
                return round(x.asNumber(GeV), 2)

            avgSMS = theoryPrediction.avgSMS
            if avgSMS is None:  # There is no commong topology
                mass = None
                widths = None
            elif self.outputFormat == "version2":
                mass = []
                widths = []
                for dIndex in avgSMS.daughterIndices(avgSMS.rootIndex):
                    daughter = avgSMS.indexToNode(dIndex)
                    mass.append([daughter.mass.asNumber(GeV)])
                    widths.append([_convWidth(daughter.totalwidth)])
                    for nodeIndex in avgSMS.dfsIndexIterator(dIndex):
                        node = avgSMS.indexToNode(nodeIndex)
                        if node.isSM:
                            continue
                        m = node.mass.asNumber(GeV)
                        mass[-1].append(m)
                        widths[-1].append(_convWidth(node.totalwidth))
            else:
                widthDict = {}
                massDict = {}
                for n in theoryPrediction.avgSMS.nodes:
                    if n.isSM or n is theoryPrediction.avgSMS.root:
                        continue
                    massDict[str(n)] = n.mass.asNumber(GeV)
                    widthDict[str(n)] = n.totalwidth.asNumber(GeV)
                mass = [(k,v) for k,v in massDict.items()]
                widths = [(k,v) for k,v in widthDict.items()]

            sqrts = expResult.globalInfo.sqrts

            r = self._round(theoryPrediction.getRValue(expected=False))
            r_expected = self._round(theoryPrediction.getRValue(
                expected=self.getTypeOfExpected()))

            resDict = {'maxcond': maxconds, 'theory prediction (fb)': self._round(value),
                       'upper limit (fb)': self._round(ul),
                       'expected upper limit (fb)': self._round(ulExpected),
                       'TxNames': sorted(txnamesDict.keys()),
                       'Mass (GeV)': mass,
                       'AnalysisID': expID,
                       'DataSetID': datasetID,
                       'AnalysisSqrts (TeV)': sqrts.asNumber(TeV),
                       'lumi (fb-1)': (expResult.globalInfo.lumi*fb).asNumber(),
                       'dataType': dataType,
                       'r': r, 'r_expected': r_expected,
                       'Width (GeV)' : widths}

            if hasattr(self, "addtxweights") and self.addtxweights:
                resDict['TxNames weights (fb)'] = txnamesDict
            llhd = theoryPrediction.likelihood()
            if llhd is not None:
                # resDict['chi2'] = self._round ( theoryPrediction.chi2 )
                resDict['likelihood'] = self._round(llhd)
                resDict['l_max'] = self._round(theoryPrediction.lmax())
                resDict['l_SM'] = self._round(theoryPrediction.lsm())
            ExptRes.append(resDict)

        return {'ExptRes': ExptRes}

    def _formatDoc(self, obj):
        """
        Format a pyslha object to be printed as a dictionary

        :param obj: pyslha object
        """

        MINPAR = dict(obj.blocks['MINPAR'].entries)
        EXTPAR = dict(obj.blocks['EXTPAR'].entries)
        mass = OrderedDict(obj.blocks['MASS'].entries.items())
        chimix = {}
        for key in obj.blocks['NMIX'].entries:
            val = obj.blocks['NMIX'].entries[key]
            if key[0] != 1:
                continue
            newkey = 'N'+str(key[0])+str(key[1])
            chimix[newkey] = val
        chamix = {}
        for key in obj.blocks['UMIX'].entries:
            val = obj.blocks['UMIX'].entries[key]
            newkey = 'U'+str(key[0])+str(key[1])
            chamix[newkey] = val
        for key in obj.blocks['VMIX'].entries:
            val = obj.blocks['VMIX'].entries[key]
            newkey = 'V'+str(key[0])+str(key[1])
            chamix[newkey] = val
        stopmix = {}
        for key in obj.blocks['STOPMIX'].entries:
            val = obj.blocks['STOPMIX'].entries[key]
            newkey = 'ST'+str(key[0])+str(key[1])
            stopmix[newkey] = val
        sbotmix = {}
        for key in obj.blocks['SBOTMIX'].entries:
            val = obj.blocks['SBOTMIX'].entries[key]
            newkey = 'SB'+str(key[0])+str(key[1])
            sbotmix[newkey] = val

        return {'MINPAR': MINPAR, 'chimix': chimix, 'stopmix': stopmix,
                'chamix': chamix, 'MM': {}, 'sbotmix': sbotmix,
                'EXTPAR': EXTPAR, 'mass': mass}

    def _formatUncovered(self, obj):
        """
        Format data of the Uncovered object containing coverage info

        :param obj: An Uncovered object to be printed.
        """

        # Number of missing topologies to be printed (ordered by cross sections)
        nprint = 10

        uncoveredDict = {}
        # First sort groups by label
        groups = sorted(obj.groups[:], key=lambda g: g.label)
        # Add summary of groups:
        for group in groups:
            sqrts = group.sqrts.asNumber(TeV)
            uncoveredDict["Total xsec for %s (fb)" % group.description] = \
                self._round(group.getTotalXSec())
            uncoveredDict["%s" % group.description] = []
            for fsSMS in group.finalStateSMS[:nprint]:
                fsSMSDict = {'sqrts (TeV)': sqrts, 'weight (fb)': self._round(fsSMS.missingX)}

                if self.outputFormat == 'version2':
                    smsStr = fsSMS.oldStr()
                    fsSMSDict['element'] = smsStr
                else:
                    smsStr = str(fsSMS)
                    fsSMSDict['SMS'] = smsStr


                if hasattr(self, "addsmslist") and self.addsmslist:
                    if self.outputFormat == "version2":
                        fsSMSDict["element IDs"] = [sms.smsID
                                                for sms in fsSMS._contributingSMS]
                    else:
                        fsSMSDict["SMS IDs"] = [sms.smsID
                                                for sms in fsSMS._contributingSMS]
                uncoveredDict["%s" % group.description].append(fsSMSDict)

        return uncoveredDict

    def _formatTheoryPredictionsCombiner(self, obj):
        """
        Format data of the TheoryPredictionsCombiner object.

        :param obj: A TheoryPredictionsCombiner object to be printed.
        """

        combRes = []  # In case we have a list of combined results in the future

        # Get list of analyses used in combination:
        expIDs = obj.analysisId()
        ul = obj.getUpperLimit()
        ulExpected = obj.getUpperLimit(expected=True)
        if isinstance(ul, unum.Unum):
            ul = ul.asNumber(fb)
        if isinstance(ulExpected, unum.Unum):
            ulExpected = ulExpected.asNumber(fb)

        r = self._round(obj.getRValue(expected=False))
        r_expected = self._round(obj.getRValue(expected=True))

        llhd = self._round(obj.likelihood())
        lmax = self._round(obj.lmax())
        lsm = self._round(obj.lsm())

        resDict = {'AnalysisID': expIDs,
                   'r': r, 'r_expected': r_expected,
                   'likelihood': llhd,
                   'l_max': lmax,
                   'l_SM': lsm}

        combRes.append(resDict)

        return {'CombinedRes': combRes}
