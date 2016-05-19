"""
.. module:: printer
   :synopsis: Facility used to print elements, theorypredictions, missing topologies et al
      in various forms

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Suchita Kulkanri <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import logging,sys,os
from smodels.theory.topology import TopologyList
from smodels.theory.element import Element
from smodels.theory.theoryPrediction import TheoryPredictionList
from smodels.experiment.txnameObj import TxName
from smodels.experiment.expResultObj import ExpResult
from smodels.tools.ioObjects import OutputStatus, ResultList
from smodels.tools.coverage import UncoveredList, Uncovered
from smodels.tools.physicsUnits import GeV, fb, TeV
from pyslha import Doc
from smodels import installation
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from collections import OrderedDict
from xml.dom import minidom
from xml.etree import ElementTree

logger = logging.getLogger(__name__)


class MPrinter(object):
    """
    Master Printer class to handle the Printers (one printer/output type)
    """
    def __init__(self,*printerList):

        self.Printers = printerList

    def addObj(self,obj,objOutputLevel=None):
        """
        Adds the object to all its Printers:
        :param obj: An object which can be handled by the Printers.
        """

        for printer in self.Printers:
            printer.addObj(obj,objOutputLevel)

    def close(self):
        """
        Close all the Printers
        """

        for printer in self.Printers:
            printer.close()

    def flush(self):
        """
        Ask all printers to write the output and clear their cache.
        """

        for printer in self.Printers:
            printer.flush()

class TextBasedPrinter(object):
    """
    Super class to handle the printing of the text-based output
    """

    def __init__(self, output, filename, outputLevel):

        self.objList = []
        self.outputList = []
        self.outputLevel = outputLevel
        self.filename = filename
        self.output = output
        self.printingOrder = []


        if filename and os.path.isfile(filename):
            logger.warning("Removing file %s" %filename)
            os.remove(filename)

    def close(self):
        """
        Closes the printer and print the objects added to the output defined
        """

        self.flush()

    def flush(self):
        """
        Print the objects added to the output defined and removes them from the printer
        """

        for objType in self.printingOrder:
            for iobj,objB in enumerate(self.objList):
                if not self.outputList[iobj]: continue  #Skip empty output
                if objType == type(objB):
                    if self.output == 'stdout':
                        sys.stdout.write(self.outputList[iobj])
                    elif self.output == 'file':
                        if not self.filename:
                            logger.error('Filename not defined for printer')
                            return False
                        with open(self.filename, "a") as outfile:
                            outfile.write(self.outputList[iobj])
                            outfile.close()
        self.objList = []
        self.outputList = []

    def addObj(self,obj,objOutputLevel=None):
        """
        Adds object to the Printer. The object will be formatted according to the outputType
        and the outputLevel. The resulting output will be stored in outputList.
        :param obj: A object to be printed. Must match one of the types defined in formatObj
        :param outputLevel: Defines object specific output level. If set to None it will use
                            the printer outputLevel value.
        :return: True if the object has been added to the output. If the object does not belong
                to printingOrder or has no output format defined, returns False.
        """

        if objOutputLevel is None: objOutputLevel = self.outputLevel
        output = self._formatObj(obj,objOutputLevel)
        if output is False:
            return False
        self.objList.append(obj)
        self.outputList.append(output)
        return True


    def _formatObj(self,obj,objOutputLevel):
        """
        Method for formatting the output depending on the type of object
        and output.
        :param obj: A object to be printed. Must match one of the types defined in formatObj
        :param outputLevel: Defines object specific output level.
        """

        typeStr = type(obj).__name__
        try:
            formatFunction = getattr(self,'_format'+typeStr)
            return formatFunction(obj,objOutputLevel)
        except:
            return False


    def _formatDoc(self,obj,objOutputLevel):

        return False

    def _formatOutputStatus(self, obj, objOutputLevel):
        """
        Format data for a OutputStatus object.

        :param obj: A OutputStatus object to be printed.
        :param outputLevel: Defines object specific output level.
        """

        outputLevel = objOutputLevel
        if not outputLevel: return None

        output = ""
        output += "Input status: " + str(obj.filestatus) + "\n"
        output += "Decomposition output status: " + str(obj.status) + " "
        output += obj.statusStrings[obj.status] + "\n"
        if obj.filestatus < 0: output += str(obj.warnings) + "\n"
        output += "# Input File: " + obj.inputfile + "\n"
        labels = obj.parameters.keys()
        labels.sort()
        # for label, par in obj.parameters.items():
        for label in labels:
            par=obj.parameters[label]
            output += "# " + label + " = " + str(par) + '\n'
        if obj.databaseVersion:
            output += "# Database version: %s\n" % obj.databaseVersion
        output += "=" * 80 + "=\n"
        return output

    def _formatTopologyList(self, obj, objOutputLevel):
        """
        Format data for a TopologyList object.

        :param obj: A TopologyList object to be printed.
        :param outputLevel: Defines object specific output level.
        """

        if not objOutputLevel: return None

        old_vertices = ""
        output = ""
        output += "   ======================================================= \n"
        output += " || \t \t\t\t\t\t\t || \n"
        output += " || \t \t Global topologies table \t \t ||\n"
        output += " || \t \t\t\t\t\t\t || \n"
        output += "   ======================================================= \n"
        for topo in obj:
            if old_vertices == str(topo.vertnumb):
                output += "\t .................................................. \n"
            else:
                output += "===================================================== \n"
                output += "Topology:\n"
                output += "Number of vertices: " + str(topo.vertnumb) + ' \n'
                old_vertices = str(topo.vertnumb)
            output += "Number of vertex parts: " + str(topo.vertparts) + '\n'
            totxsec = topo.getTotalWeight()
            output += "Total Global topology weight :\n" + totxsec.niceStr() + '\n'
            output += "Total Number of Elements: " + str(len(topo.elementList)) + '\n'
            if objOutputLevel == 2:
                for el in topo.elementList:
                    output += "\t\t "+ 73 * "." + "\n"
                    output += "\t\t Element: \n"
                    output += self._formatElement(el,1) + "\n"

        return output


    def _formatElement(self, obj, objOutputLevel):
        """
        Format data for a Element object.

        :param obj: A Element object to be printed.
        :param outputLevel: Defines object specific output level.
        """

        if not objOutputLevel: return None

        output = ""
        output +="\t\t Element ID: " + str(obj.elID)
        output += "\n"
        output += "\t\t Particles in element: " + str(obj.getParticles())
        output += "\n"
        output += "\t\t The element masses are \n"
        for i, mass in enumerate(obj.getMasses()):
            output += "\t\t Branch %i: " % i + str(mass) + "\n"
        output += "\n"
        output += "\t\t The element PIDs are \n"
        for pidlist in obj.getPIDs():
            output += "\t\t PIDs: "+ str(pidlist) + "\n"
        output += "\t\t The element weights are: \n \t\t " + obj.weight.niceStr()

        return output

    def _formatTxName(self, obj, objOutputLevel):
        """
        Format data for a TxName object.

        :param obj: A TxName object to be printed.
        :param outputLevel: Defines object specific output level.
        """

        if not objOutputLevel: return None

        output = ""
        output += "========================================================\n"
        output += "Tx Label: "+obj.txName+'\n'
        if objOutputLevel == 2:
            output += "\t -----------------------------\n"
            output += "\t Elements tested by analysis:\n"
            for el in obj._elements():
                output += "\t    " + str(el.getParticles()) + "\n"

        return output

    def _formatExpResult(self, obj, objOutputLevel):
        """
        Format data for a ExpResult object.

        :param obj: A ExpResult object to be printed.
        :param outputLevel: Defines object specific output level.
        """

        if not objOutputLevel: return None

        txnames = []
        for dataset in obj.datasets:
            for txname in dataset.txnameList:
                tx = txname.txName
                if not tx in txnames: txnames.append(tx)


        output = ""
        output += "========================================================\n"
        output += "Experimental Result ID: " + obj.globalInfo.id + '\n'
        output += "Tx Labels: " + str(txnames) + '\n'
        output += "Sqrts: " + str(obj.globalInfo.sqrts) + '\n'
        if objOutputLevel == 2:
            output += "\t -----------------------------\n"
            output += "\t Elements tested by analysis:\n"
            listOfelements = []
            for dataset in obj.datasets:
                for txname in dataset.txnameList:
                    for el in txname._topologyList.getElements():
                        if not str(el) in listOfelements: listOfelements.append(str(el))
            for el in listOfelements:
                output += "\t    " + str(el) + "\n"

        return output



    def _formatTheoryPredictionList(self, obj, objOutputLevel):
        """
        Format data for a TheoryPredictionList object.

        :param obj: A TheoryPredictionList object to be printed.
        :param outputLevel: Defines object specific output level.
        """

        if not objOutputLevel: return None

        output = ""
        for theoryPrediction in obj:
            expRes = obj.expResult
            info = expRes.globalInfo
            datasetInfo = theoryPrediction.dataset.dataInfo
            output += "\n"
            output += "---------------Analysis Label = " + info.id + "\n"
            output += "-------------------Dataset Label = " + str(datasetInfo.dataId) + "\n"
            output += "-------------------Txname Labels = " + str([str(txname) for txname in theoryPrediction.txnames]) + "\n"
            output += "Analysis sqrts: " + str(info.sqrts) + \
                    "\n"

            if theoryPrediction.mass:
                for ibr, br in enumerate(theoryPrediction.mass):
                    output += "Masses in branch %i: " % ibr + str(br) + "\n"
            if not theoryPrediction.IDs[0]==0: output += "Contributing elements: " + str(theoryPrediction.IDs) + "\n"
            if objOutputLevel == 2:
                for pidList in theoryPrediction.PIDs:
                    output += "PIDs:" + str(pidList) + "\n"
            output += "Theory prediction: " + str(theoryPrediction.value) + "\n"
            output += "Theory conditions:"
            if not theoryPrediction.conditions:
                output += "  " + str(theoryPrediction.conditions) + "\n"
            else:
                condlist = []
                for cond in theoryPrediction.conditions:
                    condlist.append(theoryPrediction.conditions[cond])
                output += str(condlist) + "\n"

            #Get upper limit for the respective prediction:
            if expRes.datasets[0].dataInfo.dataType == 'upperLimit':
                experimentalLimit = expRes.getUpperLimitFor(txname=theoryPrediction.txnames[0],mass=theoryPrediction.mass)
            elif expRes.datasets[0].dataInfo.dataType == 'efficiencyMap':
                experimentalLimit = expRes.getUpperLimitFor(dataID=theoryPrediction.dataset.dataInfo.dataId)

            output += "Experimental limit: " + str(experimentalLimit) + "\n"

        return output

    def _formatResultList(self, obj, objOutputLevel):
        """
        Format data of the ResultList object.

        :param obj: A ResultList object to be printed.
        :param outputLevel: Defines object specific output level.
        """

        if not objOutputLevel: return None

        output = ""

        output += "#Analysis  Sqrts  Cond. Violation  Theory_Value(fb)  Exp_limit(fb)  r\n\n"
        # Suggestion to obtain expected limits:
        # output += "#Analysis  Tx_Name  Sqrts  Cond. Violation  PredTheory(fb)  ULobserved (fb) ULexpected (fb)  r\n\n"
        for theoPred in obj.theoryPredictions:
            expResult = theoPred.expResult
            datasetID = theoPred.dataset.dataInfo.dataId
            dataType = expResult.datasets[0].dataInfo.dataType
            txnames = theoPred.txnames
            if dataType == 'upperLimit':
                ul = expResult.getUpperLimitFor(txname=theoPred.txnames[0],mass=theoPred.mass)
                # Suggestion to obtain expected limits:
                # expected = expResult.getUpperLimitFor(txname=theoPred.txnames[0],mass=theoPred.mass, expected=True)
                signalRegion  = '(UL)'
            elif dataType == 'efficiencyMap':
                ul = expResult.getUpperLimitFor(dataID=datasetID)
                # Suggestion to obtain expected limits:
                # expected = expResult.getUpperLimitFor(dataID=datasetID, expected=True, compute=True)
                signalRegion  = theoPred.dataset.dataInfo.dataId
            else:
                logger.error("Unknown dataType %s" %(str(dataType)))
                raise SModelSError()

            output += "%19s  " % (expResult.globalInfo.id)  # ana
            output += "%4s " % (expResult.globalInfo.sqrts/ TeV)  # sqrts
            output += "%5s " % theoPred.getmaxCondition()  # condition violation
            output += "%10.3E %10.3E " % (theoPred.value[0].value / fb, ul / fb)  # theory cross section , expt upper limit
            # Suggestion to obtain expected limits:
            # output += "%10.3E " % (expected / fb)  # expected upper limit
            output += "%10.3E\n" % obj.getR(theoPred)
            output += " Signal Region:  "+signalRegion+"\n"
            if objOutputLevel == 2:
                txnameStr = str([str(tx) for tx in txnames])
                txnameStr = txnameStr.replace("'","").replace("[", "").replace("]","")
                output += " Txnames:  " + txnameStr + "\n"
            if not theoPred == obj.theoryPredictions[-1]: output += 80 * "-"+ "\n"

        output += "\n \n"
        output += 80 * "=" + "\n"
        output += "The highest r value is = " + str(obj.getR(obj.theoryPredictions[0])) + "\n"

        return output

    def _formatUncoveredList(self, obj, objOutputLevel):
        """
        Format data of the UncoveredList object of missing topology type.

        :param obj: A UncoveredList object to be printed.
        :param outputLevel: Defines object specific output level.
        """

        if not objOutputLevel: return None

        nprint = 10  # Number of missing topologies to be printed (ordered by cross-sections)

        output = ""
        output += "\n" + 80 * "=" + "\n"
        if len(obj.topos) == 0: return output + "No missing topologies found\n"

        for topo in obj.topos:
            if topo.value > 0.: continue
            for el in topo.contributingElements:
                if not el.weight.getXsecsFor(obj.sqrts): continue
                topo.value += el.weight.getXsecsFor(obj.sqrts)[0].value.asNumber(fb)

        output += "Missing topologies with the highest cross-sections (up to " + str(nprint) + "):\n"
        output += "Sqrts (TeV)   Weight (fb)        Element description\n"

        for topo in sorted(obj.topos, key=lambda x: x.value, reverse=True)[:nprint]:
            output += "%5s %10.3E    # %45s\n" % (str(obj.sqrts.asNumber(TeV)), topo.value, str(topo.topo))
            if objOutputLevel==2:
                contributing = []
                for el in topo.contributingElements:
                    contributing.append(el.elID)
                output += "Contributing elements %s\n" % str(contributing)

        return output

    def _formatUncovered(self, obj, objOutputLevel):
        """
        Format all uncovered data
        """

        if not objOutputLevel: return None

        output = ""
        if objOutputLevel == 1:
            output += "\nTotal missing topology cross section: %10.3E\n" %(obj.getMissingXsec())
            output += "Total cross section where we are outside the mass grid: %10.3E\n" %(obj.getOutOfGridXsec())
            output += "Total cross section in long cascade decays: %10.3E\n" %(obj.getLongCascadeXsec())
            output += "Total cross section in decays with asymmetric branches: %10.3E\n" %(obj.getAsymmetricXsec())
            return output
        if objOutputLevel ==2:
            output += "\nFull information on unconstrained cross sections\n"
            output += "================================================================================\n"
            output += "Contributions outside the mass grid\n"
            for topo in obj.outsideGrid.topos:
                for el in topo.contributingElements:
                    if not el.weight.getXsecsFor(obj.sqrts): continue
                    topo.value += el.weight.getXsecsFor(obj.sqrts)[0].value.asNumber(fb)
            for topo in sorted(obj.outsideGrid.topos, key=lambda x: x.value, reverse=True):
                output += "%5s %10.3E    # %45s\n" % (str(obj.sqrts.asNumber(TeV)), topo.value, str(topo.topo))
            output += "================================================================================\n"
            output += "Long cascade decay by produced mothers\n"
            output += "Mother1 Mother2 Weight (fb)\n"
            for cascadeEntry in obj.longCascade.classes:
                output += "%s %s %10.3E # %s\n" %(cascadeEntry.motherPIDs[0], cascadeEntry.motherPIDs[1], cascadeEntry.getWeight(obj.sqrts), str(cascadeEntry.motherPIDs))
            output += "================================================================================\n"
            output += "Asymmetric branch decay by produced mothers\n"
            output += "Mother1 Mother2 Weight (fb)\n"
            for asymmetricEntry in obj.asymmetricBranches.classes:
                output += "%s %s %10.3E # %s\n" %(asymmetricEntry.motherPIDs[0], asymmetricEntry.motherPIDs[1], asymmetricEntry.getWeight(obj.sqrts),asymmetricEntry.motherPIDs)
            return output


class TxTPrinter(TextBasedPrinter):
    """
    Printer class to handle the printing of one single text output
    """
    def __init__(self, output = 'stdout', filename = None, outputLevel = 1):

        TextBasedPrinter.__init__(self, output, filename, outputLevel)
        self.printingOrder = [OutputStatus,TopologyList,Element,ExpResult,
                             TheoryPredictionList,ResultList,UncoveredList, Uncovered]

class SummaryPrinter(TextBasedPrinter):
    """
    Printer class to handle the printing of one single summary output
    """

    def __init__(self, output = 'stdout', filename = None, outputLevel = 1):

        TextBasedPrinter.__init__(self, output, filename, outputLevel)
        self.printingOrder = [OutputStatus,ResultList,UncoveredList, Uncovered]


class PyPrinter(TextBasedPrinter):
    """
    Printer class to handle the printing of one single pythonic output
    """
    def __init__(self, output = 'stdout', filename = None, outputLevel = 1):

        TextBasedPrinter.__init__(self, output, filename, outputLevel)
        self.printingOrder = [OutputStatus,ResultList,UncoveredList]

    def flush(self):
        """
        Write the python dictionaries generated by the object formatting
        to the defined output
        """

        outputDict = {}
        for objType in self.printingOrder:
            for iobj,obj in enumerate(self.objList):
                if objType == type(obj):
                    objoutput = self.outputList[iobj]
                    outputDict.update(objoutput)

        if outputDict:
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

        self.objList = []
        self.outputList = []

    def _formatOutputStatus(self, obj, objOutputLevel):
        """
        Format data for a OutputStatus object.

        :param obj: A OutputStatus object to be printed.
        :param outputLevel: Defines object specific output level.
        """

        if not objOutputLevel: return None

        infoDict = dict(obj.parameters.items())
        infoDict['input file'] = obj.inputfile
        infoDict['database version'] = obj.databaseVersion
        infoDict['smodels version'] = obj.smodelsVersion
        return infoDict

    def _formatResultList(self, obj, objOutputLevel):
        """
        Format data of the ResultList object.

        :param obj: A ResultList object to be printed.
        :param outputLevel: Defines object specific output level.
        """

        if not objOutputLevel: return None

        ExptRes = []
        for theoryPrediction in obj.theoryPredictions:
            expResult = theoryPrediction.expResult
            dataset = theoryPrediction.dataset
            expID = expResult.globalInfo.id
            datasetID = dataset.dataInfo.dataId
            dataType = dataset.dataInfo.dataType
            if dataType == 'upperLimit':
                ul = expResult.getUpperLimitFor(txname=theoryPrediction.txnames[0],
                                                mass=theoryPrediction.mass)
            elif dataType == 'efficiencyMap':
                ul = expResult.getUpperLimitFor(dataID=datasetID)
            else:
                logger.error("Unknown dataType %s" %(str(dataType)))
            value = theoryPrediction.value[0].value
            txnames = [txname.txName for txname in theoryPrediction.txnames]
            maxconds = theoryPrediction.getmaxCondition()
            mass = theoryPrediction.mass
            if mass:
                daughterMass = mass[0][-1].asNumber(GeV)
                motherMass = mass[0][0].asNumber(GeV)
            else:
                daughterMass = None
                motherMass = None
            sqrts = dataset.globalInfo.sqrts
            ExptRes.append({'maxcond': maxconds, 'theory prediction (fb)': value.asNumber(fb),
                        'upper limit (fb)': ul.asNumber(fb),
                        'TxNames': txnames,
                        'DaughterMass (GeV)': daughterMass,
                        'MotherMass (GeV)': motherMass,
                        'AnalysisID': expID,
                        'DataSetID' : datasetID,
                        'AnalysisSqrts (TeV)': sqrts.asNumber(TeV),
                        'lumi (fb-1)' : (dataset.globalInfo.lumi*fb).asNumber(),
                        'dataType' : dataType})

        return {'ExptRes' : ExptRes}


    def _formatTheoryPredictionList(self, obj, objOutputLevel):
        """
        Format a TheoryPredictionList object to a python dictionary
        :param obj: TheoryPredictionList object
        :param outputLevel: Defines object specific output level.
        :return: python dictionary
        """

        return self._formatResultList(obj,objOutputLevel)


    def _formatDoc(self, obj, objOutputLevel):
        """
        Format a pyslha object to be printed as a dictionary

        :param obj: pyslha object
        :param outputLevel: Defines object specific output level.
        """

        if not objOutputLevel: return None

        MINPAR = dict(obj.blocks['MINPAR'].entries)
        EXTPAR = dict(obj.blocks['EXTPAR'].entries)
        mass = OrderedDict(obj.blocks['MASS'].entries.items())
        chimix = {}
        for key in obj.blocks['NMIX'].entries:
            val = obj.blocks['NMIX'].entries[key]
            if key[0] != 1: continue
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

        return {'MINPAR' : MINPAR, 'chimix' : chimix, 'stopmix' : stopmix,
                'chamix' : chamix, 'MM' : {}, 'sbotmix' : sbotmix,
                'EXTPAR' : EXTPAR, 'mass' : mass}

    def _formatUncoveredList(self, obj, objOutputLevel):
        """
        Format data of the UncoveredList object of missing topology type.

        :param obj: A UncoveredList object to be printed.
        :param outputLevel: Defines object specific output level.
        """

        if not objOutputLevel: return None

        nprint = 10  # Number of missing topologies to be printed (ordered by cross-sections)

        missedTopos = []
        for topo in sorted(obj.topos, key=lambda x: x.value, reverse=True)[:nprint]:
            missed = {'sqrts (TeV)' : obj.sqrts.asNumber(TeV), 'weight (fb)' : topo.value,
                                'element' : str(topo.topo)}

            if objOutputLevel==2:
                contributing = []
                for el in topo.contributingElements:
                    contributing.append(el.elID)
                missed["element IDs"] = contributing
            missedTopos.append(missed)
        return {'Missed Topologies': missedTopos}


class XmlPrinter(PyPrinter):
    """
    Printer class to handle the printing of one single XML output
    """
    def __init__(self, output = 'stdout', filename = None, outputLevel = 1):

        TextBasedPrinter.__init__(self, output, filename, outputLevel)
        self.printingOrder = [OutputStatus,ResultList,UncoveredList]

    def convertToElement(self,pyObj,parent,tag=""):
        """
        Convert a python object (list,dict,string,...)
        to a nested XML element tree.
        :param pyObj: python object (list,dict,string...)
        :param parent: XML Element parent
        :param tag: tag for the daughter element
        """

        tag = tag.replace(" ","_").replace("(","").replace(")","")
        if not isinstance(pyObj,list) and not isinstance(pyObj,dict):
            parent.text = str(pyObj).lstrip().rstrip()
        elif isinstance(pyObj,dict):
            for key,val in pyObj.items():
                key = key.replace(" ","_").replace("(","").replace(")","")
                newElement = ElementTree.Element(key)
                self.convertToElement(val,newElement,tag=key)
                parent.append(newElement)
        elif isinstance(pyObj,list):
            parent.tag += '_List'
            for val in pyObj:
                newElement = ElementTree.Element(tag)
                self.convertToElement(val,newElement,tag)
                parent.append(newElement)

    def flush(self):
        """
        Get the python dictionaries generated by the object formatting
        to the defined output and convert to XML
        """

        outputDict = {}
        for objType in self.printingOrder:
            for iobj,obj in enumerate(self.objList):
                if objType == type(obj):
                    objoutput = self.outputList[iobj]
                    outputDict.update(objoutput)

        if outputDict:
            root = ElementTree.Element('smodelsOutput')
            self.convertToElement(outputDict,root)
            rough_xml = ElementTree.tostring(root, 'utf-8')
            nice_xml = minidom.parseString(rough_xml).toprettyxml(indent="    ")
            if self.output == 'stdout':
                sys.stdout.write(nice_xml)
            elif self.output == 'file':
                if not self.filename:
                    logger.error('Filename not defined for printer')
                    return False
                with open(self.filename, "a") as outfile:
                    outfile.write(nice_xml)
                    outfile.close()

        self.objList = []
        self.outputList = []


def printout(obj, outputLevel=1):
    """
    Simple function for printing the object to the screen
    :param obj: object to be printed
    :param outputLevel: general control for the output depth to be printed
                           (0 = no output, 1 = basic output, 2 = detailed output,...)
    """

    printer = TxTPrinter()
    printer.outputLevel = outputLevel
    printer.output = 'stdout'
    printer.addObj(obj)
    printer.close()


