"""
.. module:: printer
   :synopsis: Facility used to print elements, theorypredictions, missing topologies et al
      in various forms

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>
.. moduleauthor:: Suchita Kulkanri <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import logging,sys,os
from smodels.theory.topology import TopologyList
from smodels.theory.element import Element
from smodels.theory.theoryPrediction import TheoryPredictionList
from smodels.experiment.expResultObj import ExpResult
from smodels.experiment.databaseObj import ExpResultList
from smodels.tools.ioObjects import OutputStatus, ResultList
from smodels.tools.coverage import UncoveredList, Uncovered
from smodels.tools.physicsUnits import GeV, fb, TeV
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from collections import OrderedDict
from xml.dom import minidom
from xml.etree import ElementTree

logger = logging.getLogger("smodels")


class MPrinter(object):
    """
    Master Printer class to handle the Printers (one printer/output type)
    
    :ivar printerList: list
    """
    def __init__(self):
        
        self.name = "master"
        self.Printers = {}
    
    def setPrinterOptions(self,parser):
        """
        Define the printer types and their options.
        
        :param parser: ConfigParser storing information from the parameters file
        """
        
        #Define the printer types and the printer-specific options:
        printerTypes = parser.get("printer", "outputType").split(",")        
        for prt in printerTypes:
            if prt == 'python':
                newPrinter = PyPrinter(output = 'file')                
            elif prt == 'summary':        
                newPrinter = SummaryPrinter(output = 'file')
            elif prt == 'stdout':
                newPrinter = TxTPrinter(output = 'stdout')
            elif prt == 'log':
                newPrinter = TxTPrinter(output = 'file')
            elif prt == 'xml':
                newPrinter = XmlPrinter(output = 'file')            
            else:
                logger.warning("Unknown printer format: %s" %str(prt))
                continue
            
            #Set printer-specific options:
            if parser.has_section(prt+'-printer'):
                newPrinter.setOptions(parser.items(prt+'-printer'))
            self.Printers[prt] = newPrinter
            
                


    def addObj(self,obj):
        """
        Adds the object to all its Printers:
        
        :param obj: An object which can be handled by the Printers.
        """
        
        for prt in self.Printers.values():
            prt.addObj(obj)
            
    def setOutPutFiles(self,filename):
        """
        Set the basename for the output files. Each printer will
        use this file name appended of the respective extension 
        (i.e. .py for a python printer, .smodels for a summary printer,...)
        
        :param filename: Input file name
        """
        
        for printer in self.Printers.values():
            printer.setOutPutFile(filename)


    def flush(self):
        """
        Ask all printers to write the output and clear their cache.
        If the printers return anything other than None, 
        we pass it on.
        """
        ret = {}
        for printerType,printer in self.Printers.items():
            ret[printerType] = printer.flush()
        return ret

class BasicPrinter(object):
    """
    Super class to handle the basic printing methods
    """

    def __init__(self, output, filename):
        self.name = "basic"

        self.outputList = []
        self.filename = filename
        self.output = output
        self.printingOrder = []
        self.toPrint = []


        if filename and os.path.isfile(filename):
            logger.warning("Removing file %s" %filename)
            os.remove(filename)
            
    def setOptions(self,options):
        """
        Store the printer specific options to control the output of each printer.
        Each option is stored as a printer attribute.
        
        :param options: a list of (option,value) for the printer.
        """
        
        for opt,value in options:
            setattr(self,opt,eval(value))        
            
    def addObj(self,obj):
        """
        Adds object to the Printer. 
        
        :param obj: A object to be printed. Must match one of the types defined in formatObj

        :return: True if the object has been added to the output. If the object does not belong
                to the pre-defined printing list toPrint, returns False.
        """
        
        for iobj,objType in enumerate(self.printingOrder):
            if isinstance(obj,objType):
                self.toPrint[iobj] = obj
                return True
        return False

    def openOutFile(self, filename, mode ):
        """ creates and opens a data sink, 
            creates path if needed """
        d = os.path.dirname ( filename )
        if not os.path.exists ( d ):
            os.makedirs ( d )
            logger.info ( "creating directory %s" % d )
        return open ( filename, mode )


    def flush(self):
        """
        Format the objects added to the output, print them to the screen
        or file and remove them from the printer.
        """
        ret=""

        for iobj,obj in enumerate(self.toPrint):
                if obj is None: continue
                output = self._formatObj(obj)                
                if not output: continue  #Skip empty output                
                ret += output
                if self.output == 'stdout':
                    sys.stdout.write(output)
                elif self.output == 'file':
                    if not self.filename:
                        logger.error('Filename not defined for printer')
                        return False   
                    with self.openOutFile(self.filename, "a") as outfile:
                        outfile.write(output)
                        outfile.close()

        self.toPrint = [None]*len(self.printingOrder)  #Reset printing objects
        return ret

    def _formatObj(self,obj):
        """
        Method for formatting the output depending on the type of object
        and output.
        
        :param obj: A object to be printed. Must match one of the types defined in formatObj

        """

        typeStr = type(obj).__name__
        try:
            formatFunction = getattr(self,'_format'+typeStr)
            return formatFunction(obj)
        except AttributeError,e:
            logger.debug('Error formating object %s: \n %s' %(typeStr,e))
            return False


class TxTPrinter(BasicPrinter):
    """
    Printer class to handle the printing of one single text output
    """
    def __init__(self, output = 'stdout', filename = None):
        BasicPrinter.__init__(self, output, filename)        
        self.name = "log"
        self.printingOrder = [OutputStatus,ExpResultList,TopologyList,Element,ExpResult,
                             TheoryPredictionList,ResultList,Uncovered]
        self.toPrint = [None]*len(self.printingOrder)        
        
    def setOutPutFile(self,filename,overwrite=True):
        """
        Set the basename for the text printer. The output filename will be
        filename.log.
        
        :param filename: Base filename
        :param overwrite: If True and the file already exists, it will be removed.
        """        
        
        self.filename = filename +'.' + self.name    
        if overwrite and os.path.isfile(self.filename):
            logger.warning("Removing old output file " + self.filename)
            os.remove(self.filename)
            
    def _formatDoc(self,obj):

        return False

    def _formatOutputStatus(self, obj):
        """
        Format data for a OutputStatus object.

        :param obj: A OutputStatus object to be printed.
        """

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
        output += "=" * 80 + "\n"
        return output

    def _formatTopologyList(self, obj):
        """
        Format data for a TopologyList object.

        :param obj: A TopologyList object to be printed.
        """

        if not hasattr(self,'printdecomp') or not self.printdecomp:
            return None

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
            if not hasattr(self,'addelementinfo') or not self.addelementinfo: continue
            for el in topo.elementList:
                output += "\t\t "+ 73 * "." + "\n"
                output += "\t\t Element: \n"
                output += self._formatElement(el) + "\n"

        return output


    def _formatElement(self, obj):
        """
        Format data for a Element object.

        :param obj: A Element object to be printed.
        """

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
        output += "\t\t The element weights are: \n \t\t " + obj.weight.niceStr().replace("\n", "\n \t\t ")

        return output

    def _formatExpResultList(self, obj):
        """
        Format data for a ExpResultList object.

        :param obj: A ExpResultList object to be printed.
        """
        
        if not hasattr(self,"printdatabase") or not self.printdatabase:
            return None

        output = ""

        for expRes in obj.expResultList:
            output += self._formatExpResult(expRes)

        return output+"\n"


    def _formatExpResult(self, obj):
        """
        Format data for a ExpResult object.

        :param obj: A ExpResult object to be printed.
        """

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
        if hasattr(self,"addanainfo") and self.addanainfo:
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
    
    
    def _formatResultList(self, obj):
        """
        Format data of the ResultList object.

        :param obj: A ResultList object to be printed.
        """
        
        obj.sort()
        
        
        if hasattr(self,"printextendedresults") and self.printextendedresults:
            return self._formatExpandedResultList(obj)
                
        if hasattr(self,"expandedSummary") and not self.expandedSummary:            
            theoPredictions = [obj.theoryPredictions[0]]
        else:
            theoPredictions = obj.theoryPredictions

        output = ""

        rvalues = []
        output += "#Analysis  Sqrts  Cond_Violation  Theory_Value(fb)  Exp_limit(fb)  r  r_expected"
        output += "\n\n"
        for theoPred in obj.theoryPredictions:
            expResult = theoPred.expResult
            datasetID = theoPred.dataset.dataInfo.dataId
            dataType = expResult.datasets[0].dataInfo.dataType
            txnames = theoPred.txnames
            if dataType == 'upperLimit':
                ul = expResult.getUpperLimitFor(txname=theoPred.txnames[0],mass=theoPred.mass)
                ul_expected = None
                signalRegion  = '(UL)'
            elif dataType == 'efficiencyMap':
                ul = expResult.getUpperLimitFor(dataID=datasetID)
                ul_expected = expResult.getUpperLimitFor(dataID=datasetID, expected=True)
                signalRegion  = theoPred.dataset.dataInfo.dataId
            else:
                logger.error("Unknown dataType %s" %(str(dataType)))
                raise SModelSError()
            value = theoPred.xsection.value

            r = (value/ul).asNumber()
            if type(ul_expected)==type(None): r_expected = None
            else: r_expected = (value/ul_expected).asNumber()
            rvalues.append(r)

            output += "%19s  " % (expResult.globalInfo.id)  # ana
            output += "%4s " % (expResult.globalInfo.sqrts/ TeV)  # sqrts
            output += "%5s " % theoPred.getmaxCondition()  # condition violation
            output += "%10.3E %10.3E " % (value.asNumber(fb), ul.asNumber(fb))  # theory cross section , expt upper limit
            if r_expected: output += "%10.3E %10.3E" % (r, r_expected)
            else: output += "%10.3E  N/A" %r
            output += "\n"
            output += " Signal Region:  "+signalRegion+"\n"
            txnameStr = str([str(tx) for tx in txnames])
            txnameStr = txnameStr.replace("'","").replace("[", "").replace("]","")
            output += " Txnames:  " + txnameStr + "\n"
            if hasattr(theoPred,'expectedUL') and not theoPred.expectedUL is None:
                output += " Chi2, Likelihood = %10.3E %10.3E\n" % (theoPred.chi2, theoPred.likelihood)            
            
            if not theoPred == obj.theoryPredictions[-1]: output += 80 * "-"+ "\n"

        output += "\n \n"
        output += 80 * "=" + "\n"
        output += "The highest r value is = " + str(max(rvalues)) + "\n"

        return output



    def _formatExpandedResultList(self, obj):
        """
        Format data for a ResultList object.

        :param obj: A ResultList object to be printed.
        """

        output = ""
        for theoryPrediction in obj.theoryPredictions:
            expRes = theoryPrediction.expResult
            info = theoryPrediction.dataset.dataInfo
            output += "\n"
            output += "---------------Analysis Label = " + expRes.globalInfo.id + "\n"
            output += "-------------------Dataset Label = " + str(info.dataId).replace("None","(UL)") + "\n"
            output += "-------------------Txname Labels = " + str([str(txname) for txname in theoryPrediction.txnames]) + "\n"
            output += "Analysis sqrts: " + str(expRes.globalInfo.sqrts) + \
                    "\n"


            if theoryPrediction.mass:
                for ibr, br in enumerate(theoryPrediction.mass):
                    output += "Masses in branch %i: " % ibr + str(br) + "\n"
            if hasattr(self,"addelementinfo") and self.addelementinfo:
                if not theoryPrediction.IDs[0]==0: output += "Contributing elements: " + str(theoryPrediction.IDs) + "\n"
            for pidList in theoryPrediction.PIDs:
                output += "PIDs:" + str(pidList) + "\n"
            output += "Theory prediction: " + str(theoryPrediction.xsection.value) + "\n"
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
                upperLimit = expRes.getUpperLimitFor(txname=theoryPrediction.txnames[0],mass=theoryPrediction.mass)
            elif expRes.datasets[0].dataInfo.dataType == 'efficiencyMap':
                upperLimit = expRes.getUpperLimitFor(dataID=theoryPrediction.dataset.dataInfo.dataId)

            output += "Observed experimental limit: " + str(upperLimit) + "\n"
            if expRes.datasets[0].dataInfo.dataType == 'efficiencyMap':
                output += "Expected experimental limit: " + str(expRes.getUpperLimitFor(dataID=theoryPrediction.dataset.dataInfo.dataId,expected=True)) + "\n"
            if hasattr(theoryPrediction,'expectedUL') and not theoryPrediction.expectedUL is None:
                output += "Chi2: " + str(theoryPrediction.chi2) + "\n"
                output += "Likelihood: " + str(theoryPrediction.likelihood) + "\n"

        return output



    def _formatUncovered(self, obj):
        """
        Format all uncovered data.
        
        :param obj: Uncovered object to be printed.
        """
                
        nprint = 10  # Number of missing topologies to be printed (ordered by cross-sections)

        output = ""
        output += "\nTotal missing topology cross section: %10.3E\n" %(obj.getMissingXsec())
        output += "Total cross section where we are outside the mass grid: %10.3E\n" %(obj.getOutOfGridXsec())
        output += "Total cross section in long cascade decays: %10.3E\n" %(obj.getLongCascadeXsec())
        output += "Total cross section in decays with asymmetric branches: %10.3E\n" %(obj.getAsymmetricXsec())        
        output += "\nFull information on unconstrained cross sections\n"
        output += "================================================================================\n"
        for ix, uncovEntry in enumerate([obj.missingTopos, obj.outsideGrid]):
            if len(uncovEntry.topos) == 0:
                if ix == 0:
                    output += "No missing topologies found\n"
                else:
                    output += "No contributions outside the mass grid\n"
            else:
                for topo in uncovEntry.topos:
                    if topo.value > 0.: continue
                    for el in topo.contributingElements:
                        if not el.weight.getXsecsFor(obj.missingTopos.sqrts): continue
                        topo.value += el.weight.getXsecsFor(obj.missingTopos.sqrts)[0].value.asNumber(fb)
                if ix==0:
                    output += "Missing topologies with the highest cross-sections (up to " + str(nprint) + "):\n"
                else:
                    output += "Contributions outside the mass grid (up to " + str(nprint) + "):\n"
                output += "Sqrts (TeV)   Weight (fb)        Element description\n"        
                for topo in sorted(uncovEntry.topos, key=lambda x: x.value, reverse=True)[:nprint]:
                    output += "%5s %10.3E    # %45s\n" % (str(obj.missingTopos.sqrts.asNumber(TeV)),topo.value, str(topo.topo))
                    if hasattr(self,"addelementinfo") and self.addelementinfo:
                        contributing = []
                        for el in topo.contributingElements:
                            contributing.append(el.elID)
                        output += "Contributing elements %s\n" % str(contributing)            
            output += "================================================================================\n"
        for ix, uncovEntry in enumerate([obj.longCascade, obj.asymmetricBranches]):
            if ix==0: output += "Long cascade decay by produced mothers (up to " + str(nprint) + "):\n"
            else: output += "Asymmetric branch decay by produced mothers\n"
            output += "Mother1 Mother2 Weight (fb)\n"
            for ent in uncovEntry.getSorted(obj.sqrts)[:nprint]:
                output += "%s %s %10.3E # %s\n" %(ent.motherPIDs[0], ent.motherPIDs[1], ent.getWeight(obj.sqrts).asNumber(fb), str(ent.motherPIDs))
                if hasattr(self,"addelementinfo") and self.addelementinfo:
                    contributing = []
                    for el in ent.contributingElements:
                        contributing.append(el.elID)
                    output += "Contributing elements %s\n" % str(contributing)
            if ix==0:
                output += "================================================================================\n"
        
        return output
                      

class SummaryPrinter(TxTPrinter):
    """
    Printer class to handle the printing of one single summary output.
    It uses the facilities of the TxTPrinter.
    """

    def __init__(self, output = 'stdout', filename = None):
        TxTPrinter.__init__(self, output, filename)
        self.name = "summary"
        self.printingOrder = [OutputStatus,ResultList,UncoveredList, Uncovered]
        self.toPrint = [None]*len(self.printingOrder)
        
    
    def setOutPutFile(self,filename,overwrite=True):
        """
        Set the basename for the text printer. The output filename will be
        filename.smodels.
        :param filename: Base filename
        :param overwrite: If True and the file already exists, it will be removed.
        """        
        
        self.filename = filename +'.smodels'
        if overwrite and os.path.isfile(self.filename):
            logger.warning("Removing old output file " + self.filename)
            os.remove(self.filename)
                     


class PyPrinter(BasicPrinter):
    """
    Printer class to handle the printing of one single pythonic output
    """
    def __init__(self, output = 'stdout', filename = None):
        BasicPrinter.__init__(self, output, filename)
        self.name = "py"
        self.printingOrder = [OutputStatus,TopologyList,Element,ResultList,Uncovered]
        self.toPrint = [None]*len(self.printingOrder)
        
    def setOutPutFile(self,filename,overwrite=True):
        """
        Set the basename for the text printer. The output filename will be
        filename.py.
        :param filename: Base filename
        :param overwrite: If True and the file already exists, it will be removed.
        """        
        
        self.filename = filename +'.py'
        if overwrite and os.path.isfile(self.filename):
            logger.warning("Removing old output file " + self.filename)
            os.remove(self.filename)

    def flush(self):
        """
        Write the python dictionaries generated by the object formatting
        to the defined output
        """
        
        outputDict = {}
        for iobj,obj in enumerate(self.toPrint):
            if obj is None: continue
            output = self._formatObj(obj)                
            if not output: continue  #Skip empty output
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
        ## it is a special feature of the python printer
        ## that we also return the output dictionary
        return outputDict

    def _formatTopologyList(self, obj):
        """
        Format data for a TopologyList object.

        :param obj: A TopologyList object to be printed.
        """
                
        if not hasattr(self,'addelementlist') or not self.addelementlist:
            return None

        elements = []

        for topo in obj:
            for el in topo.elementList:
                thisEl = self._formatElement(el)
                if thisEl: elements.append(thisEl)
                
        
        return {"Element": elements}

    def _formatElement(self, obj):
        """
        Format data for a Element object.

        :param obj: A Element object to be printed.
        """

        elDic = {}
        elDic["ID"] = obj.elID
        elDic["Particles"] = str(obj.getParticles())
        elDic["Masses (GeV)"] = [[m.asNumber(GeV) for m in br] for br in obj.getMasses()]
        elDic["PIDs"] = obj.getPIDs()
        elDic["Weights (fb)"] = {}
        sqrts = [info.sqrts for info in obj.weight.getInfo()]
        allsqrts = sorted(list(set(sqrts)))
        for sqrts in allsqrts:
            xsecs = [xsec.value.asNumber(fb) for xsec in obj.weight.getXsecsFor(sqrts)]
            if len(xsecs) != 1:
                logger.warning("Element cross-sections contain multiple values for %s .\
                Only the first cross-section will be printed" %str(sqrt))
            xsecs = xsecs[0]
            sqrtsStr = 'xsec '+str(sqrts.asNumber(TeV))+' TeV'
            elDic["Weights (fb)"][sqrtsStr] = xsecs
        return elDic

    def _formatOutputStatus(self, obj):
        """
        Format data for a OutputStatus object.

        :param obj: A OutputStatus object to be printed.
        """

        infoDict = {}
        for key,val in obj.parameters.items():
            try:
                infoDict[key] = eval(val)
            except (NameError,TypeError):
                infoDict[key] = val        
        infoDict['file status'] = obj.filestatus
        infoDict['decomposition status'] = obj.status
        infoDict['warnings'] = obj.warnings
        infoDict['input file'] = obj.inputfile
        infoDict['database version'] = obj.databaseVersion
        infoDict['smodels version'] = obj.smodelsVersion
        return {'OutputStatus' : infoDict}

    def _formatResultList(self, obj):
        """
        Format data of the ResultList object.

        :param obj: A ResultList object to be printed.
        """

        obj.sort()
        
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
                ulExpected = None
            elif dataType == 'efficiencyMap':
                ul = expResult.getUpperLimitFor(dataID=datasetID)
                ulExpected = expResult.getUpperLimitFor(dataID=datasetID, expected=True).asNumber(fb)
            else:
                logger.error("Unknown dataType %s" %(str(dataType)))
                continue            
            value = theoryPrediction.xsection.value
            txnames = [txname.txName for txname in theoryPrediction.txnames]
            maxconds = theoryPrediction.getmaxCondition()
            mass = theoryPrediction.mass
            if mass:
                mass = [[m.asNumber(GeV) for m in mbr] for mbr in mass]
            else:
                mass = None
            sqrts = dataset.globalInfo.sqrts
            resDict = {'maxcond': maxconds, 'theory prediction (fb)': value.asNumber(fb),
                        'upper limit (fb)': ul.asNumber(fb),
                        'expected upper limit (fb)': ulExpected,
                        'TxNames': txnames,
                        'Mass (GeV)': mass,
                        'AnalysisID': expID,
                        'DataSetID' : datasetID,
                        'AnalysisSqrts (TeV)': sqrts.asNumber(TeV),
                        'lumi (fb-1)' : (dataset.globalInfo.lumi*fb).asNumber(),
                        'dataType' : dataType}            
            if hasattr(theoryPrediction,'expectedUL') and not theoryPrediction.expectedUL is None:
                resDict['chi2'] = theoryPrediction.chi2
                resDict['likelihood'] = theoryPrediction.likelihood                
            ExptRes.append(resDict)

        ExptRes = sorted(ExptRes, key=lambda res: [res['theory prediction (fb)'],res['TxNames'],
                                                   res['AnalysisID'],res['DataSetID']])
        

        return {'ExptRes' : ExptRes}



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

    
    def _formatUncovered(self, obj):
        """
        Format data of the Uncovered object containing coverage info

        :param obj: A Uncovered object to be printed.
        """

        nprint = 10  # Number of missing topologies to be printed (ordered by cross-sections)

        missedTopos = []
        
        
        for topo in obj.missingTopos.topos:
            if topo.value > 0.: continue
            for el in topo.contributingElements:
                if not el.weight.getXsecsFor(obj.missingTopos.sqrts): continue
                topo.value += el.weight.getXsecsFor(obj.missingTopos.sqrts)[0].value.asNumber(fb)
        obj.missingTopos.topos = sorted(obj.missingTopos.topos, 
                                        key=lambda x: [x.value,str(x.topo)], 
                                        reverse=True)        
    
        for topo in obj.missingTopos.topos[:nprint]:
            missed = {'sqrts (TeV)' : obj.sqrts.asNumber(TeV), 'weight (fb)' : topo.value,
                                'element' : str(topo.topo)}           
            if hasattr(self,"addelementlist") and self.addelementlist:
                contributing = []
                for el in topo.contributingElements:
                    contributing.append(el.elID)
                missed["element IDs"] = contributing
            missedTopos.append(missed)
            
        outsideGrid = []
        for topo in obj.outsideGrid.topos:
            if topo.value > 0.: continue
            for el in topo.contributingElements:
                if not el.weight.getXsecsFor(obj.sqrts): continue
                topo.value += el.weight.getXsecsFor(obj.sqrts)[0].value.asNumber(fb)
        obj.outsideGrid.topos = sorted(obj.outsideGrid.topos, 
                                       key=lambda x: [x.value,str(x.topo)], 
                                       reverse=True)        
        for topo in obj.outsideGrid.topos[:nprint]:
            outside = {'sqrts (TeV)' : obj.sqrts.asNumber(TeV), 'weight (fb)' : topo.value,
                                'element' : str(topo.topo)}      
            outsideGrid.append(outside)     
        
        longCascades = []        
        obj.longCascade.classes = sorted(obj.longCascade.classes, 
                                         key=lambda x: [x.getWeight(obj.sqrts),sorted(x.motherPIDs[0:2])], 
                                         reverse=True)        
        for cascadeEntry in obj.longCascade.classes[:nprint]:
            longc = {'sqrts (TeV)' : obj.sqrts.asNumber(TeV),
                     'weight (fb)' : cascadeEntry.getWeight(obj.sqrts).asNumber(fb), 
                     'mother PIDs' : sorted(cascadeEntry.motherPIDs[0:2])}        
            longCascades.append(longc)
        
        asymmetricBranches = []
        obj.asymmetricBranches.classes = sorted(obj.asymmetricBranches.classes, 
                                                key=lambda x: [x.getWeight(obj.sqrts),sorted(x.motherPIDs[0:2])],
                                                reverse=True)
        for asymmetricEntry in obj.asymmetricBranches.classes[:nprint]:
            asymmetric = {'sqrts (TeV)' : obj.sqrts.asNumber(TeV), 
                    'weight (fb)' : asymmetricEntry.getWeight(obj.sqrts).asNumber(fb),
                    'mother PIDs' : sorted(asymmetricEntry.motherPIDs[0:2])}         
            asymmetricBranches.append(asymmetric)


        return {'Missed Topologies': missedTopos, 'Long Cascades' : longCascades,
                     'Asymmetric Branches': asymmetricBranches, 'Outside Grid': outsideGrid}
    


class XmlPrinter(PyPrinter):
    """
    Printer class to handle the printing of one single XML output
    """
    def __init__(self, output = 'stdout', filename = None):
        PyPrinter.__init__(self, output, filename)
        self.name = "xml"
        self.printingOrder = [OutputStatus,TopologyList,ResultList,Uncovered]
        self.toPrint = [None]*len(self.printingOrder)

        
    def setOutPutFile(self,filename,overwrite=True):
        """
        Set the basename for the text printer. The output filename will be
        filename.xml.
        :param filename: Base filename
        :param overwrite: If True and the file already exists, it will be removed.
        """        
        
        self.filename = filename +'.xml'
        if overwrite and os.path.isfile(self.filename):
            logger.warning("Removing old output file " + self.filename)
            os.remove(self.filename)        

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
            for key,val in sorted(pyObj.items()):
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
        for iobj,obj in enumerate(self.toPrint):
            if obj is None: continue
            output = self._formatObj(obj)  # Conver to python dictionaries                        
            if not output: continue  #Skip empty output            
            outputDict.update(output)

        root = None
        #Convert from python dictionaries to xml:
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
            ret = nice_xml

        self.toPrint = [None]*len(self.printingOrder)
        return root
