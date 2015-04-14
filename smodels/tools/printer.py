"""
.. module:: printer
   :synopsis: Facility used to print elements, theorypredictions, missing topologies et al
      in various forms

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>    
.. moduleauthor:: Suchita Kulkanri <suchita.kulkarni@gmail.com>

"""

from __future__ import print_function
import logging,sys
from smodels.theory.topology import TopologyList
from smodels.theory.element import Element
from smodels.theory.theoryPrediction import TheoryPredictionList
from smodels.experiment.txnameObject import TxName
from smodels.tools.ioObjects import OutputStatus, ResultList
from smodels.tools.missingTopologies import MissingTopoList
from smodels.tools.physicsUnits import GeV, fb, TeV

logger = logging.getLogger(__name__)


class MPrinter(object):
    """
    Master Printer class to handle the Printers (one printer/output type)   
    """
    def __init__(self):
        
        self.Printers = []
            
    def addObj(self,obj):
        """
        Adds the object to all its Printers:
        :param obj: An object which can be handled by the Printers.
        """
        
        for printer in self.Printers:
            printer.addObj(obj)
            
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
    Super class to handle the printing of one text-based outputs
    """
    
    def __init__(self):
            
        self.objList = []
        self.outputList = []
        self.outputLevel = 1
        self.filename = None
        self.output = 'stdout'
        self.printingOrder = []
        
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
                if objType == type(objB):
                    if self.output == 'stdout': sys.stdout.write(self.outputList[iobj])
                    elif self.output == 'file':
                        if not self.filename:
                            logger.error('Filename not defiend for printer')
                            return False
                        with open(self.filename, "a") as outfile:
                            outfile.write(self.outputList[iobj])
                            outfile.close()
        self.objList = []
        self.outputList = []
        
    def addObj(self,obj):
        """
        Adds object to the Printer. The object will be formatted according to the outputType
        and the outputLevel. The resulting output will be stored in outputList.
        :param obj: A object to be printed. Must match one of the types defined in formatObj
        :return: True if the object has been added to the output. If the object does not belong
                to printingOrder or has no output format defined, returns False.        
        """
        
        output = self._formatObj(obj)
        if output is False:
            return False    
        self.objList.append(obj)  
        self.outputList.append(output)
        return True
    

    def _formatObj(self,obj):
        """
        Method for formatting the output depending on the type of object
        and output.
        :param obj: A object to be printed. Must match one of the types defined in formatObj
        """
        
        if isinstance(obj,TopologyList):
            return self._formatTopologyList(obj)
        elif isinstance(obj,Element):
            return self._formatElementList(obj)
        elif isinstance(obj,TxName):
            return self._formatTxName(obj)     
        elif isinstance(obj,TheoryPredictionList):
            return self._formatTheoryPredictionList(obj)
        elif isinstance(obj,OutputStatus):
            return self._formatOutputStatus(obj)   
        elif isinstance(obj,ResultList):
            return self._formatResultList(obj)
        elif isinstance(obj,MissingTopoList):
            return self._formatMissingTopoList(obj)
        else:
            return False     
        
    def _formatOutputStatus(self, obj):
        """
        Format data for a OutputStatus object. 
           
        :param obj: A OutputStatus object to be printed.     
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        output = ""
        output += "Input status: " + str(obj.filestatus) + "\n"
        output += "Decomposition output status: " + str(obj.status) + " " + obj.statusStrings[obj.status] + "\n"
        if obj.filestatus < 0: output += str(obj.warnings) + "\n"
        output += "#Input File: " + obj.inputfile + "\n"
        for label, par in obj.parameters.items(): output += "#" + label + " = " + str(par) + '\n'
        if obj.databaseVersion: output += "#Database version: %s\n" % obj.databaseVersion
        output += "================================================================================\n"

        return output

    def _formatTopologyList(self, obj):
        """
        Format data for a TopologyList object. 
           
        :param obj: A TopologyList object to be printed.
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

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
            if outputLevel == 2:
                for el in topo.elementList:
                    output += "\t\t .........................................................................\n"
                    output += "\t\t Element: \n"
                    output += self._formatElement(el) + "\n"

        return output


    def _formatElement(self, obj):
        """
        Format data for a Element object. 
           
        :param obj: A Element object to be printed.     
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        output = ""
        output += "\t\t Particles in element: " + str(obj.getParticles())
        output += "\n"
        output += "\t\t The element masses are \n"
        for i, el in enumerate(obj.getMasses()):
            output += "\t\t Branch %i: " % i + str(el) + "\n"
        output += "\t\t The element weights are: \n \t\t " + obj.weight.niceStr()

        return output

    def _formatTxName(self, obj):
        """
        Format data for a TxName object. 
           
        :param obj: A TxName object to be printed.     
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        output = ""
        output += "========================================================\n"
        output += "Tx Label: "+obj.txname+'\n'
        if outputLevel == 2:
            output += "\t -----------------------------\n"
            output += "\t Elements tested by analysis:\n"            
            for el in obj._elements():
                output += "\t    " + str(el.getParticles()) + "\n"
                
        return output


    def _formatTheoryPredictionList(self, obj):
        """
        Format data for a TheoryPredictionList object. 
           
        :param obj: A TheoryPredictionList object to be printed.     
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        output = ""

        expRes = obj.expResult
        info = expRes.info
        datasetInfo = obj.dataset.dataInfo        
        for theoryPrediction in obj:
            output += "\n"
            output += "---------------Analysis Label = " + info.id + "\n"
            output += "-------------------Dataset Label = " + str(datasetInfo.dataid) + "\n"
            output += "-------------------Txname Label = " + str(theoryPrediction.txname) + "\n"
            output += "Analysis sqrts: " + str(info.sqrts) + \
                    "\n"
            if theoryPrediction.mass:
                for ibr, br in enumerate(theoryPrediction.mass):
                    output += "Masses in branch %i: " % ibr + str(br) + "\n"
            output += "Theory prediction: " + str(theoryPrediction.value) + "\n"
            output += "Theory conditions:"
            if not theoryPrediction.conditions:
                output += "  " + str(theoryPrediction.conditions) + "\n"
            else:
                condlist = []
                for cond in theoryPrediction.conditions:
                    condlist.append(theoryPrediction.conditions[cond])
                output += str(condlist) + "\n"
            if datasetInfo.datatype == 'upper-limit':
                experimentalLimit = expRes.getUpperLimitFor(txname=theoryPrediction.txname,
                                                            mass=theoryPrediction.mass)
            elif datasetInfo.datatype == 'efficiency-map':
                experimentalLimit = expRes.getUpperLimitFor(dataID=datasetInfo.dataid)

            output += "Experimental limit: " + str(experimentalLimit) + "\n"

        return output

    def _formatResultList(self, obj):
        """
        Format data of the ResultList object.
           
        :param obj: A ResultList object to be printed.     
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        output = ""

        bestresult = obj.getBestResult()
        if obj.bestresultonly:
            output += "The result with highest R value is\n"
            obj.outputarray = [bestresult]
        output += "#Analysis  Tx_Name  Sqrts  Cond. Violation  Theory_Value(fb)  Exp_limit(fb)  r\n\n"
        for op in obj.outputarray:
            output += "%19s %16s " % (op.analysis.label.split(":")[0], op.analysis.label.split(":")[1])  # ana, topo
            output += "%4s " % (op.analysis.sqrts / TeV)  # sqrts
            output += "%5s " % op.getmaxCondition()  # condition violation
            output += "%10.3E %10.3E " % (op.value[0].value / fb, op.analysis.getUpperLimitFor(op.mass) / fb)  # theory cross section , expt upper limit
            output += "%10.3E\n" % obj.getR(op)
            if obj.describeTopo: output += "#" + str(op.analysis.constraint) + "\n"
            if not op == obj.outputarray[-1]: output += "--------------------------------------------------------------------------------\n"

        output += "\n \n"
        output += "================================================================================\n"
        output += "The highest r value is = " + str(obj.getR(bestresult)) + "\n"

        return output

    def _formatMissingTopoList(self, obj):
        """
        Format data of the MissingTopoList object.
           
        :param obj: A MissingTopoList object to be printed.     
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        nprint = 10  # Number of missing topologies to be printed (ordered by cross-sections)

        output = ""
        output += "\n================================================================================\n"
        if len(obj.topos) == 0: return output + "No missing topologies found\n"


        output += "Missing topologies with the highest cross-sections (up to " + str(nprint) + "):\n"
        output += "Sqrts (TeV)   Weight (fb)        Element description\n"

        for topo in sorted(obj.topos, key=lambda x: x.value, reverse=True)[:nprint]:
            output += "%5s %10.3E    # %45s\n" % (str(obj.sqrts / TeV), topo.value, str(topo.topo))
                
        return output 
  

class TxTPrinter(TextBasedPrinter):
    """
    Printer class to handle the printing of one single text output
    """
    def __init__(self):

        TextBasedPrinter.__init__(self)                
        self.printingOrder = [OutputStatus,TopologyList,Element,
                             TheoryPredictionList,ResultList,MissingTopoList]

class SummaryPrinter(TextBasedPrinter):
    """
    Printer class to handle the printing of one single summary output
    """
    
    def __init__(self):

        TextBasedPrinter.__init__(self)
        self.printingOrder = [OutputStatus,ResultList,MissingTopoList]
        

class PyPrinter(TextBasedPrinter):
    """
    Printer class to handle the printing of one single pythonic output
    """
    def __init__(self):

        TextBasedPrinter.__init__(self)                
        self.printingOrder = [OutputStatus,TopologyList,Element,
                             TheoryPredictionList,ResultList,MissingTopoList]
        


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


