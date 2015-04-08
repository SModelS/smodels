"""
.. module:: theory.printer
   :synopsis: Facility used in classes to derive from and be able to print
      different data types in different forms.

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>    
.. moduleauthor:: Suchita Kulkanri <suchita.kulkarni@gmail.com>

"""

from __future__ import print_function
import logging
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
    def __init__(self,outputs={'stdout' : None}, outputLevel = 1):
        """
        :param outputs: Dictionary defining the types of outputs to be printed and the respective
                        file names (for printing to the screen, set file name to None)
        :param outputLevel: general control for the output depth to be printed 
                           (0 = no output, 1 = basic output, 2 = detailed output,...)
        """               
        
        self.Printers = []
        for output,fname in outputs.items():
            printer = Printer()
            printer.outputType = output
            printer.filename = fname
            self.Printers.append(printer)
            
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
        

class Printer(object):
    """
    Printer class to handle the printing of one single output
    """
    def __init__(self):
                
        self.objList = []
        self.outputList = []
        self.outputType = 'stdout'
        self.outputLevel = 1
        self.filename = None


    def close(self):
        """
        Closes the printer and print the objects added to the output defined
        """
        
        if self.outputType == 'stdout':
            printingOrder = [OutputStatus,TopologyList,Element,
                             TheoryPredictionList,ResultList,MissingTopoList]
            for objType in printingOrder:
                for iobj,objB in enumerate(self.objList):
                    if objType == type(objB):
                        print(self.outputList[iobj])
                        
        elif self.outputType == 'summary':
            printingOrder = [OutputStatus,ResultList,MissingTopoList]
            if self.filename:                
                fout = open(self.filename,'w')
            for objType in printingOrder:
                for iobj,objB in enumerate(self.objList):
                    if objType == type(objB):
                        if self.filename:
                            fout.write(self.outputList[iobj])
                        else:
                            print(self.outputList[iobj])
            
            if self.filename: fout.close()               
            
        
    
    def addObj(self,obj):
        """
        Adds object to the Printer. The object will formatted according to the outputType
        and the outputLevel. The resulting output will be stored in outputList.
        :param obj: A object to be printed. Must match one of the types defined in formatObj
        """
        
        output = self.formatObj(obj)        
        self.objList.append(obj)  
        self.outputList.append(output) 
        
    def formatObj(self,obj):
        """
        Method for formatting the output depending on the type of object
        and output.
        :param obj: A object to be printed. Must match one of the types defined in formatObj
        """
        
        if isinstance(obj,TopologyList):
            return self.formatTopologyList(obj)
        elif isinstance(obj,Element):
            return self.formatElementList(obj)
        elif isinstance(obj,TxName):
            return self.formatTxName(obj)     
        elif isinstance(obj,TheoryPredictionList):
            return self.formatTheoryPredictionList(obj)
        elif isinstance(obj,OutputStatus):
            return self.formatOutputStatus(obj)   
        elif isinstance(obj,ResultList):
            return self.formatResultList(obj)
        elif isinstance(obj,MissingTopoList):
            return self.formatMissingTopoList(obj)


    def formatTopologyList(self, obj):
        """
        Format data for a TopologyList object. 
           
        :param obj: A TopologyList object to be printed.
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        if self.outputType != 'pyDict':
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
                        output += self.formatElement(el) + "\n"

        return output


    def formatElement(self, obj):
        """
        Format data for a Element object. 
           
        :param obj: A Element object to be printed.     
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        if self.outputType != 'pyDict':
            output = ""
            output += "\t\t Particles in element: " + str(obj.getParticles())
            output += "\n"
            output += "\t\t The element masses are \n"
            for i, el in enumerate(obj.getMasses()):
                output += "\t\t Branch %i: " % i + str(el) + "\n"
            output += "\t\t The element weights are: \n \t\t " + obj.weight.niceStr()

        return output

    def formatTxName(self, obj):
        """
        Format data for a TxName object. 
           
        :param obj: A TxName object to be printed.     
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        if self.outputType != 'pyDict':
            output = ""
            output += "========================================================\n"
            output += "Tx Label: "+obj.txname+'\n'
            if outputLevel == 2:
                output += "\t -----------------------------\n"
                output += "\t Elements tested by analysis:\n"            
                for el in obj._elements():
                    output += "\t    " + str(el.getParticles()) + "\n"
                
        return output


    def formatTheoryPredictionList(self, obj):
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
        if self.outputType != 'pyDict':
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

    def formatOutputStatus(self, obj):
        """
        Format data for a OutputStatus object. 
           
        :param obj: A OutputStatus object to be printed.     
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        if self.outputType != 'pyDict':
            output = ""
            output += "Input status: " + str(obj.filestatus) + "\n"
            output += "Decomposition output status: " + str(obj.status) + " " + obj.statusStrings[obj.status] + "\n"
            if obj.filestatus < 0: output += str(obj.warnings) + "\n"
            output += "#Input File: " + obj.inputfile + "\n"
            for label, par in obj.parameters.items(): output += "#" + label + " = " + str(par) + '\n'
            if obj.databaseVersion: output += "#Database version: %s\n" % obj.databaseVersion
            output += "================================================================================\n"

        return output

    def formatResultList(self, obj):
        """
        Format data of the ResultList object.
           
        :param obj: A ResultList object to be printed.     
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        output = ""

        if self.outputType != 'pyDict':
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

    def formatMissingTopoList(self, obj):
        """
        Format data of the MissingTopoList object.
           
        :param obj: A MissingTopoList object to be printed.     
        """

        outputLevel = self.outputLevel
        if not outputLevel: return None

        nprint = 10  # Number of missing topologies to be printed (ordered by cross-sections)

        if self.outputType != 'pyDict':
            output = ""
            output += "\n================================================================================\n"
            if len(obj.topos) == 0: return output + "No missing topologies found\n"
    
    
            output += "Missing topologies with the highest cross-sections (up to " + str(nprint) + "):\n"
            output += "Sqrts (TeV)   Weight (fb)        Element description\n"
    
            for topo in sorted(obj.topos, key=lambda x: x.value, reverse=True)[:nprint]:
                output += "%5s %10.3E    # %45s\n" % (str(obj.sqrts / TeV), topo.value, str(topo.topo))
                
        return output


def printout(obj, outputLevel=1):
    """
    Simple function for printing the object to the screen
    :param obj: object to be printed
    :param outputLevel: general control for the output depth to be printed 
                           (0 = no output, 1 = basic output, 2 = detailed output,...)
    """
        
    printer = Printer()    
    printer.outputLevel = outputLevel
    printer.outputType = 'stdout'
    printer.addObj(obj)
    printer.close()