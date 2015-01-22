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
from smodels.theory import crossSection
from smodels.experiment import smsResults
from smodels.tools.physicsUnits import GeV, fb, TeV

logger = logging.getLogger(__name__)


class Printer(object):
    """
    Printer class for defining specific print options and format for objects.    
    """
    def __init__(self):
        self.output = None


    def printout(self, target="stdout", filename="", outputLevel=1):
        """
        Print the content of the data structure to the target.

        :param target: The target to print to. Possible values: stdout, file.
           Default: stdout.
        :param filename: Filename to which the output is written
        :param outputLevel: general control for the output depth to be printed 
           (0 = no output, 1 = basic output, 2 = detailed output,...)
        :returns: None
        
        """
        # Branch, Element, Topology, TopologyList, Analysis, AnalysisList,
        # Cluster
        if outputLevel:
            self.output = self.formatData(outputLevel)

            if target == "stdout":
                print(self.output)
            elif target == "file":
                if not filename:
                    logger.error("Output filename not defined")
                    return
                f = open(filename, "a")
                f.write(self.output)
            elif target == "string":
                return self.output


    def formatData(self, outputLevel):
        """
        Format data of the derived object.
        
        Has to be implemented in the derived object. The real implementation is
        selected through dynamic binding.
        
        :raises: NotImplementedError  
              
        """
        raise NotImplementedError


    def formatTopologyListData(self, outputLevel):
        """
        Format data of to print Global topologies object.   
             
        :param outputLevel: general control for the output depth to be printed 
           (0 = no output, 1 = basic output, 2 = detailed output,...)
        
        """

        if not outputLevel: return None

        old_vertices = ""
        output = ""
        output += "   ======================================================= \n"
        output += " || \t \t\t\t\t\t\t || \n"
        output += " || \t \t Global topologies table \t \t ||\n"
        output += " || \t \t\t\t\t\t\t || \n"
        output += "   ======================================================= \n"
        for topo in self:
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
                    output += el.printout(target="string") + "\n"

        return output



    # This is my porposed format for element tabel

    def formatElementData(self, outputLevel):
        """
        Format data of to print an element object.
        
        :param outputLevel: general control for the output depth to be printed 
           (0 = no output, 1 = basic output, 2 = detailed output,...)
        """

        if not outputLevel: return None

        output = ""
        output += "\t\t Particles in element: " + str(self.getParticles())
        output += "\n"
        output += "\t\t The element masses are \n"
        for i, el in enumerate(self.getMasses()):
            output += "\t\t Branch %i: " % i + str(el) + "\n"
        output += "\t\t The element weights are: \n \t\t " + self.weight.niceStr()

        return output

    def formatULanalysisData(self, outputLevel):
        """
        Format data for a ULanalysis object.
        
        :param outputLevel: general control for the output depth to be printed 
           (0 = no output, 1 = basic output, 2 = detailed output,...)
        """

        if not outputLevel: return None

        output = ""
        output += "========================================================\n"
        output += "Analysis Name: " + self.label.split(":")[0] + '\n'
        output += "Tx Label: " + self.label.split(":")[1] + '\n'
        output += "Analysis Sqrts: " + str(self.sqrts) + '\n'
        if outputLevel == 2:
            output += "\t -----------------------------\n"
            output += "\t Elements tested by analysis:\n"
            listOfelements = []
            for el, eff in self.elementsEff.items():
                if eff > 0. and not el.getParticles() in listOfelements:
                    listOfelements.append(el.getParticles())
            for el in listOfelements:
                output += "\t    " + str(el) + "\n"

        return output




    def formatTheoryPredictionData(self, outputLevel):
        """
        Format data of a TheoryPrediction object.
        
        :param outputLevel: general control for the output depth to be printed 
           (0 = no output, 1 = basic output, 2 = detailed output,...)  
        """

        if not outputLevel: return None

        output = ""

        for theoryPrediction in self:
            output += "\n"
            output += "---------------Analysis Label = " + str(theoryPrediction.analysis.label) + "\n"
            output += "Analysis sqrts: " + str(theoryPrediction.value[0].info.label) + \
                    "\n"
            for i, el in enumerate(theoryPrediction.mass):
                output += "Masses in branch %i: " % i + str(el) + "\n"
            output += "Theory prediction: " + str(theoryPrediction.value[0].value) + \
                      "\n"
            output += "Theory conditions:"
            if not theoryPrediction.conditions:
                output += "  " + str(theoryPrediction.conditions) + "\n"
            else:
                condlist = []
                for cond in theoryPrediction.conditions:
                    condlist.append(theoryPrediction.conditions[cond])
                output += str(condlist) + "\n"
            experimentalLimit = theoryPrediction.analysis.getUpperLimitFor(
                    theoryPrediction.mass)
            output += "Experimental limit: " + str(experimentalLimit) + "\n"

        return output

    def formatStatusData(self, outputLevel):
        """
        Format data of the output status object.
        
        :param outputLevel: general control for the output depth to be printed 
           (0 = no output, 1 = basic output, 2 = detailed output,...)
        """

        if not outputLevel: return None

        output = ""
        output += "Input status: " + str(self.filestatus) + "\n"
        output += "Decomposition output status: " + str(self.status) + " " + self.statusStrings[self.status] + "\n"
        if self.filestatus < 0: output += str(self.warnings) + "\n"
        output += "#Input File: " + self.inputfile + "\n"
        for label, par in self.parameters.items(): output += "#" + label + " = " + str(par) + '\n'
        if self.databaseVersion: output += "#Database version: %s\n" % self.databaseVersion
        output += "================================================================================\n"

        return output

    def formatResultsData(self, outputLevel):
        """
        Format data of the final output object.
        
        :param outputLevel: general control for the output depth to be printed 
           (0 = no output, 1 = basic output, 2 = detailed output,...)
        """

        if not outputLevel: return None

        output = ""

        bestresult = self.getBestResult()

        if self.bestresultonly:
            output += "The result with highest R value is\n"
            self.outputarray = [bestresult]

        output += "#Analysis  Tx_Name  Sqrts  Cond. Violation  Theory_Value(fb)  Exp_limit(fb)  r\n\n"
        for op in self.outputarray:
            output += "%19s %16s " % (op.analysis.label.split(":")[0], op.analysis.label.split(":")[1])  # ana, topo
            output += "%4s " % (op.analysis.sqrts / TeV)  # sqrts
            output += "%5s " % op.getmaxCondition()  # condition violation
            output += "%10.3E %10.3E " % (op.value[0].value / fb, op.analysis.getUpperLimitFor(op.mass) / fb)  # theory cross section , expt upper limit
            output += "%10.3E\n" % self.getR(op)
            if self.describeTopo: output += "#" + str(smsResults.getConstraints(op.analysis.label.split(":")[0], op.analysis.label.split(":")[1])) + "\n"
            if not op == self.outputarray[-1]: output += "--------------------------------------------------------------------------------\n"

        output += "\n \n"
        output += "================================================================================\n"
        output += "The highest r value is = " + str(self.getR(bestresult)) + "\n"

        return output

    def formatMissingData(self, outputLevel):
        """
        Format data of missing topology list.
        
        :param outputLevel: general control for the output depth to be printed 
           (0 = no output, 1 = basic output, 2 = detailed output,...)
        """

        if not outputLevel: return None

        nprint = 10  # Number of missing topologies to be printed (ordered by cross-sections)

        output = ""
        output += "\n================================================================================\n"
        if len(self.topos) == 0: return output + "No missing topologies found\n"


        output += "Missing topologies with the highest cross-sections (up to " + str(nprint) + "):\n"
        output += "Sqrts (TeV)   Weight (fb)        Element description\n"

        for topo in sorted(self.topos, key=lambda x: x.value, reverse=True)[:nprint]:
            output += "%5s %10.3E    # %45s\n" % (str(self.sqrts / TeV), topo.value, str(topo.topo))
        return output
