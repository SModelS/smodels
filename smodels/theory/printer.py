"""
.. module:: printer
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


    def printout(self, target="stdout", filename=""):
        """
        Print the content of the data structure to the target.

        :param target: The target to print to. Possible values: stdout, file.
                       Default: stdout.
        :param filename: Filename to which the output is written
        :returns: None
        
        """
        # Branch, Element, Topology, TopologyList, Analysis, AnalysisList,
        # Cluster
        self.output = self.formatData()

        if target == "stdout":
            print(self.output)
        elif target == "file":
            if not filename: return #FIXME need error message here!
            f = open(filename,"a")
            f.write(self.output)


    def formatData(self):
        """
        Format data of the derived object.
        
        Has to be implemented in the derived object. The real implementation is
        selected through dynamic binding.
        
        :raises: NotImplementedError  
              
        """
        raise NotImplementedError


    def formatTopologyListData(self):
        """
        Format data of to print Global topologies object.
        """
        import sys
        old_vertices = ""
        output = ""
        output += "   ======================================================= \n"
        output += " || \t \t\t\t\t\t\t || \n"
        output += " || \t \t Global topologies table \t \t ||\n"
        output += " || \t \t\t\t\t\t\t || \n"
        output += "   ======================================================= \n"
        for (i,topo) in enumerate(self):
            if old_vertices == str(topo.vertnumb):
                output += "\t .................................................. \n"
            else:
                output += "===================================================== \n"    
                output += "Number of vertices: " + str(topo.vertnumb) + ' \n'
                old_vertices = str(topo.vertnumb)
            output += "\t Number of vertex parts: " + str(topo.vertparts) + '\n'
            totxsec = crossSection.XSectionList()
            for el in topo.elementList:
                totxsec.combineWith(el.weight)
            output += "\t Total Global topology weight:\n" 
            for k in totxsec.getDictionary():
                pos=k.find(" " )
                sqrts=str( float(k[:pos]) / TeV) 
                output += "\t Sqrts: " + sqrts + "\t Weight: " + str(totxsec.getDictionary()[k]) + "\n"
        return output

    def formatElementData(self):
        """
            Format data of to print an element object.
        """
        output = ""
        output += "\t Particles in topology: " + str(self.getParticles())
        output += "\n"
        output += "\t The element masses are \n"
        for i, el in enumerate(self.getMasses()):
            output += "\t Branch %i: " %i+ str(el) + "\n"
        output += "\t The element weights are: \n"
        for k in self.weight.getDictionary():
            pos=k.find(" " )
            sqrts=str( float(k[:pos]) / TeV) 
            output += "\t Sqrts: " + sqrts + "\t Weight: " + str(self.weight.getDictionary()[k]) + "\n"

        return output
    
    def formatTheoryPredictionData(self):
        """
        Format data of a TheoryPrediction object.
        
        """
        output = ""

        for theoryPrediction in self:
            output += "\n"
            output += "---------------Analysis Label = " + str(theoryPrediction.analysis.label) + "\n"
            output += "Analysis sqrts: " + str(theoryPrediction.value[0].info.label) + \
                    "\n"
            for i, el in enumerate(theoryPrediction.mass):
                output += "Masses in branch %i: " %i+ str(el) + "\n"
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

    def formatStatusData(self):
        """
        Format data of the output status object.
        """
        output = ""

        output += "SLHA status: " + str(self.slhastatus) + "\n"
        output += "Decomposition output status: "+str(self.status)+" "+self.statusStrings[self.status] + "\n" #FIXME where should status strings go?
        if self.databaseVersion: output += "Database version: %s\n" % self.databaseVersion
        if self.slhastatus < 0: output += str(self.warnings) + "\n"
        output += "================================================================================\n"

        return output

    def formatResultsData(self):
        """
        Format data of the final output object.
        """
        output = ""

        bestresult = self.getBestResult()

        if self.bestresultonly:
            output += "The result with highest R value is\n"
            self.outputarray = [bestresult]

        output += "#Analysis  Topology  Sqrts  Cond_Violation  Theory_Value(fb)  Exp_limit(fb)  r\n\n"
        for op in self.outputarray:
            output += "%19s %16s " %(op.analysis.label.split(":")[0], op.analysis.label.split(":")[1]) # ana, topo
            output += "%4s " % op.analysis.sqrts / TeV # sqrts
            output += "%5s " % op.getmaxCondition() # condition violation
            output += "%10.3E %10.3E " % (op.value[0].value / fb,op.analysis.getUpperLimitFor(op.mass) / fb) # theory cross section , expt upper limit
            output += "%10.3E\n" % self.getR(op)
            if self.describeTopo: output += "#" + str(smsResults.getConstraints(op.analysis.label.split(":")[0],op.analysis.label.split(":")[1])) + "\n"
            if not op == self.outputarray[-1]: output += "--------------------------------------------------------------------------------\n"

        output += "\n \n"
        output += "================================================================================\n"
        output += "The highest R value is r_ratio = " + str(self.getR(bestresult)) + "\n"

        return output

    def formatSLHAData(self):
        """
        Format data of the slha checks output object.
        """
        output = ""

        output += "Input file: " + self.filename + "\n"
        output += "Sigmacut: " + str(self.sigmacut* fb) + "\n"
        output += "Minmassgap: " + str(self.massgap*GeV) + "\n"
        output += "Maxcond: " + str(self.maxcond) + "\n"
        if not self.status[0] == -3:
        #cannot add this information in case the input file is not slha format
            output += "LSP PID, mass: " + str(self.findLSP(returnmass=True)) + "\n"
            output += "NLSP PID, mass: " + str(self.findNLSP(returnmass=True)) + "\n"
        output += "================================================================================\n"

        return output

    def formatMissingData(self):
        """
        Format data of missing topology list.
        """

        output = ""

        output += "\n================================================================================\n"
        output += "Missing topologies with high cross sections:\n"
        output += "Sqrts (TeV)   Weight (fb)        Topology description\n"

        for topo in sorted(self.topos, key=lambda x: x.value, reverse=True)[:10]:
            output += "%5s %10.3E    # %45s\n" % (str(self.sqrts / TeV), topo.value, str(topo.topo))
        return output


    '''def formatElementList(self):
        output = ""
        for (i,topo) in enumerate(self):
            output += '\n'
            output += "====================================================================="
            output += "====================================================================="
            output += "A new global topoloy starts here" 
            output += "====================================================================="
            output += "====================================================================="
            for j, el in enumerate(topo.elementList):
                output += "........................................................................."
                output += "........................................................................."
                output += '\n'
                output += "Particles in topology:", el.getParticles()
                output += '\n'
                output += 'The element masses are'
                for item in range(len(el.getMasses())):
                    output += "Masses branch %i:" %item, el.getMasses()[item]
                output += "\n"
                output += "The element weights are:"
                for k in el.weight.getDictionary():
                    output += "Sqrts:", k, "\t Weights:", el.weight.getDictionary()[k]
        return output'''
