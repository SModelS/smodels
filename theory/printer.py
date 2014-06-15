"""
.. module:: printer
   :synopsis: Facility used in classes to derive from and be able to print
              different data types in different forms.

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
"""

from __future__ import print_function
import logging

logger = logging.getLogger(__name__)


class Printer(object):
    """
    Printer class.
    
    """
    def __init__(self):
        self.output = None


    def printout(self, target="stdout"):
        """
        Print the content of the data structure to the target.

        :param target: The target to print to. Possible values: stdout, html.
                       Default: stdout.
        :returns: None
        
        """
        # Branch, Element, Topology, TopologyList, Analysis, AnalysisList,
        # Cluster
        self.output = self.formatData()

        if target == "stdout":
            print(self.output)
        elif target == "file":
            pass


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
        Format data of a TopologyList object.
        
        """
        from tools import smsPrettyPrinter
        from prettytable import PrettyTable
        output = ""

        printer = smsPrettyPrinter.SmsPrettyPrinter()
        evTopTable = PrettyTable(["Topology", "#Vertices", "#Insertions",
                                  "#Elements", "Sum of weights"])
        evElementTable = PrettyTable(["Topology", "Element", "Particles B[0]",
                                      "Particles B[1]", "Masses B[0]",
                                      "Masses B[1]", "Element Weight"])

        eltot = 0
        # totweight = []
        # Print Results:
        # for i in range(len(SMSTopList)):
        for i, topo in enumerate(self):
            sumw = topo.getTotalWeight().getDictionary()
            evTopTable.add_row([i, topo.vertnumb, topo.vertparts,
                                len(topo.elementList),
                                smsPrettyPrinter.wrap(printer.pformat(sumw),
                                                      width=30)])
            eltot += len(topo.elementList)

        # Print element list for Topology[i]:
            if i == 0:
                for j, el in enumerate(topo.elementList):
                    if el.getParticles() != [[['b', 'b']], [['b', 'b']]]:
                        continue
                    m1 = printer.pformat(el.getMasses()[0])
                    m2 = printer.pformat(el.getMasses()[1])
                    dc = printer.pformat(el.weight.getDictionary())
                    row = [i, j, el.getParticles()[0], el.getParticles()[1],
                           smsPrettyPrinter.wrap(m1, width=25),
                           smsPrettyPrinter.wrap(m2, width=25),
                           smsPrettyPrinter.wrap(dc, width=30)]
                    evElementTable.add_row(row)
                evElementTable.add_row(["---", "---", "---", "---", "---",
                                        "---", "---"])

        output += "\n"
        output += "Number of Global topologies = " + str(len(self)) + "\n\n"
        output += str(evTopTable) + "\n\n"
        output += "Total Number of Elements = " + str(eltot) + "\n"
        output += "Total weight = " + str(self.getTotalWeight()) + "\n"
        # output += evElementTable + "\n"

        return output


    def formatTheoryPredictionData(self):
        """
        Format data of a TheoryPrediction object.
        
        """
        output = ""

        for theoryPrediction in self:
            output += "\n"
            output += "analysis: " + str(theoryPrediction.analysis) + "\n"
            output += "mass: " + str(theoryPrediction.mass) + "\n"
            output += "theory prediction: " + str(theoryPrediction.value) + \
                      "\n"
            output += "theory conditions:\n"
            if not theoryPrediction.conditions:
                output += "  " + theoryPrediction.conditions + "\n"
            else:
                for cond in theoryPrediction.conditions:
                    output += "  " + str(theoryPrediction.conditions[cond]) + \
                              "\n"
            experimentalLimit = theoryPrediction.analysis.getUpperLimitFor(
                    theoryPrediction.mass)
            output += "experimental limit: " + str(experimentalLimit) + "\n"

        return output
