"""
.. module:: AuxiliaryFunctions
   :synopsis: Facility used in classes to derive from and be able to print
   different data types in different forms.

   .. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
"""

from __future__ import print_function
import logging

logger = logging.getLogger(__name__) # pylint: disable-msg=C0103


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
        self.output = self.prepareData()
            
        if target == "stdout":
            print(self.output)
        elif target == "file":
            pass
        
        
    def prepareData(self):
        """
        Prepare data of the derived object.
        
        Has to be implemented in the derived object. The real implementation is
        selected through dynamic binding.
        
        :raises: NotImplementedError  
              
        """
        raise NotImplementedError
    
    
    def prepareTopologyListData(self):
        """
        Prepare data of a TopologyList object.
        
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
        
        #Print element list for Topology[i]:  
            if i == 0:
                for j, el in enumerate(topo.elementList):
                    if el.getParticles() != [[['b', 'b']], [['b', 'b']]]:
                        continue
                    evElementTable.add_row([
                            i, j, el.getParticles()[0], el.getParticles()[1],
                            smsPrettyPrinter.wrap(printer.pformat(
                                    el.getMasses()[0]), width=25),
                            smsPrettyPrinter.wrap(printer.pformat(
                                    el.getMasses()[1]), width=25),
                            smsPrettyPrinter.wrap(printer.pformat(
                                    el.weight.getDictionary()), width=30)])
                evElementTable.add_row(["---", "---", "---", "---", "---",
                                         "---", "---"])  
                     
        output += "Number of Global topologies = " + str(len(self)) + "\n\n"
        output += str(evTopTable) + "\n\n"
        output += "Total Number of Elements = " + str(eltot) + "\n"
        output += "Total weight = " + str(self.getTotalWeight()) + "\n"
        # output += evElementTable + "\n"
        
        return output
    
    
    def prepareTheoryPredictionData(self):
        """
        Prepare data of a TheoryPrediction object.
        
        """
        for theoryPrediction in self:
            print("mass:", theoryPrediction.mass)
            print("theory prediction:", theoryPrediction.value)
            print("theory conditions:")
            if not theoryPrediction.conditions:
                print(theoryPrediction.conditions)
            else:
                for cond in theoryPrediction.conditions:
                    print(theoryPrediction.conditions[cond])
            print("experimental limit:",
                  theoryPrediction.analysis.getUpperLimitFor(theoryPrediction.mass))
            print("\n")
