#!/usr/bin/env python2.7
"""
.. module:: AuxiliaryFunctions
   :synopsis: Used in classes to derive from and be able to print different
   data types in different forms.

   .. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
"""

import logging
from theory import topology
from prettytable import PrettyTable
from tools import SMSPrettyPrinter
from tools.SMSPrettyPrinter import wrap
logger = logging.getLogger(__name__)


class Printer(object):
    """
    Printer class.
    
    """
    def __init__(self):
        """
        Constructor.
        
        """
        pass
    

    def printout(self, target="stdout"):
        """
        Prints the content of the data structure to the target.

        :param target: The target to print to. Possible values: stdout, html.
        Default: stdout.
        :returns: None
        
        """        
        #Branch, Element, Topology, TopologyList, Analysis, AnalysisList, Cluster
        self.output = ""
        
        if type(self) is topology.TopologyList:
            self.getTopologyListData()
        else:
            logger.error("Unknown type.")
            
        if target == "stdout":
            print self.output
        elif target == "file":
            pass
            
    
    def getTopologyListData(self):
        self.output += "printing TopologyList"        
            
        printer=SMSPrettyPrinter.SMSPrettyPrinter()
        EvTop_table = PrettyTable(["Topology","#Vertices", "#Insertions", "#Elements", "Sum of weights"])
        EvElement_table = PrettyTable(["Topology","Element","Particles B[0]","Particles B[1]", "Masses B[0]","Masses B[1]","Element Weight"])
        
        eltot = 0
        # totweight = []
        # Print Results:
        # for i in range(len(SMSTopList)):
        for i,topo in enumerate(self):
            sumw = topo.getTotalWeight().getDictionary()
            EvTop_table.add_row([i,topo.vertnumb,topo.vertparts,len(topo.ElList),wrap(printer.pformat(sumw),width=30)])
            eltot += len(topo.ElList)
        
         
              
        #Print element list for Topology[i]:  
            if i == 0:
                for j,el in enumerate(topo.ElList):
                    if el.getParticles() != [[['b','b']],[['b','b']]]:
                        continue
                    EvElement_table.add_row([i,j,el.getParticles()[0],el.getParticles()[1],wrap(printer.pformat(el.getMasses()[0]),width=25),wrap(printer.pformat(el.getMasses()[1]),width=25),wrap(printer.pformat(el.weight.getDictionary()),width=30)])
                EvElement_table.add_row(["---","---","---","---","---","---","---"])  
        
             
        print "Number of Global topologies = ",len(self)      
        print(EvTop_table)
        print "Total Number of Elements = ",eltot
        print "Total weight = ",self.getTotalWeight()
        # print(EvElement_table)
    
