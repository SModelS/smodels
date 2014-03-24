#!/usr/bin/env python2.7
"""
.. module:: AuxiliaryFunctions
   :synopsis: Used in classes to derive from and be able to print different
   data types in different forms.

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
        # Branch, Element, Topology, TopologyList, Analysis, AnalysisList,
        # Cluster
        self.output = self.prepareData()        
            
        if target == "stdout":
            print(self.output)
        elif target == "file":
            pass
        
        
    def prepareData(self):
        raise NotImplementedError
