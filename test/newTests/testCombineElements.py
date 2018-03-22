#!/usr/bin/env python

"""
.. module:: testElementCompression
   :synopsis: Tests the compression functions in smodels.theory.element
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>   
"""
import sys
sys.path.insert(0,"../../")
import unittest
import copy
import smodels.theory.element as element
from smodels.theory.topology import Topology
from smodels.theory.particleClass import Particles
from smodels.tools.physicsUnits import GeV

class CombineElementsTest(unittest.TestCase):
    
    def testCombineElements(self):   

        p1 = Particles(Z2parity='odd', label='p1', pdg=100, mass=100.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None) 
    
        el = element.Element( '[[[jet,jet]],[[jet,jet]]]' ) 
        el.branches[0].BSMparticles = [[p1]]
        el.branches[1].BSMparticles = [[p1]]       
        
        topo = Topology(el)
        
        newelement = copy.copy(el)
        p2 = Particles(Z2parity='odd', label='p2', pdg=101, mass=100.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None) 
        newelement.branches[0].BSMparticles = [[p2]]
        newelement.branches[1].BSMparticles = [[p2]]
        
        topo.addElement(newelement)
        

        self.assertEqual( len(topo.getElements()), 1)
        self.assertEqual( topo.getElements()[0], el)
        
        self.assertEqual( len(topo.getElements()[0].motherElements), 2)
        self.assertEqual( topo.getElements()[0].motherElements[0][1], el)
        self.assertEqual( topo.getElements()[0].motherElements[1][1], newelement)
           
if __name__ == "__main__":
    unittest.main()      
