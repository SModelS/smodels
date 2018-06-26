#!/usr/bin/env python

"""
.. module:: testElementCompression
   :synopsis: Tests the compression functions in smodels.theory.element
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>   
"""
import sys
sys.path.insert(0,"../../")
import unittest
import smodels.theory.element as element
from smodels.SMparticleDefinitions import SMnames
from smodels.theory.particleClass import Particles
from smodels.tools.physicsUnits import GeV

class ElementCompressionTest(unittest.TestCase):
    
    def testInvisibleCompression(self):   
    
        el = element.Element( '[[[jet],[jet]],[[nue],[numu]]]' )     
        newel = el.invisibleCompress()
        
        self.assertEqual( str(newel), '[[],[[jet],[jet]]]' )        
            
     
    def testMassCompress(self):   
    
        p1 = Particles(Z2parity='odd', label='p1', pdg=None, mass=100.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
        p2 = Particles(Z2parity='odd', label='p2', pdg=1000021, mass=100.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
        p3 = Particles(Z2parity='odd', label='p3', pdg=1, mass=102.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
        p4 = Particles(Z2parity='odd', label='p4', pdg=None, mass=102.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
        p5 = Particles(Z2parity='odd', label='p5', pdg=1, mass=50.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
        p6 = Particles(Z2parity='odd', label='p6', pdg=None, mass=60.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)    
    
        el = element.Element( '[[[jet,jet]],[[jet,jet]]]' ) 
        el.branches[0].BSMparticles = [[p1, p3],[p2, p4]]
        el.branches[1].BSMparticles = [[p6, p5]]       # what happens if first mass is smaller?
        
        newel = el.massCompress(5.*GeV)
        
        self.assertEqual( str(newel), '[[],[[jet,jet]]]' )         

           
if __name__ == "__main__":
    unittest.main()   
          
