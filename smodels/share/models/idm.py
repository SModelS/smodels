#!/usr/bin/env python3

"""
.. module:: idm 
   :synopsis: Defines the BSM particles to be used.

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

   :parameter rOdd: list of particle objects for the rOdd (Z2-odd) particles
   :parameter rEven: list of particle objects for the rEven (Z2-even) particles

   Properties not defined here and defined by the LHE or SLHA input file 
   (such as masses, width and BRs) are automatically added later.
   
   HOW TO ADD NEW PARTICLES: simply add a new entry in rOdd (rEven) if the
   particle is Z2-odd (Z2-even). For now all decays of Z2-even particles are
   ignored. Z2-odd particles are decayed assuming Z2 conservation.
   
   If you want to use slhaChecks to verify your input file (in the case of SLHA input
   only), also include the quantum numbers of the new particles below.   
"""

from smodels.theory.particle import Particle, ParticleList

####  R-odd   ##########

H0 = Particle(Z2parity='odd', label='H0', pdg=35, eCharge=0, colordim=1, spin=0)  
A0 = Particle(Z2parity='odd', label='A0', pdg=36, eCharge=0, colordim=1, spin=0)  
H = Particle(Z2parity='odd', label='H+', pdg=37, eCharge=+1, colordim=1, spin=0)  


rOdd = [H0, A0, H]
rOddC = [p.chargeConjugate() for p in rOdd]  #Define the charge conjugates


from smodels.share.models.SMparticles import SMList
rEven = SMList

#Generic BSM particles:

BSMList = rOdd + rOddC
BSMparticleList = ParticleList('BSM', BSMList)
