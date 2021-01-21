#!/usr/bin/env python3

"""
.. module:: idm
   :synopsis: Defines the BSM particles to be used. Properties not defined here and defined by the LHE or SLHA input file (such as masses, width and BRs) are automatically added later.

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>


"""

from smodels.theory.particle import Particle, MultiParticle

####  R-odd   ##########

H0 = Particle(Z2parity=-1, label='H0', pdg=35, eCharge=0, colordim=1, spin=0)
A0 = Particle(Z2parity=-1, label='A0', pdg=36, eCharge=0, colordim=1, spin=0)
H = Particle(Z2parity=-1, label='H+', pdg=37, eCharge=+1, colordim=1, spin=0)


rOdd = [H0, A0, H]
rOddC = [p.chargeConjugate() for p in rOdd]  #Define the charge conjugates

#Generic BSM particles:

BSMList = rOdd + rOddC
BSMparticleList = MultiParticle('BSM', BSMList)
