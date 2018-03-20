#!/usr/bin/env python

"""
.. module:: particles
   :synopsis: Defines the list of R-even and R-odd particles to be used.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

   :parameter rOdd: dictionary with PDG codes for the rOdd (Z2-odd) particles
   and their respective labels
   :parameter rEven: dictionary with PDG codes for the rEven (Z2-eveb)
   particles and their respective labels
   :parameter ptcDic: dictionary with inclusive labels to help defining group
   of particles in the analysis database
   
   HOW TO ADD NEW PARTICLES: simply add a new entry in rOdd (rEven) if the
   particle is Z2-odd (Z2-even). For now all decays of Z2-even particles are
   ignored. Z2-odd particles are decayed assuming Z2 conservation.
   
   If you want to use slhaChecks to verify your input file (in the case of SLHA input
   only), also include the quantum numbers of the new particles in the qNumbers dictionary below.

"""

import os
from smodels.tools.runSModelS import modelFile

if "/" in modelFile:
    import shutil
    filename=os.path.basename(modelFile)
    shutil.copy ( modelFile, filename )
    modelFile=filename

if modelFile.endswith(".py"):
    modelFile=modelFile[:-3]

from importlib import import_module
pM=import_module (modelFile, package='smodels')
        
rOdd = pM.rOdd
rEven = pM.rEven
qNumbers = pM.qNumbers
