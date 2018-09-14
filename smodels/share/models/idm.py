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


rOdd = {
         35 : "H0",
        -35 : "H0",
         36 : "A0",
        -36 : "A0",
         37 : "H+",
        -37 : "H-"
}

rEven = {25 : "higgs",
        -25 : "higgs",
         23 : "Z",
        -23 : "Z",
         22 : "photon",
        -22 : "photon",
         24 : "W+",
        -24 : "W-",
         16 : "nu",
        -16 : "nu",
         15 : "ta-",
        -15 : "ta+",
         14 : "nu",
        -14 : "nu",
         13 : "mu-",
        -13 : "mu+",
         12 : "nu",
        -12 : "nu",
         11 : "e-",
        -11 : "e+",
         4  : "c",
        -4  : "c",
         5  : "b",
        -5  : "b",
         6  : "t+",
        -6  : "t-",
         1  : "q",
         2  : "q",
         3  : "q",
         -1  : "q",
         -2  : "q",
         -3  : "q",
         21  : "g",
         -21  : "g",
         111: "pi",
         -111: "pi",
         211: "pi",
         -211: "pi" }

#Quantum numbers for the new particles.
#PDG: (spin*2, electrical charge*3, color dimension)
qNumbers={
 35:[0,0,1],
 36:[0,0,1],
 37:[0,3,1],
}
