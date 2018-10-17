#!/usr/bin/env python3

"""
.. module:: particles
   :synopsis: Defines the list of R-even and R-odd particles to be used.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

   :parameter rOdd: dictionary with PDG codes for the rOdd (Z2-odd) particles
   and their respective labels
   :parameter rEven: dictionary with PDG codes for the rEven (Z2-eveb)
   particles and their respective labels
   
   HOW TO ADD NEW PARTICLES: simply add a new entry in rOdd (rEven) if the
   particle is Z2-odd (Z2-even). For now all decays of Z2-even particles are
   ignored. Z2-odd particles are decayed assuming Z2 conservation.
   
   If you want to use slhaChecks to verify your input file (in the case of SLHA input
   only), also include the quantum numbers of the new particles in the qNumbers dictionary below.

"""


rOdd = {1000021 : "gluino",
        1000022 : "N1",
        1000023 : "N2",
        1000025 : "N3",
        1000035 : "N4",
        1000024 : "C1",
        1000037 : "C2",
        1000039 : "gravitino",
        1000001 : "squark",
        1000002 : "squark",
        1000003 : "squark",
        1000004 : "squark",
        2000001 : "squark",
        2000002 : "squark",
        2000003 : "squark",
        2000004 : "squark",
        1000005 : "sbottom",
        2000005 : "sbottom",
        1000006 : "stop",
        2000006 : "stop",
        1000011 : "slepton",
        1000013 : "slepton",
        1000015 : "stau",
        2000011 : "slepton",
        2000013 : "slepton",
        2000015 : "stau",
        1000012 : "sneutrino",
        1000014 : "sneutrino",
        1000016 : "sneutrino",
        2000012 : "sneutrino",
        2000014 : "sneutrino",
        2000016 : "sneutrino",
       -1000021 : "gluino",
       -1000022 : "N1",
       -1000023 : "N2",
       -1000025 : "N3",
       -1000035 : "N4",
       -1000024 : "C1",
       -1000037 : "C2",
       -1000039 : "gravitino",
       -1000001 : "squark",
       -1000002 : "squark",
       -1000003 : "squark",
       -1000004 : "squark",
       -2000001 : "squark",
       -2000002 : "squark",
       -2000003 : "squark",
       -2000004 : "squark",
       -1000005 : "sbottom",
       -2000005 : "sbottom",
       -1000006 : "stop",
       -2000006 : "stop",
       -1000011 : "slepton",
       -1000013 : "slepton",
       -1000015 : "stau",
       -2000011 : "slepton",
       -2000013 : "slepton",
       -2000015 : "stau",
       -1000012 : "sneutrino",
       -1000014 : "sneutrino",
       -1000016 : "sneutrino",
       -2000012 : "sneutrino",
       -2000014 : "sneutrino",
       -2000016 : "sneutrino",
        1000111 : "he+",
       -1000111 : "he-",
        1000255 : "s0",
       -1000255 : "s0"
}

rEven = {25 : "higgs",
        -25 : "higgs",
         35 : "H0",
        -35 : "H0",
         36 : "A0",
        -36 : "A0",
         37 : "H+",
        -37 : "H-",
         23 : "Z",
        -23 : "Z",
         22 : "photon",
        -22 : "photon",
         24 : "W+",
        -24 : "W-",
         18 : "nu",
        -18 : "nu",
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
 1000024:[1,3,1],
 1000037:[1,3,1],
 1000022:[1,0,1],
 1000023:[1,0,1],
 1000025:[1,0,1],
 1000035:[1,0,1],
 1000021:[1,0,8],
 1000011:[0,-3,1],
 2000011:[0,-3,1],
 1000013:[0,-3,1],
 2000013:[0,-3,1],
 1000015:[0,-3,1],
 2000015:[0,-3,1],
 1000012:[0,0,1],
 1000014:[0,0,1],
 1000016:[0,0,1],
 1000002:[0,2,3],
 2000002:[0,2,3],
 1000001:[0,-1,3],
 2000001:[0,-1,3],
 1000004:[0,2,3],
 2000004:[0,2,3],
 1000003:[0,-1,3],
 2000003:[0,-1,3],
 1000006:[0,2,3],
 2000006:[0,2,3],
 1000005:[0,-1,3],
 2000005:[0,-1,3],
 2000004:[0,2,3],
 1000111:[1,3,1],
 1000255:[0,0,1]
}


