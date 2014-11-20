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
   ignored. Z2-odd particles are decayed assuming Z2 convervation.

"""

particles = {
#R-odd particles:
#         PDG  label,2*spin,3*Qcharge,SU3charge,R-parity                
        1000021 : ["gluino",1,0,8,-1],
        1000022 : ["N1",1,0,1,-1],
        1000023 : ["N2",1,0,1,-1],
        1000025 : ["N3",1,0,1,-1],
        1000035 : ["N4",1,0,1,-1],
        1000024 : ["C1",1,3,1,-1],
        1000037 : ["C2",1,3,1,-1],
        1000039 : ["gravitino",3,0,1,-1],
        1000001 : ["squark",0,-1,3,-1],
        1000002 : ["squark",0,2,3,-1],
        1000003 : ["squark",-1],
        1000004 : ["squark",-1],
        2000001 : ["squark",-1],
        2000002 : ["squark",-1],
        2000003 : ["squark",-1],
        2000004 : ["squark",-1],
        1000005 : ["sbottom",-1],
        2000005 : ["sbottom",-1],
        1000006 : ["stop",-1],
        2000006 : ["stop",-1],
        1000011 : ["slepton",-1],
        1000013 : ["slepton",-1],
        1000015 : ["stau",-1],
        2000011 : ["slepton",-1],
        2000013 : ["slepton",-1],
        2000015 : ["stau",-1],
        1000012 : ["sneutrino",-1],
        1000014 : ["sneutrino",-1],
        1000016 : ["sneutrino",-1],
        2000012 : ["sneutrino",-1],
        2000014 : ["sneutrino",-1],
        2000016 : ["sneutrino",-1],
#R-even particles:
#         PDG  label,2*spin,3*Qcharge,SU3charge,R-parity   
        25 : ["higgs",1],
        -25 : ["higgs",1],
        35 : ["H0",1],
        -35 : ["H0",1],
        36 : ["A0",1],
        -36 : ["A0",1],
        37 : ["H+",1],
        -37 : ["H-",1],
        23 : ["Z",1],
        -23 : ["Z",1],
        22 : ["photon",1],
        -22 : ["photon",1],
        24 : ["W+",1],
        -24 : ["W-",1],
        16 : ["nu",1],
        -16 : ["nu",1],
        15 : ["ta-",1],
        -15 : ["ta+",1],
        14 : ["nu",1],
        -14 : ["nu",1],
        13 : ["mu-",1],
        -13 : ["mu+",1],
        12 : ["nu",1],
        -12 : ["nu",1],
        11 : ["e-",1],
        -11 : ["e+",1],
        5  : ["b",1],
        -5  : ["b",1],
        6  : ["t+",1],
        -6  : ["t-",1],
        1  : ["jet",1],
        2  : ["jet",1],
        3  : ["jet",1],
        4  : ["jet",1],
        21 : ["jet",1],
        -1  : ["jet",1],
        -2  : ["jet",1],
        -3  : ["jet",1],
        -4  : ["jet",1],
        -21 : ["jet",1]}

#Particle dictionary
ptcDic = {"e"  : ["e+", "e-"],
          "mu" : ["mu+", "mu-"],
          "ta" : ["ta+", "ta-"],
          "l+" : ["e+", "mu+"],
          "l-" : ["e-", "mu-"],
          "l"  : ["e-", "mu-", "e+", "mu+"],
          "W"  : ["W+", "W-"],
          "t"  : ["t+", "t-"],
          "L+" : ["e+", "mu+", "ta+"],
          "L-" : ["e-", "mu-", "ta-"],
          "L"  : ["e+", "mu+", "ta+", "e-", "mu-", "ta-"]}
