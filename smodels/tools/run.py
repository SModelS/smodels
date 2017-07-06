#!/usr/bin/python

from __future__ import print_function
from smodels.tools.statistics import UpperLimitComputer
from smodels.tools.physicsUnits import fb
import random

def run():
    computer = UpperLimitComputer ( 50000, 20.5 / fb, .95 )
    for i in range(100):
        nsig_,nobs_,nb_,deltab_,eff_=1,15,17.5,3.2,0.00454755
        # deltab_ = random.uniform ( 15., 25. )
        ul = computer.ulSigmaTimesEpsilon ( nobs_, nb_, deltab_ )
        uls = computer.ulSigma ( [nobs_], [nb_], [[deltab_**2]], [eff_] )
        print ( ul/eff_, uls )

run()
