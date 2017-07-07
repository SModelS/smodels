#!/usr/bin/python

from __future__ import print_function
from smodels.tools.statistics import UpperLimitComputer
from smodels.tools.physicsUnits import fb
import math
import time
import random


def run():
    f=open("bru.txt","w")
    computer = UpperLimitComputer ( 10000, 20.5 / fb, .95 )
    for i in range(200):
        # nobs_,nb_,deltab_,eff_=15,17.5,3.2,0.00454755
        nb_ = random.uniform ( 15., 25. )
        nobs_ = int ( random.gauss ( nb_, math.sqrt ( nb_ ) ) )
        deltab_ = random.gauss ( math.sqrt ( nb_ ), math.sqrt ( nb_/4. ) )
        if deltab_ < 0.: deltab_ = math.sqrt ( nb_ )
        eff_ = random.uniform ( 1e-5, 1e-1 )
        t0=time.time()
        ul = computer.ulSigmaTimesEpsilon ( nobs_, nb_, deltab_ )
        t1=time.time()-t0
        uls = computer.ulSigma ( [nobs_], [nb_], [[deltab_**2]], [eff_] )
        t2=time.time()-t1-t0
        line = "%d %f %f %f %f %f %f %f %f\n" % ( i, nb_, nobs_, deltab_, eff_, ul.asNumber(fb)/eff_, uls.asNumber(fb), t1, t2 )
        print ( line )
        f.write ( line )
    f.close()

run()
