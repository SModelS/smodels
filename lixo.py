#!/usr/bin/env python


#from smodels.tools.statistics import getUL,getPValue,computeCLInterval
import sys
from numpy import sqrt,inf
from scipy import stats,special,integrate,optimize
from smodels.tools.statistics import getUL,getPValue


Nobs = 102
Nbg = 97.935
NbgErr = 0.005

#print computeCLInterval(Nobs, Nbg, 1.)

#sys.exit()

x = getUL(Nobs,Nbg,NbgErr)
print x,getPValue(x,Nobs,Nbg,NbgErr)
#print getPValue(10.,Nobs,Nbg,NbgErr)



