#!/usr/bin/env python


#from smodels.tools.statistics import getUL,getPValue,computeCLInterval
import sys
from numpy import sqrt,inf
from scipy import stats,special,integrate,optimize
from smodels.tools.statistics import getUL,getPValue


Nobs = 24
Nbg = 3.3
NbgErr = 1.5

#print computeCLInterval(Nobs, Nbg, 1.)

#sys.exit()

x = getUL(Nobs,Nbg,NbgErr)
print x,getPValue(x,Nobs,Nbg,NbgErr)
#print getPValue(10.,Nobs,Nbg,NbgErr)



