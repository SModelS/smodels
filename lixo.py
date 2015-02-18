#!/usr/bin/env python


#from smodels.tools.statistics import getUL,getPValue,computeCLInterval
import sys
from numpy import sqrt,inf
from scipy import stats,special,integrate,optimize
from smodels.tools.statistics import getUL


Nobs = 2
Nbg = 6.
NbgErr = 4.1

#print computeCLInterval(Nobs, Nbg, 1.)

#sys.exit()

x = getUL(Nobs,Nbg,NbgErr)
print x
#print getPValue(10.,Nobs,Nbg,NbgErr)



