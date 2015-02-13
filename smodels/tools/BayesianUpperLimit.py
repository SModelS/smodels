from math import sqrt, exp, log
import math
import numpy as np

def upperLimit ( nev, sac, xbg, sbg, cl=.95, prec=-1., smax=0. ):
  """ Computes the upper limit.

    :param nev: number of observed events
    :param sac: relative uncertainty in acceptance
    :param xbg: expected background
    :param sbg: uncertainty in background
    :param  cl: desired CL
    :param smax: upper limit of integration
    :param prec: integration step size
  """
  if smax < 1e-5 : 
    smax = 10. * nev
  if prec < 0.:
    prec=xbg*.001
  if xbg < 0.:
    print "[BayesianUpperLimit] error: cannot deal with negative expected background"
    return 0.

  if sac < 0. or sac > 1.:
    print "[BayesianUpperLimit] error: signal acceptance must be between 0 and 1"
    return 0.

  if nev < 0:
    print "[BayesianUpperLimit] error: negative number of observed events."
    return 0.

  if sbg < 0.:
    print "[BayesianUpperLimit] error: negative sigma on bkg expectation!"
    return 0.

  if ( cl < 0. or cl > 1. ):
    print "[BayesianUpperLimit] error, confidence limit must be between 0 and 1"

  if ( smax < 2. * nev or smax > 100. * nev ):
    print "[BayesianUpperLimit] warning: strange choice for smax, smax=",\
          smax, ", nev=", nev
  if ( prec < 0. ):
    print "[BayesianUpperLimit] error: negative precision."
    return 0.

  if ( prec > .3 * nev ):
    print "[BayesianUpperLimit] error: precision too low, prec=",prec," nev=", nev
    return 0.

  xevmax = smax
  dxev = prec
  bsum=0.
  xev=dxev/2.
  blist=[0]*10000
  xlist=[0]*10000

  nlist =0
  ## FIXME you are here
  while True:
    xlike = _blike ( nev,sac,xbg,sbg,xev)
    if ( math.isinf(xlike)) or math.isnan(xlike):
      if ( nlist==0 ):
        print "[BayesianUpperLimit] first likelihood is nan/inf! return -1!"
        return -1.
      print
      print "(D=" << nev << ", s=" << xev << ", l=" << xlike
      sys.exit(0)
    xlist[nlist]=xev
    blist[nlist]=xlike
    bsum+=xlike
    if False: 
      print "(D=" << nev << ", s=" << xev << ", l=" << xlike << ")\033[1A"
    xev+=dxev
    if ( blist[nlist]/blist[0] < 1e-6 ): break
    if ( xev > xevmax ): break
    nlist+=1
  if False: print

  icl=0
  bint=0.
  bcl=0.
  for i in range(nlist):
    if ( ( bint < cl * bsum ) and ( (bint+blist[i])>cl*bsum) ):
      icl=i
      bcl=bint
    bint+=blist[i]
  plim=xlist[icl]+( xlist[icl+1]-xlist[icl] ) * (cl*bsum-bcl)/blist[icl+1]
  return plim

def _factorial ( n ):
    if n<0:
        return 1.
    if n<9:
        return [1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880. ][n]
    return 2.506628*sqrt(n)*n**n*exp(-n)*(1.+1./12./n )


def _blike ( nev, sac, xbg, sbg, xev ):
    xxbg=xbg
    ssbg=sbg
    ssac = sac
    xxev=xev
    xint=0.

    nmax=2000
    for i in range(nmax):
        yybg=-1.
        yyev=-1.
        while ( (yybg < 0.) or (yyev < 0.) ):
            a,b =np.random.normal(), np.random.normal()
            yybg = xxbg + a * ssbg
            yyev = xxev * ( 1. + b * ssac )

        yyex=yybg + yyev

        xxx = exp (  nev * log ( yyex ) - yyex - math.lgamma(nev+1 ) )
        if math.isinf (xxx) or math.isnan ( xxx ):
            print "[blike] xxx=", xxx, " yyex=", yyex, " nev=", nev, ", nev!=", _factorial(nev)
        xint +=xxx
    return xint/float(nmax)
