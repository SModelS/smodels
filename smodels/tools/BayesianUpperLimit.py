from math import sqrt, exp, log
import math
import numpy as np

def upperLimit ( nev, sac, xbg, sbg, cl=.95, prec=None, smax=None ):
  """ Computes the upper limit.

    :param nev: number of observed events
    :param sac: relative uncertainty in acceptance
    :param xbg: expected background
    :param sbg: uncertainty in background
    :param  cl: desired CL
    :param prec: integration step size
    :param smax: upper limit of integration
  """
  if smax == None: 
    smax = 10. * ( nev + 1 )
  if prec == None:
    prec=xbg*.003
    if prec > 200:
        prec=xbg*0.00001
  if xbg < 0.:
    print "[BayesianUpperLimit] error: cannot deal with negative expected background"
    return float('nan')

  if sac < 0. or sac > 1.:
    print "[BayesianUpperLimit] error: signal acceptance must be between 0 and 1"
    return float('nan')

  if nev < 0:
    print "[BayesianUpperLimit] error: negative number of observed events."
    return float('nan')

  if sbg < 0.:
    print "[BayesianUpperLimit] error: negative sigma on bkg expectation!"
    return float('nan')

  if ( cl < 0. or cl > 1. ):
    print "[BayesianUpperLimit] error, confidence limit must be between 0 and 1"
    return float('nan')

  if ( smax < 2. * nev or smax > 100. * nev ):
    print "[BayesianUpperLimit] warning: strange choice for smax, smax=",\
          smax, ", nev=", nev
  if ( prec < 0. ):
    print "[BayesianUpperLimit] error: negative precision."
    return float('nan')

  if ( prec > .3 * nev and nev>0 ):
    print "[BayesianUpperLimit] error: precision too low, prec=",prec," nev=", nev
    return float('nan')

  bsum=0.
  xev=prec/2.
  blist=[0]*10000
  xlist=[0]*10000

  nlist =0

  while True:
    xlike = _blike ( nev,sac,xbg,sbg,xev)
    if ( math.isinf(xlike)) or math.isnan(xlike):
      if ( nlist==0 ):
        print "[BayesianUpperLimit] first likelihood is nan/inf! return nan!"
        return float('nan')

      print
      print "(D=" << nev << ", s=" << xev << ", l=" << xlike
    if nlist>=len(blist):
        print "[BayesianUpperLimit] need to extend array"
        for i in range(10000):
            blist.append(0)
            xlist.append(0)
    xlist[nlist]=xev
    blist[nlist]=xlike
    bsum+=xlike
    if False: 
      print "(D=" << nev << ", s=" << xev << ", l=" << xlike << ")\033[1A"
    xev+=prec
    if ( blist[nlist]/blist[0] < 1e-6 ): break
    if ( xev > smax ): break
    nlist+=1
  if False: print

  icl=0
  #  find place just below threshold
  bint=0.
  bcl=0.
  for i in range(nlist):
    if ( ( bint < cl * bsum ) and ( (bint+blist[i])>cl*bsum) ):
      icl=i
      bcl=bint
    bint+=blist[i]

  # interpolate linearly
  plim=xlist[icl]+( xlist[icl+1]-xlist[icl] ) * (cl*bsum-bcl)/blist[icl+1]
  return plim

def _blike ( nev, sac, xbg, sbg, xev ):
    """ return likelihood to observe nev events given expected background xbg, 
        error on background sbg, expected number of signal events xev, 
        error on signal acceptance sac """
    xint=0.

    ## perform double Gaussian integral by Monte Carlo

    nmax=2000
    for i in range(nmax):
        # pick expected background and signal from Gaussian
        bg, sig = -1., -1.
        while bg < 0. or sig < 0.:
            a,b =np.random.normal(), np.random.normal()
            bg = xbg + a * sbg
            sig = xev * ( 1. + b * sac )

        # total expected
        ex= bg + sig

        # value of integrand, poisson likelihood value
        # xxx = e (-ex) * ex^nev / nev! 
        xxx = exp (  nev * log ( ex ) - ex - math.lgamma(nev+1 ) )
        if math.isinf (xxx) or math.isnan ( xxx ):
            print "[blike] xxx=", xxx, " yyex=", yyex, " nev=", nev, ", nev!=", math.gamma(nev+1)
        xint +=xxx
    # print "[_blike] returns",xint/float(nmax)
    return xint/float(nmax)
