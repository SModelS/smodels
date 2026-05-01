#!/usr/bin/env python3

"""
.. module:: basicStats
   :synopsis: a module thats collects various statistical codes, that \
   are shared among other submodules.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
"""

from scipy import stats
from smodels.base.smodelsLogging import logger
import numpy as np
from smodels.statistics.exceptions import SModelSStatisticsError as SModelSError
from typing import Text, Union, Tuple, Optional
from collections.abc import Callable

__all__ = [ "NllEvalType", "CLsfromNLL", "determineBrentBracket",
            "exponentiateNLL" ]

from enum import Enum

def exponentiateNLL( nll : Optional[float],
        doIt : bool = True ) -> float:
    """ if doIt, then compute likelihood from nll,
    else return nll. Should become obsolete longterm
    """
    if nll == None:
        return None
        #if doIt:
        #    return 0.0
        #return 9000.0
    if doIt:
        return np.exp(-nll )
    return nll

class NllEvalType(Enum):
    """ an enum to account for the different types of likelihood values: observed,
    a priori evaluationType, a posteriori evaluationType """
    observed = 0
    aposteriori = 1
    apriori = 2

    @classmethod
    def init ( cls, evaluationType : Union[str,bool] ):
        """ get evaluationtype either from a string (e.g. 'posteriori') or a bool
            (true is priori, false is observed)
        """
        evaluationType = str(evaluationType).lower().replace("_","")
        if evaluationType in [ "posteriori", "aposteriori", "posterior" ]:
            return cls.aposteriori
        if evaluationType in [ "apriori", "prior", "priori", "true" ]:
            return cls.apriori
        if evaluationType in [ "false", "observed", "obs" ]:
            return cls.observed
        raise SModelSError ( f"NllEvalType {evaluationType} unknown" )

    def __hash__(self):
        return hash(self.value)

    def __eq__ ( self, other ):
        if type ( other ) == NllEvalType:
            return super().__eq__ ( other  )
        if type ( other ) in [ bool, str ]:
            return super().__eq__ ( NllEvalType.init ( other ) )
        raise SModelSError ( f"comparing a NllEvalType with {other}({type(other)})" )

## convenience
observed, aposteriori, apriori = NllEvalType.observed, NllEvalType.aposteriori, \
                                     NllEvalType.apriori

def clsType ( CLs : float, return_type : str, cl : float = 0.95 ) -> float:
    assert return_type in [ "1-CLs", "CLs", "CLs-alpha", "alpha-CLs" ], \
        f"return_type {return_type} unknown"
    if return_type == "1-CLs":
        return 1.0 - CLs
    if return_type == "CLs":
        return CLs
    if return_type == "CLs-alpha":
        return CLs - 1 + cl
    return 1 - cl - CLs

def CLsfromNLL(
        nllA: float, nll_minA: float, nll: float, nll_min: float,
        big_muhat : bool, return_type: Text = "CLs-alpha",
        nSigma : int = 0 ) -> float:
    """
    compute CLs (or similar) from the NLLs
    :param nllA: negative log likelihood for Asimov data
    :param nll_minA: negative log likelihood at muhat for Asimov data
    :param nll: negative log likelihood
    :param nll_min: negative log likelihood at muhat
    :param big_muhat: true if muhat>mu
    :param return_type: (Text) one of "CLs-alpha", "alpha-CLs", "1-CLs", or "CLs":
    ========== =======================
    CLs-alpha: returns CLs - 0.05
    alpha-CLs: returns 0.05 - CLs
    1-CLs:     returns 1-CLs value
    CLs:       returns "raw" CLs value
    ========== =======================
    :param nSigma: compute CLs not for central value but for this number of
    sigmas (positive or negative), see Equations 86 - 89 in the CCGV paper.

    :returns: Cls-type value, see above
    """
    assert return_type in ["CLs-alpha", "alpha-CLs", "1-CLs", "CLs"], \
           f"Unknown return type: {return_type}."
    qmu = 0.0 if ( nll < nll_min or big_muhat ) else 2 * (nll - nll_min)
    sqmu = np.sqrt(qmu)
    qA = 2 * (nllA - nll_minA)
    if qA < 0.0:
        qA = 0.0
    sqA = np.sqrt(qA)
    if qA >= qmu:
        CLsb = 1.0 - stats.norm.cdf(sqmu + nSigma )
        CLb = stats.norm.cdf(sqA - sqmu + nSigma )
    else:
        CLsb = 1. if qA == 0. else \
            1. - stats.norm.cdf((qmu + qA) / (2 * sqA) + nSigma )
        CLb = 1. if qA == 0. else \
            1. - stats.norm.cdf((qmu - qA) / (2 * sqA) + nSigma )

    CLs = CLsb / CLb if CLb > 0 else 0.0
    return clsType ( CLs, return_type, 0.95 )

def CLsWithErrorsfromNLL(
        nllA: float, nll_minA: float, nll: float, nll_min: float,
        s_nllA: float, s_nll_minA : float, s_nll: float, s_nll_min : float,
        big_muhat : bool, return_type: Text = "CLs-alpha",
        nSigma : int = 0 ) -> Tuple[float,float]:
    """
    compute CLs (or similar) from the NLLs
    :param nllA: negative log likelihood for Asimov data
    :param nll_minA: negative log likelihood at muhat for Asimov data
    :param nll: negative log likelihood
    :param nll_min: negative log likelihood at muhat
    :param s_nll*: the nlls' errors
    :param big_muhat: true if muhat>mu
    :param return_type: (Text) one of "CLs-alpha", "alpha-CLs", "1-CLs", or "CLs":
    ========== =======================
    CLs-alpha: returns CLs - 0.05
    alpha-CLs: returns 0.05 - CLs
    1-CLs:     returns 1-CLs value
    CLs:       returns "raw" CLs value
    ========== =======================
    :param nSigma: compute CLs not for central value but for this number of
    sigmas (positive or negative), see Equations 86 - 89 in the CCGV paper.

    :returns: Cls-type value, see above
    """
    assert return_type in ["CLs-alpha", "alpha-CLs", "1-CLs", "CLs"], \
           f"Unknown return type: {return_type}."
    qmu = 0.0 if ( nll < nll_min or big_muhat ) else 2 * (nll - nll_min)
    s_qmu = 0.0
    if nll >= nll_min and not big_muhat:
        var_qmu = 4 * ( s_nll**2 + s_nll_min**2 )

    sqrt_qmu = np.sqrt(qmu)
    var_sqrt_qmu = 0.
    if sqrt_qmu > 0.:
        var_sqrt_qmu = var_qmu / ( 4 * qmu )

    qA = 2 * (nllA - nll_minA)
    var_qA = 4 * ( s_nllA**2 + s_nll_minA**2 )

    if qA < 0.0:
        qA = 0.0
        var_qA = 0.

    sqrt_qA = np.sqrt(qA)
    var_sqrt_qA = 0.
    if sqrt_qA > 0.:
        var_sqrt_qA = var_qA / ( 4 * qA )

    if qA >= qmu:
        CLsb = 1.0 - stats.norm.cdf(sqrt_qmu + nSigma )
        var_CLsb = stats.norm.pdf(sqrt_qmu + nSigma )**2 * ( var_sqrt_qmu )
        x = sqrt_qA - sqrt_qmu + nSigma
        var_x = var_sqrt_qA + var_sqrt_qmu
        CLb = stats.norm.cdf( x )
        var_CLb = stats.norm.pdf( x )**2 * var_x
    else:
        s_CLsb, s_CLb = 0., 0.
        CLsb, CLb = 1., 1.
        var_CLsb, var_CLb = 0., 0. # FIXME is that a good choice?
        if qA != 0.:
            x = (qmu + qA) / (2 * sqrt_qA) + nSigma
            var_x = var_qmu  / ( 4* qA ) + var_sqrt_qA / 4
            CLsb = 1. - stats.norm.cdf( x )
            var_CLsb = stats.norm.pdf( x )**2 * var_x
            x2 = ( qmu - qA) / (2 * sqrt_qA ) + nSigma
            # s_x2 = s_sqrt_qmu / ( 2*sqrt_qA ) + s_sqrt_qA / 2.
            # (var_x2) is the same as var_x
            CLb = 1. - stats.norm.cdf( x2 )
            var_CLb = stats.norm.pdf ( x2 )**2 * var_x

    CLs = 0.
    s_CLs = 0.
    if CLb > 0:
        CLs = CLsb / CLb
        if abs(CLs)>1e-8: ##  CLs is close to zero, so this wont work
            s_CLs = float ( np.sqrt ( var_CLsb / CLb**2 + var_CLb * CLsb**2 / (CLb**4) ))
    # s_CLs is the same for them all
    ret = float ( clsType ( CLs, return_type, 0.95 ) )
    return ( ret, s_CLs )

def findRoot ( func : Callable, lower_bound : float, upper_bound : float,
        args : tuple = (), rtol : float = 8.881784197001252e-16,
        xtol : float = 2e-12 ) -> float:
    """ find the root of the function "func", within [lower_bound,upper_bound].
    We first try the faster toms748. If that doesnt converge,
    we fall back to brentq.
    :param func: the callable function to find the root for
    :param lower_bound: lower end of bracket
    :param upper_bound: upper end of bracket
    :param rtol: rtol for both toms748 and brentq
    :param xtol: xtol for both toms748 and brentq
    :returns: place where function evaluates to zero
    """
    logger.debug("Starting root finding")
    from scipy import optimize
    root = optimize.toms748( func, lower_bound, upper_bound, args=args,
                             rtol=rtol, xtol=xtol, full_output=True )
    if root[1].converged:
        root = root[0]
    else:
        root = optimize.brentq( func, lower_bound, upper_bound, args=args,
                                rtol=rtol, xtol=xtol )
    return root

def determineBrentBracket( mu_hat : float, sigma_mu : float ,
            rootfinder : Callable, allowNegative : bool = True,
            args : dict = {}, verbose :  bool = True ) -> Tuple:
    """find a, b for brent bracketing

    :param mu_hat: mu that maximizes likelihood
    :param sigm_mu: error on mu_hat (not too reliable)
    :param rootfinder: function that finds the root (usually root_func)
    :param allowNegative: if False, then do not allow a or b to become negative
    :returns: the interval a,b
    """
    msg = logger.error if verbose else logger.debug
    sigma_mu = max(sigma_mu, 0.5)  # there is a minimum on sigma_mu
    sigma_mu = min(sigma_mu, 100.) # there is a maximum on sigma_mu
    # the root should be roughly at mu_hat + 2*sigma_mu
    a = mu_hat + 1.5 * sigma_mu
    ra = rootfinder(a,**args)
    ntrials = 20
    i = 0
    foundExtra = False
    avalues = [ 0.0, 1.0, -1.0, 3.0, -3.0, 10.0, -10.0, 0.1, -0.1,
                0.01, -0.01, .001, -.001, 100., -100., .0001, -.0001 ]
    if not allowNegative:
        avalues = [0.0, 1.0, 3.0, 10.0, 0.1, 0.01, .001, 100., .0001 ]
    while ra < 0.0:
        # if this is negative, we move it to the left
        i += 1
        a -= (i**2.0) * sigma_mu
        if not allowNegative and a < 0:
            a = 0
        ra = rootfinder(a,**args)
        if i > ntrials or a < -10000.0 or ra is None:
            for a in avalues:
                ra = rootfinder(a,**args)
                if ra is None: # if cls computation failed, try with next a value
                    continue
                if ra > 0.0:
                    foundExtra = True
                    break
            if not foundExtra:
                msg(
                    f"cannot find an a that is left of the root. last attempt, a={a:.2f}, root = {ra:.2f}."
                )
                msg(f"mu_hat={mu_hat:.2f}, sigma_mu={sigma_mu:.2f}")
                raise SModelSError(
                    f"cannot find an a that is left of the root. last attempt, a={a:.2f}, root = {ra:.2f}."
                )
    i = 0
    foundExtra = False
    b = mu_hat + 2.5 * sigma_mu
    rb = rootfinder(b,**args)
    bvalues = [1.0, 0.0, 3.0, -1.0, 10.0, -3.0, 0.1, -0.1, -10.0,
               100.0, -100.0, 1000.0, .01, -.01, .001, -.001, 10000.0,
               100000.0, 1000000.0 ]
    if not allowNegative:
        bvalues = [ 1.0, 0.0, 3.0, 10.0, 0.1, 100.0, 1000.0, .01, .001,
                    10000.0, 100000.0, 1000000.0 ]
    while rb > 0.0:
        # if this is positive, we move it to the right
        i += 1
        b += (i**2.0) * sigma_mu
        rb = rootfinder(b,**args)
        if rb is None: # if cls computation failed, try with a bigger b value
            continue
        closestr, closest = float("inf"), None
        memoize = {}
        if i > ntrials or ( b < 0 and not allowNegative ):
            for b in bvalues:
                rb = rootfinder(b,**args)
                memoize[b]=float(rb)
                if rb is None: # if cls computation failed, try with next b value
                    continue
                if rb < 0.0:
                    foundExtra = True
                    break
                if rb < closestr:
                    closestr = rb
                    closest = b
            if not foundExtra:
                msg(f"cannot find a b that is right of the root (i.e. rootfinder(b) < 0).")
                msg(f"memoize {memoize}" )
                msg(f"closest to zero rootfinder({closest})={closestr}")
                msg(f"mu_hat was at {mu_hat:.2f} sigma_mu at {sigma_mu:.2f}")
                raise SModelSError( f"cannot find b right of the root" )
    return a, b

def deltaChi2FromLlhd(likelihood):
    """compute the delta chi2 value from a likelihood (convenience function)"""
    if likelihood == 0.0:
        return 1e10  # a very high number is good
    elif likelihood is None:
        return None

    return -2.0 * np.log(likelihood)

