#!/usr/bin/env python

"""
.. module:: likelihoods
   :synopsis: Obtain likelihood distribution functions from known variables.
              Use observed and expected upper limits as input.
              Compute a chi2 with a theory cross section x Br prediction.
              Source: arXiv:1202.3415.

.. moduleauthor:: Jory Sonneveld <jory@opmijnfiets.nl>

"""

from scipy.stats import poisson
import scipy.stats as stats
import math
import logging
import scipy.optimize as opt
from unum import Unum
from smodels.tools.physicsUnits import pb, fb
import numpy as np

#########################################################
### Upper limits from the mean and standard deviation ###
###               The inverse function                ###
#########################################################
def obs_exp(mean, std, start=0.0):
    """(float, float) -> float, float
    Computes the observed and expected upper limits from a mean
    and standard devation (std).
    This is the inverse of the mean and standard deviation
    computed from upper limits following the method described in
    1202.3415.
    Returns the observed and expected upper limits.

    Use that the confidence level for the limit is equal to
    _cl_from_erfs(mean, std, limit) = 0.95
    and that the mean is the
    mu = upper limit / theory prediction
    and
    mu_exp = 1.96 * sigma_exp
    and
    sigma_exp=sigma_obs. See arXiv:1202.3415.

    Examples of fsolve (note the negative starting point):
    >>> def poly(x, y): return x**2 + y**3 - 1
    >>> opt.fsolve(poly, 2, args=(-2))
    array([ 3.])
    >>> opt.fsolve(poly, 2, args=(-4))
    array([ 8.06225775])
    >>> opt.fsolve(poly, 0, args=(-2))
    array([ 3.])
    >>> opt.fsolve(poly, 100, args=(-2))
    array([ 3.])
    >>> opt.fsolve(poly, -100, args=(-2))
    array([-3.])

    Examples of mean_std and obs_exp:
    >>> mean, std = mean_std(11)
    >>> obs, exp = obs_exp(mean, std)
    >>> (round(obs), round(exp))
    (11.0, 11.0)
    >>> mean, std = mean_std(11, 10, lumi=20)
    >>> obs, exp = obs_exp(mean, std)
    >>> (round(obs), round(exp))
    (11.0, 10.0)
    >>> mean, std = mean_std(11, lumi=20)
    >>> obs, exp = obs_exp(mean, std)
    >>> (round(obs), round(exp))
    (11.0, 11.0)
    >>> mean, std = mean_std(11, 10, lumi=20)
    >>> obs, exp = obs_exp(mean, std)
    >>> (round(obs), round(exp))
    (11.0, 10.0)
    """

    # Don't do anything if the mean or standard devation are None:
    if mean == None or std == None:
        return None, None

    # Compute the observed upper limit by requiring 95% C.L.:
    obs =  opt.fsolve(_cl_observed, start, args=(mean, std, 0.95))[0]

    # Remember that the expected upper limit is equal to 1.96 * std
    exp = 1.96 * std

    return obs, exp


#########################################################
###############    Helper functions     #################
#########################################################
def _poisson_cls(s, b, d, pvalue_desired=0.05):
    '''(float, float, int) -> float
    Calculated the Poisson expected p-value according
    to the CLs technique for one channel.
    Return this minus the desired confidence level

    >>> round(_poisson_cls(math.log(1/0.05), 1, 0, 0), 2)
    0.05
    >>> val = round(math.exp(-1) * 1.5, 2)
    >>> valp = round(_poisson_cls(1, 1, 1, 0), 2)
    >>> val == valp
    True

    '''
    CLs_plus_b = poisson.cdf(d, s+b)
    CLb = poisson.cdf(d, b)
    return CLs_plus_b/float(CLb) - pvalue_desired

def _simple_signal_calculator(background, data):
    '''(float, int) -> float
    Take the number of background events (can be
    a decimal), the number of observed data events
    (must be an integer)
    and calculate the signal given that it was
    Poisson distributed and yields a 0.95 %
    confidence level limit.
    '''
    # solve poisson CLs for 95 % CL, start at s=b:
    signal = opt.fsolve(_poisson_cls, data-background,\
            args=(background, data, 0.05))
    sig = signal[0]
    return sig

def _in_picobarn(quantity):
    """(float/Unum/int) -> float
    Return the given quantity
        1) in picobarn if it is a Unum instance
        2) as is if it is not a Unum instance
    """
    if isinstance(quantity, Unum):
        return quantity.asNumber(pb)
    return quantity

def _in_inverse_picobarn(quantity):
    """(float/Unum/int) -> float
    Return the given quantity
        1) in inverse picobarn if it is a Unum instance
        2) as is if it is not a Unum instance
    """
    if isinstance(quantity, Unum):
        return quantity.asNumber(pb**(-1))
    return quantity


#########################################################
###############    Confidence levels   ##################
#########################################################

def _cl_observed(obs, mean, std, cl_desired = 0.95):
    '''(float, float, float):
    Simple function to reverse order of arguments taken by
    _cl_from_erfs so that it can be neatly passed to fsolve
    with the mean as a variable.
    '''
    return _cl_from_erfs(mean, std, obs, cl_desired)

def _cl_from_erfs(mean, std, obs, cl_desired=0.0):
    """(float, float, float) -> float

    Take the mean, standard deviation std, and observed
    upper limit obs; use then equation 3.22 on page 11
    of 1202.3415 to compute the confidence level, which is
    0.95 on this page.

    The mean and standard deviation are the parameters
    of the Gaussian likelihood approximation.

    Returns the confidence level. This should be a value
    between 0 and 1.

    If a desired confidence level is known, the outcome will
    be subtracted from the result so as to return 0 if this
    level is reached. This can be of help in computing roots
    to inversely find the mean.

    >>> (math.erf(-2) + math.erf(2)) / (1 + math.erf(2))
    0.0
    >>> _cl_from_erfs(2, 1/math.sqrt(2), 0)
    0.0
    >>> denom = math.sqrt(2)
    >>> (math.erf(1/denom))/(1+math.erf(1/denom))
    0.4057132913274699
    >>> _cl_from_erfs(1, 1, 1)
    0.4057132913274699

    """
    # See equation 3.22 on page 11 of 1202.3415 for details.
    denom = std * math.sqrt(2)
    cl = math.erf((obs - mean)/denom) + math.erf(mean/denom)
    division = (1 + math.erf(mean/denom))
    if division == 0:
        return None
    return cl / division - cl_desired


#########################################################
####   Computing the mean and standard deviation  #######
#########################################################
def _compute_mean_std(obs, exp=None, obs_err = 0.0, exp_err = 0.0, lumi = 0.0,
        lumi_err = 0.0, start=0.0, tol=0.05):
    """(float, [float, float, float, float, float, float, float])
                -> float, float

    Compute the mean and standard deviation for a Gaussian
    likelihood approximation from the given observed and
    expected upper limits obs and exp.

    The upper limit cross sections are unitless but must be in the
    same units;
    the luminosity lumi must be given in these same units inverted.

    Take the expected upper limit on the cross section exp to be equal to
    the observed upper limit obs when this quantity is absent.

    To solve the error function expressions for the mean, use the starting
    value start.
    When more than one solution is found, make sure that one is returned
    that yields a 95% confidence level within the level of tolerance tol.

    Return the mean and the standard deviation.

    Future: Return errors on mean and standard deviation.

    >>> _compute_mean_std(0, 0, 0)
    (None, None)
    >>> mean, std = mean_std(1, 1.96, lumi=20000)
    >>> (round(mean, 5), std)
    (-2.18733, 1.0)

    """

    # If no expected upper limit is given, take it to be equal
    # to the observed upper limit:
    if type(exp) == type(None):
        exp = obs

    # We approximate std = std_obs = std_exp
    std = exp / 1.96 # standard deviation
    # Compute here in the future the error on std.

    # Take care not to divide by zero in the error function:
    # There one finds a division by std_obs = std_exp = exp/1.96
    if exp == 0 or std == 0:
        logging.warning('Zero expected upper limit or std.' +
                ' Will not calculate' +
                ' mean and standard deviation.' +
                ' obs=' + str(obs) + ', exp='  + str(exp) + ', lumi='
                + str(lumi))
        return None, None

    # The mean can be positive or negative.
    sols = opt.fsolve(_cl_from_erfs, start,\
            args=(std, obs, 0.95), full_output = 1)

    # Do not return the last found solution if fsolve failed:
    if sols[-2] != 1:
        return None, None

    # Take vector solution from the successful fsolve result:
    means = sols[0]

    # This vector is one-dimensional:
    mean = means[0]
    # In the future, compute here the also the error on the mean.
    # In case more than one solution was found, cycle through until
    # the observed and expected upper limits can be obtained from them:
    if len(means) > 1:
        for mu in means:
            if _check_obs_exp_from_mean_std(mean, std, obs, exp, tol):
                mean = mu

    return mean, std

def mean_std(obs, exp=None, obs_err = 0.0, exp_err = 0.0, lumi = 20000.0,\
        lumi_err = 0.0, start=0.0, tol=0.1):
    """(float, [float, float, float, float, float, float, float])
                    -> float, float

    The cross sections are unitless but must be in the same units;
    the luminosity lumi must be given in these same units inverted.

    Take the expected upper limit on the cross section exp to be equal to
    the observed upper limit obs when this quantity is absent.

    To solve the error function expressions for the mean, use the starting
    value start.

    Check that the desired confidence level (0.95) and original observed
    (obs) and expected (exp) limits are obtained from the new mean and
    standard deviation within a degree of tolerance tol.

    Return the mean and the standard deviation:
    The mean determines the maximum of the Gaussian approximation of the
    likelihood. It measures the deviation of the oberved events from the
    background compared to that of signal and background.

    >>> mean_std(0, 0)
    (None, None)
    >>> mean, std = mean_std(1, 1.96, lumi=20000) # pb, pb, pb-1
    >>> (round(mean, 5), std)
    (-2.18733, 1.0)

    """

    # If no expected upper limit is given, take it to be equal
    # to the observed upper limit:
    if type(exp) == type(None):
        exp = obs

    mean, std = _compute_mean_std(obs, exp=exp, obs_err = obs_err,\
            exp_err = exp_err, lumi = lumi, lumi_err = lumi_err,\
            start= start, tol=tol)
    if mean == None or std == None:
        return None, None

    # Check that mean, std
    #       1. yield a 95% confidence level and
    #       2. yield the original observed and expected upper limits;
    # if not, return None types:
    if not _check_cl95_from_mean_std(mean, std, obs, tol=0.01) or\
            not _check_obs_exp_from_mean(mean, std, obs, exp, tol=tol):
        # Change the start value to a larger or negative value if necessary:
        start = obs if obs >= exp else -obs
        # Compute the mean and std again:
        mean, std = _compute_mean_std(obs, exp=exp, obs_err = obs_err,\
                exp_err = exp_err, lumi = lumi, lumi_err = lumi_err,\
                start=start, tol=tol)
        # Issue a warning if the 95% confidence level still cannot be reached:
        if not _check_cl95_from_mean_std(mean, std, obs, tol=0.01):
            logging.warning('Desired confidence level not reached.' +\
                    ' Could not obtain mean and std for obs=' + str(obs) +\
                    ', exp=' + str(exp) + ', cl=' +\
                    str(_cl_from_erfs(mean, std, obs)))
            return None, None
        # Issue a warning if the original observed limit still
        # cannot be obtained:
        if not _check_obs_exp_from_mean(mean, std, obs, exp, start=start,
                tol=tol):
            logging.warning('Desired observed and/or expected upper limit' +\
                    ' could not be obtained from mean=' + str(mean) +\
                    ', std=' + str(std) + '; original obs=' + str(obs) +\
                    ', original exp=' + str(exp) +\
                    ', computed from mean, std were (obs, exp):' +\
                    str(obs_exp(mean, std)))
            return None, None

    # Check the validity of the Gaussian approximation:
    _check_gaussian_approx(mean, std, lumi)

    ## NB How to get error on mean?
    # Want to return errors also!
    return mean, std



def mn_sd(obs, exp=None, obs_err = 0.0, exp_err = 0.0, lumi=20000.0,\
        lumi_err = 0.0):
    """(float/Unum, float/Unum[, float/Unum, float/Unum, float/Unum,
            float/Unum]) -> float/Unum, float/Unum

    Wrapper for mean_std.
    Take upper limit(s) and theory prediction, as well as the luminosity,
    and calculate the mean, standard deviation, and thereafter chi2.

    Make sure input is all in same units if no units given.

    """
    # Convert units to unitless numbers if necessary:
    obs = _in_picobarn(obs)
    obs_err = _in_picobarn(obs_err)
    exp = _in_picobarn(exp)
    exp_err = _in_picobarn(exp_err)
    lumi = _in_inverse_picobarn(lumi)
    lumi_err = _in_inverse_picobarn(lumi_err)

    if lumi == 0:
        # needed for estimate of validity of approximation
        lumi = _in_picobarn(20.0*fb)
    mean, std = mean_std(obs, exp, obs_err, exp_err, lumi, lumi_err)

    if mean == None or std == None:
        return None, None

    return (mean*pb, std*pb)


def _mean_std_direct(lumi, b, n):
    """(float, int, int) -> float, float
    Calculate mean and standard deviation directly
    from the numbers of events observed (n), predicted
    by theory (s), and predicted by background (b) as
    they are by construction in 1202.3415.
    """
    mean = (n - b)/(1.0*lumi)
    std = math.sqrt(n)/(1.0*lumi)
    return mean, std



#########################################################
####################    Checks   ########################
#########################################################

def _check_obs_exp_from_mean(mean, std, obs, exp, start=0.0, tol=0.1):
    """(float, float, float, float, float) -> bool
    Take the mean and standard deviation std and check
    if they yield the observed upper limit obs within
    a certain tolerance tol.

    """
    o, e = obs_exp(mean, std, start)
    if (1 - tol)*obs < o < (1 + tol)*obs and (1 - tol)*exp < e < (1 + tol)*exp:
        return True
    return False


def _check_cl95_from_mean_std(mean, std, obs, tol=0.01):
    """(flaot, float, float) -> bool
    Take the mean and standard deviation std and check
    whether they yield a 95% confidence level for the given
    observed upper limit obs within the tolerance tol.
    """

    if 0.95 - tol < round(_cl_from_erfs(mean, std, obs), 2) < 0.95 + tol:
        return True
    return False

def _background_fluctuations(mean, std, lumi=20000.0):
    """(float, float, float) -> float
    Return the fluctuations compared to the background.
    that is, (nobs - nb) / nb.

    """
    return 1/((std**2*lumi/mean) - 1)

def _n_observed(std, lumi=20000.0):
    """
    Return the number of observed events as
    defined by the mean and standard deviation
    of the Gaussian approximation of the likelihood.
    """
    return std**2 * lumi**2


def _check_gaussian_approx(mean, std, lumi=20000.0):
    """(float, float, float) -> bool
    Check that there were only small fluctuations compared to the background:
    This means (n_obs - b)/b <= 2sigma.

    Check also that the number of observed events is high enough for the
    Gaussian approximation:
    This means n_obs > 10.
    """
    # If lumi was not given, it may be zero, which always leads to
    # large fluctuations from the background.
    # In this case, set luminosity to 20, assuming inverse femtobarn.
    if lumi == 0.0:
        lumi = _in_inverse_picobarn(20.0*fb)
    small_fluct = True

    # Check fluctuations compared to background:
    # Note that
    # (n_obs - b)/b <= 2sigma
    # can be rewritten as
    # 1/((std**2*lumi/mean) -1) <= 2sigma
    # because std = sqrt(nobs)/lumi and mean=nobs-b/lumi by construction
    if abs(1/((std**2*lumi/mean) - 1)) > abs(2*std):
        nobs = std**2 * lumi**2
        b = nobs - mean*lumi
        # Issue a warning if large fluctuations compared to background occur:
        logging.warning('The fluctuations compared to the background were' +
                ' large (> 20%).\nstd=' + str(std) +', lumi=' + str(lumi)
                + ', mean='+ str(mean) + ', so nobs=' + str(nobs) + ', b=' +
                str(b) + '.\nThe Gaussian approximation of the likelihood' +
                ' used for computing chisquares may not be valid here.')
        small_fluct = False

    # Check number of observed events:
    # Note that
    # n_obs >= 10
    # can be rewritten as
    # std**2 * lumi >= 10,
    # because std = sqrt(nobs)/lumi by construction.
    if std**2 * lumi**2 < 10:
        nobs = std**2 * lumi**2
        # Issue a warning if the number of observed events is low:
        logging.warning('The number of observed events is estimated to be' +
                ' low here (< 10 events). The Gaussian approximation of the' +
                ' likelihood used for computing chisquares may not be valid' +
                ' here.' +
                ' \nstd=' + str(std) +', lumi=' + str(lumi)
                + ', mean='+ str(mean) + ', so nobs=' + str(nobs))
        small_fluct = False
    return small_fluct


#########################################################
##########    Computing likelihoods and chi2  ###########
#########################################################
def gaussian_likelihood(obs, exp=None, obs_err = 0.0, exp_err = 0.0,\
        lumi=20000.0, lumi_err = 0.0, x=np.linspace(0, 10, num=1000)):
    '''(float, float, [float, float, float, array]) -> array

    Make sure values have same cross section/inverse cross section units.

    This returns the gaussian distribution with mean and standard deviation
    obtained from mean_std above.
    For more info see this function.

    The distribution returned is a scipy.stats.norm.pdf, a probability
    density function from a normal (Gaussian) with the given mean and
    standard deviation obtained from mean_std.
    It returns the function over the given x-range.

    '''
    mean, std = mean_std(obs, exp, obs_err, exp_err, lumi, lumi_err)
    return stats.norm.pdf(x, loc=mean, scale=std)


def chi_squared(mean, std, theo):
    """(float, float, float) -> float

    Take mean and standard deviation std of a Gaussian,
    and return the mininum chi squared corresponding to
    this distribution.

    Check that mean, std are 73, 221 for obs=484 fb,
    exp=433 fb, lumi = 20 fb-1.
    Then check that for theo=344fb is (344-73)^2/221^2
    = 1.5

    >>> chi_squared(1, 1, 1)
    0.0
    >>> chi_squared(1, 0.5, 2)
    4.0
    >>> mean, std = mean_std(484, 433, lumi=20) # fb, fb, fb-1
    >>> round(chi_squared(mean, std, 344),1)
    1.5
    """
    # As in fittino: chi2 = (O_measured - O_predicted)**2/sigma**2
    return (float(theo) - mean)**2/float(std)**2

def chi2(theo, obs, exp=None, obs_err=0.0, exp_err=0.0,\
        lumi=20000.0, lumi_err=0.0):
    '''(float/Unum, float/Unum[, float/Unum, float/Unum, float/Unum,
                float/Unum])->float)

    Take the observed upper limit obs on the cross section as
    given by experiment and the theory prediction theo of the cross
    section, as well as the expected upper limit exp on the cross
    section. From that, calculate the parameters for a Gaussian
    likelihood distribution for this upper limit.

    Wrapper for chi_squared and mean_std functions.
    Take observed (and expected) upper limit(s) as well any given
    luminosity lumi, and calculate the mean and standard deviation.
    After that, compute the chi2 with the mean, standard deviation,
    and theory prediction.

    Make sure input is all in same units if no units given.

    Check again that mean, std are 73, 221 for obs=484 fb,
    exp=433 fb, lumi = 20 fb-1.
    Then check that for theo=344fb is (344-73)^2/221^2
    = 1.5

    >>> round(chi2(344, 484, 433),1)
    1.5
    >>> round(chi2(344*fb, 484*fb, 433*fb), 1)
    1.5

    '''
    if not isinstance(theo, Unum) and theo == None:
        return None
    if not isinstance(obs, Unum) and obs == None:
        return None
    

    mean, std = mn_sd(obs, exp, obs_err, exp_err, lumi, lumi_err)

    # Convert units to unitless numbers if necessary:
    theo = _in_picobarn(theo)

    # If no mean can be found, return None:
    if not isinstance(mean, Unum) and mean == None:
        return None

    return chi_squared(mean.asNumber(pb), std.asNumber(pb), theo)



if __name__ == '__main__':

    import doctest
    doctest.testmod()

