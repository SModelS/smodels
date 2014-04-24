"""
.. module:: experiment.smsInterpolation
   :synopsis: TODO: write synopsis
   
   smsInterpolation is called by smsResults.getSmartUpperLimit. UpperLimit
   takes arbitrary input masses and checks if there is a corresponding upper
   limit for the given analysis and topology. The upper limit is returned in
   'pb'. If several histograms with different x-values are available, an
   interpolation is performed.

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""

from . import smsResults
from . import smsHelpers
import numpy as np
from scipy.interpolate import griddata
from tools.physicsUnits import rmvunit
from tools.physicsUnits import addunit
import logging

logger = logging.getLogger(__name__) # pylint: disable-msg=C0103


def upperLimit(analysis, topology, masses, run=None):
    """
    Return upper limit for analysis-topology for given masses. 
    
    :param masses: list of masses, with (mother, intermediate(s), LSP). For
    intermediate masses: if possible do interpolation over upper limits for
    different x-values. If interpolation is not possible: check if masses are
    comparable to the assumptions in the histogram.
    
    """
    d = smsResults.getaxes(analysis, topology)
    if not run:
        run = smsHelpers.getRun(analysis)
    if not d:
        logger.error("%s/%s not found." % (analysis, topology))
        return None
    if len(masses) == 2 and not d[0]['mz']:
        return smsResults.getUpperLimit(analysis, topology,
                                        masses[_getAxis('x', d[0]['axes'])],
                                        masses[_getAxis('y', d[0]['axes'])],
                                        interpolate=True)
    if len(masses) == 2 and d[0]['mz']:
        logger.error("Need intermediate mass input for %s/%s." \
                         % (analysis, topology))
        return None
    if len(masses) > 2 and not d[0]['mz']:
        logger.error("No intermediate mass in %s/%s." % (analysis, topology))
        return None
    if len(masses) > 3 or len(d[0]['mz']) > 1:
        logger.error("More than one intermediate mass in %s/%s. Cannot find "
                     "upper limit for topologies with more than one "
                     "intermediate mass." % (analysis, topology))
        return None
    if len(masses) > 2 and len(d) == 1:
        if _compareMasses(masses, d[0]):
            logger.error("Only one histogram available for %s/%s, cannot "
                         "interpolate for intermediate mass." \
                         % (analysis, topology))
            return smsResults.getUpperLimit(analysis, _getHistName(topology,
                                         d[0]['mz'][0]),
                                         masses[_getAxis('x', d[0]['axes'])],
                                         masses[_getAxis('y', d[0]['axes'])],
                                         interpolate=True)
        return None
    return _doGridData(analysis, topology, masses, d, run)


def _getHistName(topo, mz):
    """
    Build histogram name for given topology and mz information.
    
    :param mz: given in the axes-information
    
    """
    if mz == None:
        return topo
    elif 'D' in mz:
        return topo + 'D' + mz.split('=')[1]
    else:
        return topo + mz


def _doGridData(analysis, topology, masses, dPar, run=None):
    """
    Create np.array and uses scipy.griddata function for analysis-topology.
    
    """
    masslist = []
    ullist = []

    for ds in dPar:
        if not ds['mz']:
            logger.error("No information on intermediate mass available for "
                         "%s/%s." % (analysis, topology))
            return None

        d = None
        l = None
        m1 = None

        if ds['mz'][0].find('D') > -1:
            d = float(ds['mz'][0].split('=')[1])
        elif ds['mz'][0].find('LSP')>-1:
            l = float(ds['mz'][0].split("LSP")[1])
        elif ds['mz'][0].find('M1')>-1:
            m1 = float(ds['mz'][0].split("M1")[1])

        ulDict = smsHelpers.getUpperLimitDictionary(analysis,
                                        _getHistName(topology, ds['mz'][0]), run)
        for x in ulDict:
            for y in ulDict[x]:

                if d:
                    massv = [0., 0., 0.]
                    massv[_getAxis('x', ds['axes'])] = x
                    massv[_getAxis('y', ds['axes'])] = y
                    if massv[_getAxis('x', ds['mz'][0])] == 0.:
                        massv[_getAxis('x', ds['mz'][0])] = massv[_getAxis('y',
                                                            ds['mz'][0])] + d
                    if massv[_getAxis('y', ds['mz'][0])] == 0.:
                        massv[_getAxis('y', ds['mz'][0])] = massv[_getAxis('x',
                                                            ds['mz'][0])] - d

                elif l:
                    massv = [0., 0., l]
                    massv[_getAxis('x', ds['axes'])] = x
                    massv[_getAxis('y', ds['axes'])] = y

                elif m1:
                    massv = [m1, 0., 0.]
                    massv[_getAxis('x', ds['axes'])] = x
                    massv[_getAxis('y', ds['axes'])] = y


                else:
                    massv = [x, _getxval(x, y, ds['mz'][0], mass=True), y]

                masslist.append(massv)
                ullist.append(ulDict[x][y])

    p = np.array(masslist)
    v = np.array(ullist)

    mx = rmvunit(masses[0], "GeV")
    my = rmvunit(masses[1], "GeV")
    mz = rmvunit(masses[2], "GeV")

    r = griddata(p, v, (mx, my, mz), method = "linear")

    if np.isnan(r):
        logger.error("Masses out of range for %s/%s (no extrapolation)" \
                      % (analysis, topology))
        return None

    return addunit(float(r), 'pb')


def _getAxis(w, a):
    """
    Find according index in the masses-list for w == x, y.
    
    Use the axes-information a == (mx - my).
    
    """
    ml = []
    ml.append(a.find('M1'))
    ml.append(a.find('M2'))
    ml.append(a.find('M0'))
    if w == 'x':
        return _getIndex(ml, second=True)
    if w == 'y':
        return _getIndex(ml)


def _getxval(mx, my, mz, mass=False):
    """
    Calculate x-value for one point.
    
    If mass == True is selected, return intermediate mass instead of x-value.
    
    :param mx: Mother-mass
    :param my: LSP-mass
    :param mz: information on the intermediate mass as given in the
    axes-information.
    
    """
    mx = rmvunit(mx, "GeV")
    my = rmvunit(my, "GeV")
    mz = rmvunit(mz, "GeV")
    if mz.find('x') == -1 and mz.find('C') == -1 and mz.find('y') == -1:
        xfac = float(mz)/100
        if mass:
            return xfac*mx + (1 - xfac)*my
        return float(mz)/100
    if mz.find('x') > -1:
        tx = float(mz[mz.find('x') + 1:mz.find('x') + 4])
    else:
        tx = None
    # if mz.find('y') > -1: ty = float(mz[mz.find('y') + 1])
    # else: ty = None
    if mz.find('C') > -1:
        c = float(mz.split("C")[1])
    else:
        c = None
    z = 0.
    if tx:
        z += tx*my/100
    # if ty: z += ty*mx
    if c:
        z += c
    if mass:
        return z
    xval = (z - my)/(mx - my)
    return xval
    

def _compareMasses(masses, d):
    """
    Check if input masses are comparable to masses in the histogram.
    
    Masses are comparable corresponding to the information given in
    axes-dictionary d.
    
    """
    # Check if histogram axes are M1, M0, return 1 if x-value of histogram is
    # comparable to x value for given masses, 0 if not
    try: 
        x1 = _getxval(masses[0], masses[-1], d['mz'][0])
        x2 = float(rmvunit(masses[1], "GeV") - rmvunit(masses[-1], "GeV"))/ \
                (rmvunit(masses[0], "GeV") - rmvunit(masses[-1], "GeV"))
        if abs(x1 - x2)/(x1 + x2) < 0.1:
            return True
        else:
            return None
    except: # TODO: which exception?
        # Check if histogram for fixed LSP mass, return 1 if my is comparable
        # to LSP mass of the histogram, 0 if not
        if d['mz'][0].find('LSP') > -1:
            mlsp = float(d['mz'][0].split("LSP")[1])
            if abs(mlsp - rmvunit(masses[-1], "GeV"))/ \
                    (mlsp + rmvunit(masses[-1], "GeV")) < 0.1:
                return True
            else:
                return None
        # Check for fixed deltaM
        elif d['mz'][0].find('D') > -1: 
            ml = []
            ml.append(d['mz'][0].find('M1'))
            ml.append(d['mz'][0].find('M2'))
            ml.append(d['mz'][0].find('M0'))
            deltam = float(d['mz'][0].split('=')[1])
            deltain = rmvunit(masses[_getIndex(ml, second=True)], "GeV") - \
                    rmvunit(masses[_getIndex(ml)], "GeV")
            if deltain < 0:
                return None
            if abs(deltain-deltam)/(deltain+deltam) < 0.1:
                return True
            else:
                return None
        # Check for fixed m_mother
        elif d['mz'][0].find('M1') > -1:
            mmother = float(d['mz'][0].split("M1")[1])
            if abs(mmother - rmvunit(masses[0], "GeV")) / \
                    (mmother + rmvunit(masses[0], "GeV")) < 0.1:
                return True
            else:
                return None
            

def _getIndex(ls, second=False):
    """
    Find index of list element with maximum value.
    
    If the last element is the maximum, return -1. If second == True, find
    second largest list element.
    
    """
    ind = np.argmax(ls)
    if second:
        v = 0
        for i in range(len(ls)):
            if not i == ind:
                if ls[i] >= v:
                    v = ls[i]
                    ind2 = i
        if ind2 == len(ls)-1:
            return -1
        return ind2
    if ind == len(ls)-1:
        return -1
    return ind
