"""
.. module:: auxiliaryFunctions
   :synopsis: A collection of functions used to evaluate fuzzy the conditions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory import crossSection
from smodels.tools.physicsUnits import pb, GeV, fb
import numpy as np
from scipy import stats
from collections import Iterable
import copy
from smodels.theory.exceptions import SModelSTheoryError as SModelSError

from smodels.tools.smodelsLogging import logger

def massPosition(mass, txdata):
    """ Give mass position in upper limit space.    
    Use the analysis experimental limit data. 
    :param txdata: TxNameData object holding the data and interpolation   
    """
    xmass = txdata.getValueFor(mass)
    if type(xmass) != type(1.*pb):
        return None
    xmass = xmass / fb
    return xmass.asNumber()


def distance(xmass1, xmass2):
    """
    Define distance between two mass positions in upper limit space.
    The distance is defined as d = 2*|xmass1-xmass2|/(xmass1+xmass2).
    
    
    :parameter xmass1: upper limit value (in fb) for the mass1
    :parameter xmass2: upper limit value (in fb) for the mass2
    :returns: relative mass distance in upper limit space     
    """
    if xmass1 is None or xmass2 is None:
        return None
    distanceValue = 2.*abs(xmass1 - xmass2) / (xmass1 + xmass2)
    if distanceValue < 0.:
        # Skip masses without an upper limit
        return None
    
    return distanceValue


def massAvg(massList, method='weighted', weights=None):
    """
    Compute the average mass of massList according to method.

    If method=weighted but weights were not properly defined,
    switch method to harmonic.    
    If massList contains a zero mass, switch method to mean.
    
    :parameter method: possible values: harmonic, mean, weighted
    :parameter weights: weights of elements (only for weighted average)
    
    """
    if not massList:
        return massList
    if massList.count(massList[0]) == len(massList):
        return massList[0]

    if method == 'weighted' and (not weights or len(weights) != len(massList)):
        method = 'harmonic'
    flatList = [ mass / GeV for mass in _flattenList(massList)]
    if method == 'harmonic' and 0. in flatList:
        method = 'mean'

    for mass in massList:
        if len(mass) != len(massList[0]) \
                or len(mass[0]) != len(massList[0][0]) \
                or len(mass[1]) != len(massList[0][1]):
            logger.error('Mass shape mismatch in mass list:\n' + str(mass) +
                         ' and ' + str(massList[0]))
            raise SModelSError()

    avgmass = copy.deepcopy(massList[0])
    for ib, branch in enumerate(massList[0]):
        for ival in enumerate(branch):
            vals = [ float(mass[ib][ival[0]] / GeV) for mass in massList]
            if method == 'mean':
                avg = np.mean(vals)
            elif method == 'harmonic':
                avg = stats.hmean(vals)
            elif method == 'weighted':
                weights = [ float(weight) for weight in weights ]
                avg = np.average(vals,weights=weights)                
            avgmass[ib][ival[0]] = float(avg)*GeV

    return avgmass


def cSim(*weights):
    """
    Define the auxiliar similar function.
    
    Return the maximum relative difference between any element weights of the
    list, normalized to [0,1].
    
    :returns: XSectionList object with the values for each label.
    
    """
    for weight in weights:
        if type(weight) != type(crossSection.XSectionList()):
            logger.error("Trying to evaluate non-xsection objects")
            raise SModelSError()

    # Make sure both xsec lists have the same entries (add zero xsecs for the
    # missing entries)
    infoList = []
    for weight in weights:
        for info in weight.getInfo():
            if not info in infoList:
                infoList.append(info)
    zeros = crossSection.XSectionList(infoList)
    for weight in weights:
        weight.combineWith(zeros)

    # Evaluate the inequality for each cross-section info
    result = crossSection.XSectionList()
    for info in infoList:
        res = 0.
        xsecRes = crossSection.XSection()
        xsecRes.info = info
        for weightA in weights:
            for weightB in weights:
                a = weightA.getXsecsFor(info.label)[0].value / fb
                b = weightB.getXsecsFor(info.label)[0].value / fb
                if a + b == 0.:
                    continue
                res = max(res, abs(a - b) / abs(a + b))
        xsecRes.value = res
        result.add(xsecRes)

    return result


def cGtr(weightA, weightB):
    """
    Define the auxiliary greater function.
    
    Return a number between 0 and 1 depending on how much it is violated
    (0 = A > B, 1 = A << B).
    
    :returns: XSectioList object with the values for each label.
    
    """
    if type(weightA) != type(crossSection.XSectionList()) or \
            type(weightB) != type(crossSection.XSectionList()):
        logger.error("Trying to evaluate non-xsection objects")
        raise SModelSError()

    # Make sure both xsec lists have the same entries (add zero xsecs for the
    # missing entries)
    infoList = weightA.getInfo()
    for info in weightB.getInfo():
        if not info in infoList:
            infoList.append(info)
    if not infoList:
        # If there are no cross-sections, can not evaluate
        return 'N/A'
    zeros = crossSection.XSectionList(infoList)
    weightA.combineWith(zeros)
    weightB.combineWith(zeros)

    # Evaluate the inequality for each cross-section info
    result = crossSection.XSectionList()
    for info in infoList:
        a = weightA.getXsecsFor(info.label)[0].value / fb
        b = weightB.getXsecsFor(info.label)[0].value / fb
        xsecRes = crossSection.XSection()
        xsecRes.info = info
        if a + b == 0.:
            xsecRes.value = 'N/A'
        else:
            xsecRes.value = (abs(a - b) - (a - b)) / (2.*(a + b))
        result.add(xsecRes)

    return result


def _flattenList(inlist, dims=None):
    """
    Flatten a multi-dimensional nested list.
    
    Output ordering: [first level objects, second level objects, ...].    
    
    If dims == [], include dimensions of nested list to it. This is useful when
    comparing lists).
    
    """
    flat = []
    for item in inlist:
        if isinstance(item, Iterable) and not isinstance(item, str ):
            if not dims is None:
                dims.append(len(item))
            for x in _flattenList(item, dims):
                flat.append(x)
        else:
            flat.append(item)
    return flat
    
def index_bisect(inlist, el):
    """
    Return the index where to insert item el in inlist.
    inlist is assumed to be sorted and a comparison function (lt or cmp)
    must exist for el and the other elements of the list.
    If el already appears in the list, inlist.insert(el) will
    insert just before the leftmost el already there.  
    """

    lo = 0    
    hi = len(inlist)
    while lo < hi:
        mid = (lo+hi)//2
        if inlist[mid] < el: lo = mid+1
        else: hi = mid
    return lo
