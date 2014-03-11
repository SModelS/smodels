"""
.. module:: theory.auxiliaryFunctions
   :synopsis: A collection of functions used to evaluate fuzzy the conditions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from functools import wraps
import crossSection
import copy
import itertools
from experiment import LimitGetter
from tools.PhysicsUnits import addunit, rmvunit
import numpy as np
from scipy import stats
from collections import Iterable
import logging

logger = logging.getLogger(__name__)


def memoize(func):
    """
    A wrapping to store in cache the results of massPosition, since this is
    an expensive function.
    
    """    
    cache = {}
    @wraps(func)
    def wrap(*args):
        """
        missing
        
        """
        if str(args) not in cache:
            cache[str(args)] = func(*args)
        return cache[str(args)]
    return wrap


@memoize
def massPosition(mass, analysis):
    """
    Gives the mass position in upper limit space, using the analysis
    experimental limit data. If nounit=True, the result is given as number
    assuming fb units.
    
    """    
    xmass = LimitGetter.GetPlotLimit(mass, analysis, complain=False)
    if type(xmass) != type(addunit(1.,'pb')):
        return None
    xmass = rmvunit(xmass,'fb')    
    return xmass


def distance(xmass1, xmass2):
    """
    Definition of distance between two mass positions.
    
    """
    if xmass1 is None or xmass2 is None:
        return None
    distance = 2.*abs(xmass1-xmass2)/(xmass1+xmass2)
    if distance < 0.:
        return None # Skip masses without an upper limit
    return distance


def massAvg(massList, method='harmonic'):
    """
    Computes the average mass of massList, according to method (harmonic or
    mean). If massList contains a zero mass, switch method to mean.
    
    """

    # Checks:    
    if not massList:
        return massList
    if len(massList) == 1:
        return massList[0]               
    flatList = [rmvunit(mass, 'GeV') for mass in flattenList(massList)]
    if method == 'harmonic' and 0. in flatList:
        method = 'mean'
    
    for mass in massList:
        if len(mass) != len(massList[0]) \
                or len(mass[0]) != len(massList[0][0]) \
                or len(mass[1]) != len(massList[0][1]):  
            logger.error('Mass shape mismatch in mass list:\n'+str(mass)
                         +' and '+str(massList[0]))
            return False
    
    avgmass = massList[0][:]
    for ib, branch in enumerate(massList[0]):
        for ival in enumerate(branch):
            vals = [rmvunit(mass[ib][ival[0]],'GeV') for mass in massList]
            if method == 'mean':
                avg = np.mean(vals)
            elif method == 'harmonic':
                avg = stats.hmean(vals)
            avgmass[ib][ival[0]] = addunit(float(avg),'GeV')   
    return avgmass


def Csim(*weights):
    return cSim(*weights)


def cSim(*weights):
    """
    Defines the auxiliar similar function.
    
    Returns the maximum relative difference between any element weights of the
    list, normalized to [0,1].
    
    Returns a XSectioList object with the values for each label.
    
    """    
    for weight in weights:
        if type(weight) != type(crossSection.XSectionList()):
            logger.error("Trying to evaluate non-xsection objects")
            return False

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
    result =  crossSection.XSectionList()
    for info in infoList:
        res = 0.       
        xsecRes = crossSection.XSection()
        xsecRes.info = info     
        for weightA in weights:
            for weightB in weights:
                a = rmvunit(weightA.getXsecsFor(info.label)[0].value,'fb')
                b = rmvunit(weightB.getXsecsFor(info.label)[0].value,'fb')
                if a + b == 0.:
                    continue
                res = max(res, abs(a-b)/abs(a+b))                  
        xsecRes.value = res
        result.add(xsecRes)
        
    return result 


def Cgtr(weightA, weightB):
    return cGtr(weightA, weightB)


def cGtr(weightA, weightB):
    """
    Defines the auxiliary greater function returns a number  between 0 and 1
    depending on how much it is violated (0 = A > B, 1 = A << B).
    
    Returns a XSectioList object with the values for each label.
    
    """    
    if type(weightA) != type(crossSection.XSectionList()) or \
            type(weightB) != type(crossSection.XSectionList()):
        logger.error("Trying to evaluate non-xsection objects")
        return False

    # Make sure both xsec lists have the same entries (add zero xsecs for the
    # missing entries)   
    infoList = weightA.getInfo()
    for info in weightB.getInfo():
        if not info in infoList:
            infoList.append(info)
    if not infoList:
        return 'N/A'    # If there are no cross-sections, can not evaluate   
    zeros = crossSection.XSectionList(infoList)    
    weightA.combineWith(zeros)
    weightB.combineWith(zeros)

    # valuate the inequality for each cross-section info    
    result =  crossSection.XSectionList()        
    for info in infoList:
        a = rmvunit(weightA.getXsecsFor(info.label)[0].value,'fb')
        b = rmvunit(weightB.getXsecsFor(info.label)[0].value,'fb')
        xsecRes = crossSection.XSection()
        xsecRes.info = info                        
        if a + b == 0.:
            xsecRes.value = 'N/A'
        else:
            xsecRes.value = (abs(a-b) - (a-b))/(2.*(a+b))
        result.add(xsecRes)
        
    return result

    
def flattenList(inlist, dims=None):
    """
    An auxliary function to completely flatten a multi-dimensional nested
    list. The output ordering is: [first level objects, second level objects,
    ...].    
    
    If dims = [], include dimensions of nested list to it (useful when 
    comparing lists).
    
    """
    flat = []
    for item in inlist:
        if isinstance(item, Iterable) and not isinstance(item, basestring):
            if not dims is None:
                dims.append(len(item))
            for x in flattenList(item, dims):
                flat.append(x)
        else:        
            flat.append(item)             
    return flat    
    
    if type(inlist) != type(list()):
        return inlist
    try:
        flat = list(itertools.chain(*inlist))
    except:
        flat = copy.deepcopy(inlist)
    go = True
    while go:
        go = False
        try:
            flat = list(itertools.chain(*flat))
        except:
            pass
        if len(flat) > 0:    
            for ival, val in enumerate(flat):
                if type(val) ==  type(list()):
                    go = True
                    try:
                        val = list(itertools.chain(*val))
                    except:
                        pass
                    flat.pop(ival)                
                    flat.extend(val)
    return flat
