#!/usr/bin/env python

"""
.. module:: auxiliaryFunctions
    :synopsis: A collection of functions used to evaluate fuzzy the conditions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from functools import wraps
import crossSection
import logging
import copy, itertools
from ParticleNames import Reven, PtcDic
from Experiment import LimitGetter
from Tools.PhysicsUnits import addunit,rmvunit
import numpy as np
from scipy import stats
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def memoize(func):
    """ A wrapping to store in cache the results of massPosition, since this is an expensive function """
    
    cache = {}
    @wraps(func)
    def wrap(*args):
        if str(args) not in cache: cache[str(args)] = func(*args)
        return cache[str(args)]
    return wrap

@memoize
def massPosition(mass,Analysis):
    """ gives the mass position in upper limit space, using the analysis experimental limit data.
    If nounit=True, the result is given as number assuming fb units """
    
    xmass = LimitGetter.GetPlotLimit(mass,Analysis,complain=False)
    if type(xmass) != type(addunit(1.,'pb')): return None
    xmass = rmvunit(xmass,'fb')    
    return xmass


def distance(xmass1,xmass2):
    """ Definition of distance between two mass positions."""

    if xmass1 is None or xmass2 is None: return None
    distance = 2.*abs(xmass1-xmass2)/(xmass1+xmass2)
    if distance < 0.: return None #Skip masses without an upper limit
    return distance


def massAvg(massList,method='harmonic'):
    """Computes the average mass of massList, according to method (harmonic or mean).
    If massList contains a zero mass, switch method to mean"""

#Checks:    
    if not massList: return massList
    if len(massList) == 1: return massList[0]               
    flatList = [rmvunit(mass,'GeV') for mass in flattenList(massList)]
    if method == 'harmonic' and 0. in flatList: method = 'mean'
    
    for mass in massList:
        if len(mass) != len(massList[0]) or len(mass[0]) != len(massList[0][0]) or len(mass[1]) != len(massList[0][1]):  
            logger.error('[massAvg]: mass shape mismatch in mass list:\n'+str(mass)+' and '+str(massList[0]))
            return False
    
    avgmass = copy.deepcopy(massList[0])
    for ib,branch in enumerate(massList[0]):
        for ival in enumerate(branch):
            vals = [rmvunit(mass[ib][ival[0]],'GeV') for mass in massList]
            if method == 'mean': avg = np.mean(vals)
            elif method == 'harmonic': avg = stats.hmean(vals)
            avgmass[ib][ival[0]] = addunit(float(avg),'GeV')   
    return avgmass    


def Csim(*weights):
    """ Defines the auxiliar similar function
    returns the maximum relative difference between any element weights
    of the list, normalized to [0,1].
    Returns a XSectioList object with the values for each label"""
    
    for weight in weights:
        if type(weight) != type(crossSection.XSectionList()):
            logger.error("[Csim]: Trying to evaluate non-xsection objects")
            return False

#Make sure both xsec lists have the same entries (add zero xsecs for the missing entries)
    infoList = []
    for weight in weights:
        for info in weight.getInfo():
            if not info in infoList: infoList.append(info)    
    zeros = crossSection.XSectionList(infoList)   
    for weight in weights: weight.combineWith(zeros)
 
#Evaluate the inequality for each cross-section info   
    result =  crossSection.XSectionList()
    for info in infoList:
        res = 0.       
        xsecRes = crossSection.XSection()
        xsecRes.info = info     
        for weightA in weights:
            for weightB in weights:
                a = rmvunit(weightA.getXsecsFor(info.label)[0].value,'fb')
                b = rmvunit(weightB.getXsecsFor(info.label)[0].value,'fb')
                if a + b == 0.: continue
                res = max(res,abs(a-b)/abs(a+b))                  
        xsecRes.value = res
        result.add(xsecRes)
  
    return result 


def Cgtr(weightA,weightB):
    """Defines the auxiliary greater function returns a number 
    between 0 and 1 depending on how much it is violated (0 = A > B, 1 = A << B).
    Returns a XSectioList object with the values for each label"""
    
    if type(weightA) != type(crossSection.XSectionList()) or type(weightB) != type(crossSection.XSectionList()):
        logger.error("[Cgtr]: Trying to evaluate non-xsection objects")
        return False

#Make sure both xsec lists have the same entries (add zero xsecs for the missing entries)   
    infoList = weightA.getInfo()
    for info in weightB.getInfo():
        if not info in infoList:infoList.append(info)
    if not infoList: return 'N/A'   #If there are no cross-sections, can not evaluate   
    zeros = crossSection.XSectionList(infoList)    
    weightA.combineWith(zeros)
    weightB.combineWith(zeros)

#Evaluate the inequality for each cross-section info    
    result =  crossSection.XSectionList()        
    for info in infoList:
        a = rmvunit(weightA.getXsecsFor(info.label)[0].value,'fb')
        b = rmvunit(weightB.getXsecsFor(info.label)[0].value,'fb')
        xsecRes = crossSection.XSection()
        xsecRes.info = info                        
        if a + b == 0.: xsecRes.value = 'N/A'
        else: xsecRes.value = (abs(a-b) - (a-b))/(2.*(a+b))
        result.add(xsecRes)       
        
    return result

def elementsInStr(instring):
    """ Parses instring (or a list of strings) and return a list of elements (in string format) appearing in instring"""
  
    if type(instring) == type('st'): outstr = instring
    elif type(instring) == type([]):
        outstr = ""
        for st in instring:
            if type(st) != type('st'):
                logger.error("[elementsInStr]: Input must be a string or a list of strings")
                return False        
            outstr += st    #Combines list of strings in a single string
  
    elements = []
    outstr = outstr.replace(" ","").replace("'","")
    elStr = ""
    nc = 0
#Parses the string and looks for matching ['s and ]'s, when the matching is complete, store element    
    for c in outstr:
        delta = 0
        if c == '[':  delta = -1
        elif c == ']': delta = 1
        nc += delta
        if nc != 0: elStr += c
        if nc == 0 and delta != 0:
            elements.append(elStr+c)
            elStr = ""
            #Syntax checks:
            ptclist = elements[-1].replace(']',',').replace('[',',').split(',')
            for ptc in ptclist:
                if not ptc: continue
                if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                    logger.error("[elementsInStr]: Unknown particle "+ptc)
                    return False

#Check if there are not unmatched ['s and/or ]'s in the string:        
    if nc != 0:
        logger.error("[elementsInStr]: wrong input (incomplete elements?) "+instring)
        return False     
      
    return elements

    
def flattenList(inlist):
    """An auxliary function to completely flatten a multi-dimensional nested list.
    The output ordering is: [first level objects, second level objects,...] """
    
    if type(inlist) != type(list()): return inlist
    try: flat = list(itertools.chain(*inlist))
    except: flat = copy.deepcopy(inlist)
    go = True
    while go:
        go = False
        try: flat = list(itertools.chain(*flat))
        except: pass
        if len(flat) > 0:    
            for ival,val in enumerate(flat):
                if type(val) ==  type(list()):
                    go = True
                    try: val = list(itertools.chain(*val))
                    except: pass
                    flat.pop(ival)                
                    flat.extend(val)
    return flat