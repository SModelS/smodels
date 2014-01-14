#!/usr/bin/env python

"""
.. module:: auxiliaryFunctions
    :synopsis: A collection of functions used to evaluate fuzzy the conditions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from Tools.PhysicsUnits import rmvunit
import crossSection, element
import logging
import copy, itertools
from ParticleNames import Reven, PtcDic
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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
    result =  crossSection.XSectionList(infoList)
    for weight in weights: weight.combineWith(zeros)
    for info in infoList:
        res = 0.        
        xsecRes = result.getXsecsFor(info.label).XSections[0]
        for weightA in weights:
            for weightB in weights:
                a = rmvunit(weightA.getXsecsFor(info.label).XSections[0].value,'fb')
                b = rmvunit(weightB.getXsecsFor(info.label).XSections[0].value,'fb')
                if a + b == 0.: continue
                res = max(res,abs(a-b)/abs(a+b))                  
        xsecRes.value = res
  
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
    
    zeros = crossSection.XSectionList(infoList)
    result =  crossSection.XSectionList(infoList)
    weightA.combineWith(zeros)
    weightB.combineWith(zeros)
    for info in infoList:
        a = rmvunit(weightA.getXsecsFor(info.label).XSections[0].value,'fb')
        b = rmvunit(weightB.getXsecsFor(info.label).XSections[0].value,'fb')
        xsecRes = result.getXsecsFor(info.label).XSections[0]
        if a + b == 0.: xsecRes.value = 'N/A'
        else: xsecRes.value = (abs(a-b) - (a-b))/(2.*(a+b))
        
    return result

def elementsInStr(instring):
    """ Parses instring (or a list of strings) and return a list of elements (in string format) appearing in instring"""
  
    if type(instring) == type('st'):
        outstr = copy.deepcopy(instring)
    elif type(instring) == type([]):
        outstr = ""
    for st in instring:
        if type(st) != type('st'):
            print "getlements: Input must be a string or a list of strings"
            return False    
        outstr += st
  
    elements = []
    outstr = outstr.replace(" ","")
    while "[[[" in outstr:  #String has element
        st = outstr[outstr.find("[[["):outstr.find("]]]")+3] #Get duplet
        element = element.Element ( st )
        ptclist = element.allParticles()
        #Syntax checks:
        for ib in range(2):
            for ptcL in ptclist[ib]:
                for ptc in ptcL:
                    if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                        print "EvalRes: Unknown particle ",ptc
                        return False
        outstr = outstr.replace(st,"")  # delete element
        elements.append(st)   #Store elements
      
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