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
from smodels.experiment.finalStateParticles import finalStates


#Get all finalStateLabels
finalStateLabels = finalStates.getValuesFor('label')


def distance(ul1, ul2, el1=None, el2=None):
    """
    Define distance between two values in upper limit space.
    The distance is defined as d = 2*|ul1-ul2|/(ul1+ul2).
    If el1 and el2 are Element objects, the distance
    includes the relative differences between their masses
    and widths: d = max([d_ul,d_mass,d_width])
    
    
    :parameter ul1: upper limit value (in fb) for element1
    :parameter ul2: upper limit value (in fb) for element2
    :parameter el1: Element1
    :parameter el2: Element2

    :returns: relative distance
    """
    if ul1 is None or ul2 is None:
        return None
    ulDistance = 2.*abs(ul1 - ul2)/(ul1 + ul2)
    if ulDistance < 0.:
        # Skip masses without an upper limit
        return None
    if (not el1 is None) and (not el2 is None):
        mass1 = [x.asNumber(GeV) for x in _flattenList(el1.getMasses())]
        mass2 = [x.asNumber(GeV) for x in _flattenList(el2.getMasses())]
        width1 = [x.totalwidth.asNumber(GeV) for x in _flattenList(el1.getBSMparticles())]
        width2 = [x.totalwidth.asNumber(GeV) for x in _flattenList(el1.getBSMparticles())]
        massDistance = max([0.]+[2*abs(m1-mass2[i])/(m1+mass2[i]) for i,m1 in enumerate(mass1) if m1 != mass2[i]])
        widthDistance = max([0.]+[2*abs(w1-width2[i])/(w1+width2[i]) for i,w1 in enumerate(width1) if w1 != width2[i]])
    else:
        massDistance = 0.
        widthDistance = 0.

    distanceValue = max([ulDistance,massDistance,widthDistance])
    
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

    # Evaluate the inequality for each cross section info
    result = crossSection.XSectionList()
    for info in infoList:
        res = 0.
        xsecRes = crossSection.XSection()
        xsecRes.info = info
        for weightA in weights:
            for weightB in weights:
                a = weightA.getXsecsFor(info.label)[0].value.asNumber ( fb )
                b = weightB.getXsecsFor(info.label)[0].value.asNumber ( fb )
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
        # If there are no cross sections, can not evaluate
        return 'N/A'
    zeros = crossSection.XSectionList(infoList)
    weightA.combineWith(zeros)
    weightB.combineWith(zeros)

    # Evaluate the inequality for each cross section info
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


def elementsInStr(instring,removeQuotes=True):
    """
    Parse instring and return a list of elements appearing in instring.
    instring can also be a list of strings.
    
    :param instring: string containing elements (e.g. "[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]")
    :param removeQuotes: If True, it will remove the quotes from the particle labels.
                         Set to False, if one wants to run eval on the output.
    
    :returns: list of elements appearing in instring in string format
    
    """
    outstr = ""
    if isinstance(instring,str):
        outstr = instring
    elif isinstance(instring,list):
        for st in instring:
            if not isinstance(st,str):
                logger.error("Input must be a string or a list of strings")
                raise SModelSError()
            # Combine list of strings in a single string
            outstr += st
    else:
        raise SModelSError ( "syntax error in constraint/condition: ``%s''." \
              "Check your constraints and conditions in your database." % str(instring) )

    elements = []
    outstr = outstr.replace(" ", "")
    if removeQuotes:
        outstr = outstr.replace("'", "")
    elStr = ""
    nc = 0
    # Parse the string and looks for matching ['s and ]'s, when the matching is
    # complete, store element
    for c in outstr:
        delta = 0
        if c == '[':
            delta = -1
        elif c == ']':
            delta = 1
        nc += delta
        if nc != 0:
            elStr += c
        if nc == 0 and delta != 0:
            elements.append(elStr + c)
            elStr = ""
            # Syntax checks
            ptclist = elements[-1].replace(']', ',').replace('[', ',').\
                    split(',')       
            for ptc in ptclist:
                ptc = ptc.replace("'","")
                if not ptc:
                    continue          
                if not ptc in finalStateLabels:
                    raise SModelSError("Unknown particle. Add " + ptc + " to finalStateParticles.py")

    # Check if there are not unmatched ['s and/or ]'s in the string
    if nc != 0:
        raise SModelSError("Wrong input (incomplete elements?) " + instring)

    return elements


def sortParticleList(ptcList):
    """
    sorts a list of particle or particle list objects by their label
    :param ptcList: list to be sorted containing particle or particle list objects
    :return: sorted list of particles
    """
    
    newPtcList = sorted(ptcList, key=lambda x: x.label) 
    return newPtcList    


def getValuesForObj(obj, attribute):
    """
    Loops over all attributes in the object and in its attributes
    and returns a list of values for the desired attribute:
    
    :param obj: Any object with a __dict__ attribute
    :param attribute: String for the desired attribute
    
    :return: List with unique attribute values. If the attribute is not found, returns empty list.
    """
    
    values = []
    try:
        objDict = obj.__dict__.items()
    except:
        return values        
    
      
    for attr,value in objDict:
        if attribute == attr:
            values += [value]
        elif isinstance(value,Iterable):
            values += [getValuesForObj(v,attribute) for v in value]
        else:
            values += getValuesForObj(value,attribute)
    
    values =  list(filter(lambda a: (not isinstance(a,list)) or a != [], values))
    values = _flattenList(values)
    uniqueValues = [v for n,v in enumerate(values) if v not in values[:n]]
    
    return uniqueValues


def getAttributesFrom(obj):
    """
    Loops over all attributes in the object and return a list
    of the attributes.
    
    :param obj: Any object with a __dict__ attribute
    
    :return: List with unique attribute labels.
    """
    
    attributes = []
    try:
        objDict = obj.__dict__.items()
    except:
        return attributes        
    
      
    for attr,value in objDict:
        attributes.append(attr)
        if isinstance(value,list):
            attributes += [getAttributesFrom(v) for v in value]
        elif isinstance(value,dict):
            attributes += [getAttributesFrom(v) for v in value.values()]
        else:
            attributes += getAttributesFrom(value)
    
    attributes =  list(filter(lambda a: a != [], attributes))    
    
    return list(set(_flattenList(attributes)))

