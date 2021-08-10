"""
.. module:: auxiliaryFunctions
   :synopsis: A collection of functions used to evaluate fuzzy the conditions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory import crossSection
from smodels.tools.physicsUnits import fb, GeV
import unum
try:
    from collections.abc import Iterable
except:
    from collections import Iterable
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
import numpy as np



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
        weight += zeros

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
    weightA += zeros
    weightB += zeros

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

def removeUnits(value,standardUnits):
    """
    Remove units from unum objects. Uses the units defined
    in physicsUnits.standard units to normalize the data.

    :param value: Object containing units (e.g. [[100*GeV,100.*GeV],3.*pb])
    :param standardUnits: Unum unit or Array of unum units defined to
                          normalize the data.
    :return: Object normalized to standard units (e.g. [[100,100],3000])
    """

    if isinstance(standardUnits,unum.Unum):
        stdunits = [standardUnits]
    else:
        stdunits = standardUnits

    if isinstance(value,list):
        return [removeUnits(x,stdunits) for x in value]
    if isinstance(value,tuple):
        return tuple([removeUnits(x,stdunits) for x in value])
    elif isinstance(value,dict):
        return dict([[removeUnits(x,stdunits),removeUnits(y,stdunits)] for x,y in value.items()])
    elif isinstance(value,unum.Unum):
        #Check if value has unit or not:
        if not value._unit:
            return value.asNumber()
        #Now try to normalize it by one of the standard pre-defined units:
        for unit in stdunits:
            y = (value/unit).normalize()
            if not y._unit:
                return value.asNumber(unit)
        raise SModelSError("Could not normalize unit value %s using the standard units: %s"
                       %(str(value),str(standardUnits)))
    else:
        return value

def addUnit(obj,unit):
    """
    Add unit to object.
    If the object is a nested list, adds the unit to all of its elements.

    :param obj: Object without units (e.g. [[100,100.]])
    :param unit: Unit to be added to the object (Unum object, e.g. GeV)
    :return: Object with units (e.g. [[100*GeV,100*GeV]])
    """

    if isinstance(obj,list):
        return [addUnit(x,unit) for x in obj]
    elif isinstance(obj,tuple):
        return tuple([addUnit(x,unit) for x in obj])
    elif isinstance(obj,dict):
        return dict([[addUnit(x,unit),addUnit(y,unit)] for x,y in obj.items()])
    elif isinstance(obj,(float,int,unum.Unum)):
        return obj*unit
    else:
        return obj

def flattenArray(objList):
    """
    Flatten any nested list to a 1D list.

    :param objList: Any list or nested list of objects (e.g. [[[(100.,1e-10),100.],1.],[[200.,200.],2.],..]

    :return: 1D list (e.g. [100.,1e-10,100.,1.,200.,200.,2.,..])
    """

    ret = []

    for obj in objList:
        if isinstance(obj, Iterable) and not isinstance(obj, str):
            ret.extend(flattenArray(obj))
        else:
            ret.append(obj)
    return ret

def reshapeList(objList,shapeArray):
    """
    Reshape a flat list according to the shape of shapeArray.
    The number of elements in shapeArray should equal the length of objList.

    :param objList: 1D list of objects (e.g. [200,100,'*',50,'*'])
    :param shapeArray: Nested array (e.g. [[float,float,'*',float],'*'])

    :return: Array with elements from objList shaped according to shapeArray
             (e.g. [[200.,100.,'*',50],'*'])
    """

    if isinstance(shapeArray,list):
        return [reshapeList(objList,v) for v in shapeArray]
    else:
        return objList.pop(0)

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

    # Check if there are not unmatched ['s and/or ]'s in the string
    if nc != 0:
        raise SModelSError("Wrong input (incomplete elements?) " + instring)

    return elements

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
    except AttributeError:
        return values


    for attr,value in objDict:
        if attribute == attr:
            values += [value]
        elif isinstance(value,Iterable):
            values += [getValuesForObj(v,attribute) for v in value if not v is obj]
        elif not value is obj:
            values += getValuesForObj(value,attribute)

    values =  list(filter(lambda a: (not isinstance(a,list)) or a != [], values))
    values = flattenArray(values)
    uniqueValues = [v for n,v in enumerate(values) if v not in values[:n]]

    return uniqueValues

def getAttributesFrom(obj,skipIDs=[]):
    """
    Loops over all attributes in the object and return a list
    of the attributes.

    :param obj: Any object with a __dict__ attribute
    :param skipIDs: List of object ids. Any object which has its id on the list
                    will be ignored (useful to avoid recursion).

    :return: List with unique attribute labels.
    """

    if id(obj) in skipIDs or isinstance(obj,(tuple,list,float,int,unum.Unum)):
        return []
    else:
        newSkipIDs = skipIDs[:] + [id(obj)]
    attributes = []
    try:
        objDict = obj.__dict__.items()
    except AttributeError:
        return attributes


    for attr,value in objDict:
        attributes.append(attr)
        if isinstance(value,list):
            attributes += [getAttributesFrom(v,newSkipIDs) for v in value]
        elif isinstance(value,dict):
            attributes += [getAttributesFrom(v,newSkipIDs) for v in value.values()]
        else:
            attributes += getAttributesFrom(value,newSkipIDs)

    attributes =  list(filter(lambda a: a != [], attributes))

    return list(set(flattenArray(attributes)))

def roundValue(value,nround=-1):
    """
    Round a value to nround significant digits. If the input value is not
    a float or Unum object, it is zero or infinity, nothing is done.

    :param nround: number of significant digits

    :return: rounded value
    """

    if nround <= 0:
        return value
    elif not isinstance(value,(float,unum.Unum)):
        return value

    #Remove units, if it is a unum object
    if isinstance(value,unum.Unum):
        if not value.asNumber():
            return value
        elif np.isinf(value.asNumber()):
            return value
        unit = value/value.asNumber()
        v = value.asNumber()
    else:
        unit = 1.0
        v = value

    #Round value:
    v_rounded = round(v, nround-int(np.floor(np.log10(abs(v))))-1)
    v_rounded *= unit
    return v_rounded

def average(values,weights=None,nround=-1):
    """
    Compute the weighted average of a list of objects.
    All the objects must be of the same type.
    If all objects are equal returns the first entry of the list.
    Only the average of ints, floats and Unum objects or nested lists
    of these can be computed. If the average can not be computed
    returns None.

    :param values: List of objects of the same type
    :param weights: Weights for computing the weighted average. If None it will assume
                    unit weights.
    :param nround: If greater than zero and the returning attibute is numeric, will round it
                      to this number of significant digits.

    """

    if weights is None:
        weights = [1.]*len(values)
    if not isinstance(values,list):
        raise SModelSError("Values must be a list of objects")
    if not isinstance(weights,list):
        raise SModelSError("Weights must be a list of objects")
    if any(not isinstance(w,(float,int)) for w in weights):
        raise SModelSError("Weights must be a list of integers or floats")
    if len(values) != len(weights):
        return SModelSError("Values and weights must have the same length")

    if any(type(v) != type(values[0]) for v in values):
        raise SModelSError("Can not compute average of distinct types of objects")
    if isinstance(values[0],(float,int,unum.Unum)):
        total = values[0]*weights[0]
        for i,v in enumerate(values[1:]):
            total += v*weights[i+1]
        total = total/sum(weights)
        return roundValue(total,nround)
    elif isinstance(values[0],list):
        ndim = len(values[0])
        if any(len(v) != ndim for v in values):
            raise SModelSError("Can not compute average of lists of distinct lengths")
        res = []
        for i in range(ndim):
            try:
                res.append(average([v[i] for v in values],weights,nround))
            except SModelSError:
                return None
        return res
    else:
        if all(values[0] is v for v in values):
            return roundValue(values[0],nround)
        if all(values[0] == v for v in values):
            return roundValue(values[0],nround)
        return None

def rescaleWidth(width):
    """
    The function that is applied to all widths to
    map it into a better variable for interpolation.
    It grows logarithmically from zero (for width=0.)
    to a large number (machine dependent) for width = infinity.

    :param width: Width value (in GeV) with or without units

    :return x: Coordinate value (float)
    """

    if isinstance(width,unum.Unum):
        w = width.asNumber(GeV)
    else:
        w = width

    minWidth = 1e-30 #Any width below this can be safely considered to be zero
    maxWidth = 1e50 #Any width above this can be safely considered to be infinity
    w = (min(w,maxWidth)/minWidth) #Normalize the width and convert it to some finite number (if not finite)
    if w < 1e-10: #The log function misbehaves for very small values of w (step behavior), so we use log(1+x) = x for x << 1
        return w
    else:
        return np.log(1+w)

def unscaleWidth(x):
    """
    Maps a coordinate value back to width (with GeV unit).
    The mapping is such that x=0->width=0 and x=very large -> width = inf.

    :param x: Coordinate value (float)

    :return width: Width value (in GeV) with unit
    """

    minWidth = 1e-30 #Any width below this can be safely considered to be zero
    maxWidth = 1e50 #Any width above this can be safely considered to be infinity
    with np.errstate(over='ignore'): #Temporarily disable overflow error message
        #The small increase in x is required to enforce unscaleWidth(widthToCoordinae(np.inf)) = np.inf
        width = minWidth*(np.exp(x)-1)
        if width > maxWidth:
            width = np.inf
    return width*GeV

def removeInclusives(massArray,shapeArray):
    """
    Remove all entries corresponding to '*' in shapeArray.
    If shapeArray contains entries = '*', the corresponding entries
    in value will be removed from the output.

    :param massArray: Array to be formatted (e.g. [[200.,100.],[200.,100.]] or [[200.,'*'],'*'],0.4])
    :param shapeArray: Array with format info (e.g. ['*',[float,float]])

    :return: formatted array (e.g. [[200.,100.]] or [[200.]],0.4])
    """

    if shapeArray == '*':
        return None
    elif isinstance(massArray,list):
        if len(shapeArray) != len(massArray):
            raise SModelSError("Input value and data shape mismatch (%s,%s)"
                               %(len(shapeArray),len(massArray)))
        else:
            return [removeInclusives(xi,shapeArray[i]) for i,xi in enumerate(massArray)
                    if not removeInclusives(xi,shapeArray[i]) is None]
    else:
        return massArray

def addInclusives(massList,shapeArray):
    """
    Add entries corresponding to '*' in shapeArray.
    If shapeArray contains entries = '*', the corresponding entries
    in value will be added from the output.

    :param massList: 1D array of floats. Its dimension should be equal to the number
                  of non "*" items in shapeArray (e.g. [200.0,100.0])
    :param shapeArray: 1D array containing the data type and '*'. The
                       values of data type are ignored (e.g. [float,'*',float,'*','*']).

    :return: original array with '*' inserted at the correct entries.
    """

    if isinstance(shapeArray,list):
        return [addInclusives(massList,v) for v in shapeArray]
    elif shapeArray != '*':
        return massList.pop(0)
    else:
        return shapeArray
