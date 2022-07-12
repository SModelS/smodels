"""
.. module:: auxiliaryFunctions
   :synopsis: A collection of functions used to evaluate fuzzy the conditions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.tools.physicsUnits import standardUnits, GeV
import unum
import re
try:
    from collections.abc import Iterable
except (ImportError, ModuleNotFoundError):
    from collections import Iterable

from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
import numpy as np
from collections import OrderedDict, deque


def cSim(*weights):
    """
    Define the auxiliar similar function.

    Return the maximum relative difference between any element weights of the
    list, normalized to [0,1].

    :returns: List of values.

    """

    if any(not isinstance(w, (float, unum.Unum)) for w in weights):
        logger.error("Trying to evaluate with invalid objects %s")
        raise SModelSError()

    x = np.array([removeUnits(w) for w in weights], dtype=float)
    # Convert x to a matrix
    xM = np.reshape(x, (len(x), 1))
    # Compute the differences between all values
    d = np.abs(xM - xM.transpose())
    # Compute the average between all values
    s = np.abs(xM + xM.transpose())
    # Compute the relative difference for
    # the entries which have a non-zero average
    delta = d[np.nonzero(s)]/s[np.nonzero(s)]
    if len(delta) > 0:
        return delta.max()
    else:
        return 0.


def cGtr(weightA, weightB):
    """
    Define the auxiliary greater function.

    Return a number between 0 and 1 depending on how much it is violated
    (0 = A > B, 1 = A << B).

    :returns: XSectioList object with the values for each label.

    """

    if any(not isinstance(w, (float, unum.Unum)) for w in [weightA, weightB]):
        logger.error("Trying to evaluate with invalid objects")
        raise SModelSError()

    # Evaluate the inequality for each cross section info
    a = removeUnits(weightA)
    b = removeUnits(weightB)
    if a + b == 0.:
        return None
    else:
        result = (abs(a - b) - (a - b)) / (2.*(a + b))
    return result


def removeUnits(value, stdUnits=standardUnits, returnUnit=False):
    """
    Remove units from unum objects. Uses the units defined
    in physicsUnits.standard units to normalize the data.

    :param value: Object containing units (e.g. [[100*GeV,100.*GeV],3.*pb])
    :param standardUnits: Unum unit or Array of unum units defined to
                          normalize the data.
    :param returnUnit: If True also resturns the unit corresponding to the returned
                       value.

    :return: Object normalized to standard units (e.g. [[100,100],3000]).
             If returnUnit = True, returns a tuple with the value and its unit (e.g. 100,GeV).
             For unitless values return 1.0 as the unit.
    """

    if isinstance(stdUnits, unum.Unum):
        stdunits = [stdUnits]
    else:
        stdunits = stdUnits

    if isinstance(value, list):
        return [removeUnits(x, stdunits) for x in value]
    if isinstance(value, tuple):
        return tuple([removeUnits(x, stdunits) for x in value])
    elif isinstance(value, dict):
        return dict([[removeUnits(x, stdunits), removeUnits(y, stdunits)]
                     for x, y in value.items()])
    elif isinstance(value, unum.Unum):
        # Check if value has unit or not:
        if not value._unit:
            return value.asNumber()
        # Now try to normalize it by one of the standard pre-defined units:
        for unit in stdunits:
            y = (value/unit).normalize()
            if not y._unit:
                val = value.asNumber(unit)
                if returnUnit:
                    return val, unit
                else:
                    return val
        raise SModelSError("Could not normalize unit value %s using the standard units: %s"
                           % (str(value), str(standardUnits)))
    else:
        if returnUnit:
            return value, 1.0
        else:
            return value


def addUnit(obj, unit):
    """
    Add unit to object.
    If the object is a nested list, adds the unit to all of its elements.

    :param obj: Object without units (e.g. [[100,100.]])
    :param unit: Unit to be added to the object (Unum object, e.g. GeV)
    :return: Object with units (e.g. [[100*GeV,100*GeV]])
    """

    if isinstance(obj, list):
        return [addUnit(x, unit) for x in obj]
    elif isinstance(obj, tuple):
        return tuple([addUnit(x, unit) for x in obj])
    elif isinstance(obj, dict):
        return dict([[addUnit(x, unit), addUnit(y, unit)]
                     for x, y in obj.items()])
    elif isinstance(obj, (float, int, unum.Unum)):
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


def reshapeList(objList, shapeArray):
    """
    Reshape a flat list according to the shape of shapeArray.
    The number of elements in shapeArray should equal the length of objList.

    :param objList: 1D list of objects (e.g. [200,100,'*',50,'*'])
    :param shapeArray: Nested array (e.g. [[float,float,'*',float],'*'])

    :return: Array with elements from objList shaped according to shapeArray
             (e.g. [[200.,100.,'*',50],'*'])
    """

    if isinstance(shapeArray, list):
        return [reshapeList(objList, v) for v in shapeArray]
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
        if inlist[mid] < el:
            lo = mid+1
        else:
            hi = mid
    return lo


def elementsInStr(instring):
    """
    Parse instring and return a list of elements appearing in instring.
    instring can also be a list of strings.

    :param instring: string containing elements (e.g. "[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]")

    :returns: list of elements appearing in instring in string format

    """
    outstr = ""
    if isinstance(instring, str):
        outstr = instring
    elif isinstance(instring, list):
        for st in instring:
            if not isinstance(st, str):
                logger.error("Input must be a string or a list of strings")
                raise SModelSError()
            # Combine list of strings in a single string
            outstr += st
    else:
        raise SModelSError("syntax error in constraint/condition: ``%s''."
              "Check your constraints and conditions in your database." % str(instring))

    outstr = outstr.replace(" ", "")

    if 'PV' in outstr and '>' in outstr:
        if '{' not in outstr or '}' not in outstr:
            raise SModelSError("Elements in %s should be enclosed by curly brackets" % outstr)
        elements = re.findall(r"\{(.*?)\}", outstr)
        return elements

    elements = []
    elStr = ""
    bCounter = 0
    # Parse the string and looks for matching ['s and ]'s, when the matching is
    # complete, store element
    for c in outstr:
        if c == '[':
            bCounter += 1
        elif c == ']':
            bCounter -= 1

        if bCounter != 0:
            elStr += c
        elif elStr:
            elements.append(elStr + c)
            elStr = ""

    # Check if there are not unmatched ['s and/or ]'s in the string
    if bCounter != 0:
        raise SModelSError("Wrong input (incomplete elements?) " + instring)

    # Make sure single quotes are used for the particle strings:
    newElements = []
    for el in elements:
        # Remove all single quotes
        el = el.replace("'", "")
        # Get particles without quotes
        ptcList = el.replace(']', '').replace('[', '').split(',')
        ptcList = [ptc for ptc in ptcList if ptc]
        newEl = el[:]
        iptc = 0
        # Replace particle strings by their values with single quotes
        for ptc in ptcList:
            # Search for ptc only starting after the last replacement
            iptc += newEl[iptc:].find(ptc)
            newEl = newEl[:iptc] + newEl[iptc:].replace(ptc, "'%s'" % ptc, 1)
            iptc += len(ptc) + 1  # Update the position
        newElements.append(newEl)

    return newElements


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

    for attr, value in objDict:
        if attribute == attr:
            values += [value]
        elif isinstance(value, Iterable):
            values += [getValuesForObj(v, attribute) for v in value if v is not obj]
        elif value is not obj:
            values += getValuesForObj(value, attribute)

    values = list(filter(lambda a: (not isinstance(a, list)) or a != [], values))
    values = flattenArray(values)
    uniqueValues = [v for n, v in enumerate(values) if v not in values[:n]]

    return uniqueValues


def bracketToProcessStr(stringSMS, finalState=None, intermediateState=None):
    """
    :parameter stringSMS: string describing the SMS in bracket notation
                         (e.g. [[[e+],[jet]],[[e-],[jet]]])

    :parameter finalState: list containing the final state labels for each branch
                         (e.g. ['MET', 'HSCP'] or ['MET','MET'])
    :parameter intermediateState: nested list containing intermediate state labels
                                     for each branch  (e.g. [['gluino'], ['gluino']])
    """

    # Make sure inclusive symbols come with single quotes:
    branches = eval(stringSMS.replace("*", "'*'").replace("''", "'"))
    # Make replacements to take care of inclusive objects:
    newBranches = []
    for br in branches:
        if br == ['*']:
            newBranches.append('InclusiveNode')
        else:
            newBranches.append([[ptc.replace('*', 'anySM')
                                 for ptc in vt] for vt in br])

    branches = newBranches
    # Get ordered list of BSM states:
    bsmStates = []
    bsmIndices = []
    iptc = 1
    for ibr, br in enumerate(branches):
        bsmStates.append([])
        bsmIndices.append([])
        if br == 'InclusiveNode':
            bsmStates[ibr].append(br)
            bsmIndices[ibr].append(iptc)
            iptc += 1
        else:
            for idec, dec in enumerate(br):
                if intermediateState is None:
                    bsmStates[ibr].append('anyBSM')
                else:
                    bsmStates[ibr].append(intermediateState[ibr][idec])
                bsmIndices[ibr].append(iptc)
                iptc += 1
        if finalState is None:
            bsmStates[ibr].append('MET')
        else:
            bsmStates[ibr].append(finalState[ibr])

    # Get PV:
    pvStrs = []
    for ibr, br in enumerate(branches):
        if not bsmIndices[ibr]:
            pvStrs.append(bsmStates[ibr][0])
        else:
            pvStrs.append(bsmStates[ibr][0]+'(%i)' % bsmIndices[ibr][0])

    pv = '(PV > '+','.join(pvStrs)+')'
    edges = [pv]
    for ibr, br in enumerate(branches):
        # If there is a single BSM state do nothing
        if len(bsmStates[ibr]) < 2:
            continue

        for ibsm, bsm in enumerate(bsmStates[ibr][:-1]):
            iptc = bsmIndices[ibr][ibsm]
            mom = bsm+'(%i)' % iptc
            if bsm == 'InclusiveNode':
                smDaughters = 'anySM'
            else:
                smDaughters = ','.join(br[ibsm])
            bsmDaughter = bsmStates[ibr][ibsm+1]
            if ibsm+1 < len(bsmIndices[ibr]):
                bsmDaughter += '(%i)' % (bsmIndices[ibr][ibsm+1])

            daughters = smDaughters+','+bsmDaughter
            edges.append('('+mom+' > '+daughters+')')

    procStr = ','.join(edges)
    return procStr


def getAttributesFrom(obj, skipIDs=[]):
    """
    Loops over all attributes in the object and return a list
    of the attributes.

    :param obj: Any object with a __dict__ attribute
    :param skipIDs: List of object ids. Any object which has its id on the list
                    will be ignored (useful to avoid recursion).

    :return: List with unique attribute labels.
    """

    if id(obj) in skipIDs or isinstance(obj, (tuple, list, float, int, unum.Unum)):
        return []
    else:
        newSkipIDs = skipIDs[:] + [id(obj)]
    attributes = []
    try:
        objDict = obj.__dict__.items()
    except AttributeError:
        return attributes

    for attr, value in objDict:
        attributes.append(attr)
        if isinstance(value, list):
            attributes += [getAttributesFrom(v, newSkipIDs) for v in value]
        elif isinstance(value, dict):
            attributes += [getAttributesFrom(v, newSkipIDs) for v in value.values()]
        else:
            attributes += getAttributesFrom(value, newSkipIDs)

    attributes = list(filter(lambda a: a != [], attributes))

    return list(set(flattenArray(attributes)))


def rescaleWidth(width):
    """
    The function that is applied to all widths to
    map it into a better variable for interpolation.
    It grows logarithmically from zero (for width=0.)
    to a large number (machine dependent) for width = infinity.

    :param width: Width value (in GeV) with or without units

    :return x: Coordinate value (float)
    """

    if isinstance(width, unum.Unum):
        w = width.asNumber(GeV)
    else:
        w = width

    minWidth = 1e-30  # Any width below this can be safely considered to be zero
    maxWidth = 1e50  # Any width above this can be safely considered to be infinity
    w = (min(w, maxWidth)/minWidth)  # Normalize the width and convert it to some finite number (if not finite)
    if w < 1e-10:  # The log function misbehaves for very small values of w (step behavior), so we use log(1+x) = x for x << 1
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

    minWidth = 1e-30  # Any width below this can be safely considered to be zero
    maxWidth = 1e50  # Any width above this can be safely considered to be infinity
    with np.errstate(over='ignore'):  # Temporarily disable overflow error message
        # The small increase in x is required to
        # enforce unscaleWidth(widthToCoordinae(np.inf)) = np.inf
        width = minWidth*(np.exp(x)-1)
        if width > maxWidth:
            width = np.inf
    return width*GeV


def removeInclusives(massArray, shapeArray):
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
    elif isinstance(massArray, list):
        if len(shapeArray) != len(massArray):
            raise SModelSError("Input value and data shape mismatch (%s,%s)"
                               % (len(shapeArray), len(massArray)))
        else:
            return [removeInclusives(xi, shapeArray[i]) for i, xi in enumerate(massArray)
                    if not removeInclusives(xi, shapeArray[i]) is None]
    else:
        return massArray


def addInclusives(massList, shapeArray):
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

    if isinstance(shapeArray, list):
        return [addInclusives(massList, v) for v in shapeArray]
    elif shapeArray != '*':
        return massList.pop(0)
    else:
        return shapeArray


def sortParticleList(ptcList):
    """
    sorts a list of particle or particle list objects by their label
    :param ptcList: list to be sorted containing particle or particle list objects
    :return: sorted list of particles
    """

    newPtcList = sorted(ptcList, key=lambda x: x.label)
    return newPtcList


def maximal_matching(left, right, edges):
    """
    Computes the maximal matching from left nodes to right nodes.
    The maximal matching is the maximal number of left nodes which can be
    connected to the right nodes without any node belonging to more than one edge.
    Adpated from networkx.algorithms.bipartite.matching.hopcroft_karp_matching.

    :param left: List of left nodes
    :param right: List of right nodes
    :param edges: Nested dictionary with left nodes as keys and macthing right nodes as values
                  (e.g. {nL1 : {nR2 : {}, nR3 : {}}, nL2 : {nR2 : {}, nR1 : {}},... })
    """

    INFINITY = float("inf")
    # Initialize the "global" variables that maintain state during the search.
    leftmatches = {v: None for v in left}
    rightmatches = {v: None for v in right}
    distances = {}
    queue = deque()

    def breadth_first_search():
        for v in left:
            if leftmatches[v] is None:
                distances[v] = 0
                queue.append(v)
            else:
                distances[v] = INFINITY
        distances[None] = INFINITY
        while queue:
            v = queue.popleft()
            if distances[v] < distances[None]:
                for u in edges[v]:
                    if distances[rightmatches[u]] is INFINITY:
                        distances[rightmatches[u]] = distances[v] + 1
                        queue.append(rightmatches[u])
        return distances[None] is not INFINITY

    def depth_first_search(v):
        if v is not None:
            for u in edges[v]:
                if distances[rightmatches[u]] == distances[v] + 1:
                    if depth_first_search(rightmatches[u]):
                        rightmatches[u] = v
                        leftmatches[v] = u
                        return True
            distances[v] = INFINITY
            return False
        return True

    # Implementation note: this counter is incremented as pairs are matched but
    # it is currently not used elsewhere in the computation.
    num_matched_pairs = 0
    while breadth_first_search():
        for v in left:
            if leftmatches[v] is None:
                if depth_first_search(v):
                    num_matched_pairs += 1

    # Strip the entries matched to `None`.
    leftmatches = {k: v for k, v in leftmatches.items() if v is not None}

    return leftmatches
