"""
.. module:: expAuxiliaryFunctions
   :synopsis: A collection of functions used to evaluate fuzzy the conditions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.base.physicsUnits import standardUnits, GeV
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.base.smodelsLogging import logger

import unum
import re
import numpy as np
from collections import OrderedDict, deque

try:
    from collections.abc import Iterable
except (ImportError, ModuleNotFoundError):
    from collections import Iterable

minWidth = 1e-30  # Any width below this can be safely considered to be zero
maxWidth = 1e50  # Any width above this can be safely considered to be infinity


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
    in base.physicsUnits.standard units to normalize the data.

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


def smsInStr(instring):
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
            raise SModelSError(f"Elements in {outstr} should be enclosed by curly brackets")
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
            newEl = newEl[:iptc] + newEl[iptc:].replace(ptc, f"'{ptc}'", 1)
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


def bracketToProcessStr(stringSMS, 
                        finalState=None, intermediateState=None,
                        returnNodeDict=False):
    """
    :parameter stringSMS: string describing the SMS in bracket notation
                         (e.g. [[[e+],[jet]],[[e-],[jet]]])

    :parameter finalState: list containing the final state labels for each branch
                         (e.g. ['MET', 'HSCP'] or ['MET','MET'])
    :parameter intermediateState: nested list containing intermediate state labels
                                     for each branch  (e.g. [['gluino'], ['gluino']])
    :param returnNodeDict: If True, return a dictionary mapping the nested bracket
                           indices to the particle nodes 
                           ({(branchIndex,vertexIndex) : nodeIndex})
                                     
    :return: process string in new format (str) and dictionary nodes dictionary (if returnNodeDict=True) 
    """
    
    # Make sure all particle labels are enclosed by single quotes:
    stringSMS = stringSMS.replace("'","").replace(" ","")
    # Add quotes
    addQuotes = lambda match: "[%s]" %(",".join(["'%s'" %ptc for ptc in match.group(1).split(',')]))
    # Find all strings surrounded by brackets and replace
    # them by the string with quotes added
    stringSMS = re.sub(r"\[([^\[\]]+)\]", addQuotes, stringSMS)

    # Convert string to nested list:
    smsList = eval(stringSMS)

    # Collect all BSM states in each branch
    # and convert inclusive branches
    branches = []
    bsmStates = []
    for ibr,br in enumerate(smsList):
        bsmStates.append([])
        if br == ['*']:
            branches.append([['*anySM']])
            bsmStates[ibr] = ['InclusiveNode']
        else:
            branches.append([sorted([ptc if ptc != '*' else  'anySM'
                                 for ptc in vt]) for vt in br])
            if intermediateState is None:
                bsmStates[ibr] = ['anyBSM']*len(br)
            else:
                bsmStates[ibr] = intermediateState[ibr][:]

        if finalState is None:
            bsmStates[ibr].append('MET')
        else:
            bsmStates[ibr].append(finalState[ibr])

    # Create decay string for primary vertex:
    bsmNodesDict = {}
    decayStrings = []
    daughtersStr = []
    momStr = 'PV(0)'
    nodeIndex = 1
    for ibr,bsmBranch in enumerate(bsmStates):
        bsmNodesDict[(ibr,0)] = nodeIndex
        nodeIndex += 1
        daughtersStr.append('%s(%i)' %(bsmBranch[0],bsmNodesDict[(ibr,0)]))
    decStr = f"({momStr} > {','.join(daughtersStr)})"
    decayStrings.append(decStr)

    # Create decay string for all branches:
    smNodesDict = {}
    for ibr,bsmBranch in enumerate(bsmStates):
        for idec,bsm in enumerate(bsmBranch[:-1]):        
            # Define node index for BSM mom
            if (ibr,idec) not in bsmNodesDict:
                bsmNodesDict[(ibr,idec)] = nodeIndex
                nodeIndex += 1
            # Create string with node index for BSM mom
            momStr = '%s(%i)' %(bsm,bsmNodesDict[(ibr,idec)])
            daughtersStr = []
            # Define node index for BSM daughter 
            if (ibr,idec+1) not in bsmNodesDict:
                bsmNodesDict[(ibr,idec+1)] = nodeIndex
                nodeIndex += 1
            # Create string with node index for mom            
            daughtersStr.append('%s(%i)' %(bsmBranch[idec+1],bsmNodesDict[(ibr,idec+1)]))
            # Define node indices for SM daughters:        
            for ism,_ in enumerate(branches[ibr][idec]):
                smNodesDict[(ibr,idec,ism)] = nodeIndex
                nodeIndex += 1
            # Create string with node index for mom            
            daughtersStr += ['%s(%i)' %(sm,smNodesDict[(ibr,idec,ism)]) 
                             for ism,sm in enumerate(branches[ibr][idec])]
            # Create decay string:
            decStr = f"({momStr} > {','.join(daughtersStr)})"
            decayStrings.append(decStr)

    processStr = ', '.join(decayStrings)
    
    if returnNodeDict:
        return processStr,bsmNodesDict
    else:
        return processStr


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

    w = (min(w, maxWidth)/minWidth)  # Normalize the width and convert it to some finite number (if not finite)
    if w < 1e-10:  # The log function misbehaves for very small values of w (step behavior), so we use log(1+x) = x for x << 1
        return w
    else:
        return np.log(1+w)


def unscaleWidth(x):
    """
    Maps a coordinate value back to width.
    The mapping is such that x=0->width=0 and x=very large -> width = inf.

    :param x: Coordinate value (float)

    :return width: Width value without units
    """

    with np.errstate(over='ignore'):  # Temporarily disable overflow error message
        # The small increase in x is required to
        # enforce unscaleWidth(widthToCoordinae(np.inf)) = np.inf
        width = minWidth*(np.exp(x)-1)
        if width > maxWidth:
            width = np.inf
    return width


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


def cleanWalk ( topdir ):
    """ perform os.walk, but ignore all hidden files and directories """
    import os
    ret = []
    for root, d_, f_ in os.walk ( topdir ):
        isHidden=False
        tokens = root.split("/")
        for token in tokens:
            if len(token)>0 and token[0]==".":
                isHidden=True
                break
        if isHidden:
            continue
        dirs,files = [],[]
        for d in d_:
            if not d[0]==".":
                dirs.append ( d )
        for f in f_:
            if not f[0]==".":
                files.append ( f )
        ret.append ( [ root, dirs, files ] )
    return ret

def concatenateLines ( oldcontent ):
    """ of all lines in the list "oldcontent", concatenate the ones
        that end with '\'  or ',' """
    content=[] ## concatenate lines that end with "," or "\"
    tmp=""
    import re
    for i,line in enumerate ( oldcontent ):
        p1 = line.find ( "#" )
        if p1 > -1:
            line = line[:p1]
        tmp+=line.strip()
        ## if next line starts with tab or whitespace or "}",
        ## merge the lines
        if i < len(oldcontent)-1 and re.match("[ \t}]",oldcontent[i+1] ):
            # if next line starts with white space, we add also
            continue
        if tmp != "" and tmp[-1] not in [ ",", '\\' ]:
            content.append ( tmp )
            tmp=""
        if tmp != "" and tmp[-1] == '\\':
            tmp=tmp[:-1] # remove trailing \ (but keep trailing ,)
    return content
