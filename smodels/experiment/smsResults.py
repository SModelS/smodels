"""
.. module:: experiment.smsResults
   :synopsis: Centralized facility to access the SMS results.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

from smodels.tools.physicsUnits import TeV, pb, fb, GeV
## from smodels.tools import rcFile
from smodels.experiment import smsHelpers
from smodels.experiment.experimentExceptions import MetaInfoError
from smodels.experiment.smsHelpers import getPath
from smodels.tools.uniqueLogFilter import UniqueFilter
import logging

logger = logging.getLogger(__name__)
logger.addFilter(UniqueFilter())

allresults = {}
constraints = {}
conditions = {}


def getAllResults(sqrts=None, experiment=None):
    """
    Get all analyses and topologies that have results.

    """
    import os
    key = str(sqrts)+str(experiment)
    if key in allresults:
        return allresults[key]
    sqrtsList = smsHelpers.sqrts
    experimentList = smsHelpers.experiments
    if sqrts:
        sqrtsList = [sqrts]
    if experiment:
        experimentList = [experiment]
    ret = {}
    for s in sqrtsList:
        for e in experimentList:
            if not os.path.exists("%s/%s/%s/" % (smsHelpers.base, s, e)):
                logger.warning("Expected directory %s/%s was not found in the database at %s" % (s, e, smsHelpers.base))
                continue
            dirs = os.listdir("%s/%s/%s" % (smsHelpers.base, s, e))
            for ana in dirs:
                if os.path.exists("%s/%s/%s/%s/info.txt" % (smsHelpers.base, s, e, ana)):
                    topos = getTopologies(ana)
                    ret[ana] = topos
    allresults[key] = ret
    return ret


def getTopologies(analysis):
    """
    Get all topologies of an analysis with constraints.

    """
    path = smsHelpers.getPath(analysis)
    x = getConstraints(analysis, path=path)
    return x.keys()


def getConstraints(analysis, topology="all", path=None):
    """
    Get constraints of an analysis.

    :returns: dictionary of constraints, if topology == "all"; single
              constraint for the passed topology, if only one topology is passed; 
              None if non-existent;

    """
    key = analysis + topology + str(path)
    if key in constraints:
        return constraints[key]
    path = smsHelpers.getPath(analysis, path)
    ret = smsHelpers.getLines(analysis, path, "constraint")
    if topology == "all":
        constraints[key] = ret
        return ret
    if not topology in ret:
        constraints[key] = None
        return None
    constraints[key] = ret[topology]
    return ret[topology]


def getBranchCondition(anaName, txName, path=None):
    """
    Get the branch condition for an analysis.

    :returns: string containing the branch condition (e.g. equal branch masses)
    """
        
    path = smsHelpers.getPath(anaName, path)
    ret = smsHelpers.getLines(anaName, path, "branchcondition")
    if not txName in ret: return None
    else: return ret[txName]


def getSqrts(analysis, path=None):
    """ get the center-of-mass energy of the analysis.
    """
    sqrts = smsHelpers.getMetaInfoField(analysis, "sqrts", path)
    try:
        return float(sqrts) * TeV

    except ValueError:
        try:
            return eval(sqrts)
        except:
            pass
    return sqrts


def getConditions(analysis, topology="all", fuzzy=True, path=None):
    """
    Get conditions of an analysis.

    :returns: dictionary of conditions, if topology == "all"; single condition
              for the passed topology, if only one topology is passed; None if
              non-existent.

    """
    key = analysis + topology + str(fuzzy) + str(path)
    if key in conditions:
        return conditions[key]
    path = smsHelpers.getPath(analysis, path)
    if fuzzy:
        ret = smsHelpers.getLines(analysis, path, "fuzzycondition")
    else:
        ret = smsHelpers.getLines(analysis, path, "condition")
    if topology == "all":
        conditions[key] = ret
        return ret
    if not topology in ret:
        conditions[key] = None
        return None
    conditions[key] = ret[topology]
    return ret[topology]


def getaxes(analysis, topology=None, path=None):
    """Get information about the histogram axes for an analysis.

    For each topology list of dictionary, each dictionary corresponds to one
    histogram. The key axes gives string (mx-my), the key mz gives information
    on other masses, if you supply a topology, returns list for this topology
    only.

    """
    if not _exists(analysis, topology=None):
        return None
    try:
        st = smsHelpers.getMetaInfoField(analysis, "axes", path)
    except MetaInfoError:
        logger.error("Meta info field 'axes' does not exist in %s.", analysis)
        st = None
    if not st:
        if not topology:
            # Cannot return default without info on topology
            return None
        # If there is no information about the axes, return the default
        return [{'axes': 'M1-M0', 'mz': None}]
    st = st.split(',')
    d = {}
    for i in range(len(st)):
        l = st[i].split(':')
        nm = l[0].replace(" ", "")
        d[nm] = []
        m = l[1].split('-')
        for j in range(len(m)):
            n = m[j].split()
            if len(n) == 2:
                d[nm].append({'axes': n[0] + '-' + n[1], 'mz': None})
            else:
                d[nm].append({'axes': n.pop(0) + '-' + n.pop(0), 'mz': n})

    if topology:
        topology = topology.replace(" ", "")
        if not d or not topology in d:
            # Return None, if topology does not exist
            return None
        else:
            return d[topology]

    return d

def setBase(base):
    """
    Set the base directory of the database.
    
    """
    smsHelpers.base = base

def getBase():
    """
    Return the base directory of the database.
    
    """
    return smsHelpers.base

def getURL(analysis, path=None):
    """
    Get the URL of an analysis.
    
    """
    logger.warning ("getURL is deprecated")
    return smsHelpers.getMetaInfoField(analysis, "url", path)


def hasURL(analysis, path=None):
    """
    Check if URL of an analysis exists."""
    logger.warning("hasURL is deprecated")
    return smsHelpers.getMetaInfoField(analysis, "url", path)


def getPAS(analysis, path=None):
    """
    Get the PAS of an analysis.
    
    """
    logger.warning("getPAS is deprecated")
    return smsHelpers.getMetaInfoField(analysis, "pas", path)


def getJournal(analysis, path=None):
    """
    Get the journal of an analysis.
    
    """
    logger.warning ("getJournal is deprecated")
    return smsHelpers.getMetaInfoField(analysis, "journal", path)


def getLumi(analysis, path=None):
    """
    Get the integrated luminosity for an analysis.
    
    """
    lumifb = smsHelpers.getMetaInfoField(analysis, "lumi", path)
    try:
        return float(lumifb) / fb
    except ValueError:
        try:
            return eval(lumifb)
        except:
            pass
    return lumifb


def isPrivate(analysis, path=None):
    """
    Check if analysis is flagged as private.
    """
    field=smsHelpers.getMetaInfoField(analysis, "private", path)
    if field==None: return False
    return bool(int(field))

def isSuperseded (analysis, path=None):
    """
    check if analysis is superseded, if yes,
    return analysis name of newer analysis
    """
    return smsHelpers.getMetaInfoField (analysis, "superseded_by", path)


def getExperiment(analysis, path=None):
    """
    Check if path is to ATLAS directory, else return CMS.
    
    """
    logger.warning("getExperiment is deprecated")
    path1 = getPath(analysis, path)
    if path1.find("ATLAS") > -1:
        return "ATLAS"
    return "CMS"


def getComment(analysis, path=None):
    """
    Get the comment of an analysis.
    
    """
    logger.warning ("getComment is deprecated")
    return smsHelpers.getMetaInfoField(analysis, "comment", path)


def considerExperiment(experiment):
    """
    Define the experiment to be considered.
    
    """
    smsHelpers.experiments = [experiment]

def considerSqrts(sqrts):
    """
    Define the center of mass energies to be considered (as strings in list format)
    """
    smsHelpers.sqrts = [sqrts]


def _exists(analysis, topology, path=None):
    """
    Check if the dictionary 'limit_topo' in <path>/sms.py exists.

    For topologies with intermediate masses, check if all dictionaries listed
    in the axes-information exist. If topology == None, check if
    <path>/sms.py exists.

    """
    path2 = smsHelpers.getPath(analysis, path)
    if not topology:
        import os
        base = smsHelpers.base
        pydict = "%s/%s/%s/sms.py" % (base, path2, analysis)
        if os.path.exists(pydict):
            return True
        else:
            return False
    axes = getaxes(analysis, topology)
    if not axes:
        return False
    hasDict = smsHelpers.hasDictionary(analysis, path2)
    if not hasDict:
        return False
    return smsHelpers.hasDictionary(analysis, path2, topology)


def getUpperLimit(analysis, topology, mx=None, my=None, path=None,
                  interpolate=False, expected=False):
    """
    Get the upper limit for path/analysis/topology.

    :returns: None, if it does not exist; entire dictionary, if mx and my are
              None; upper limit at mx/my, if mx and my are floats;

    """
    path = smsHelpers.getPath(analysis, path)
    if smsHelpers.hasDictionary(analysis, path):
        return getUpperLimitFromDictionary(analysis, topology, mx, my, path,
                                           interpolate=interpolate,
                                           expected=expected)
    logger.warning("No upper limits found for %s", analysis)
    return None


def getUpperLimitFromDictionary(analysis, topology, mx=None, my=None,
                                path=None, png=None, interpolate=False,
                                expected=False):
    """
    Get an upper limit from the python dictionary.

    """
    dictionary = smsHelpers.getUpperLimitDictionary(analysis, topology, path,
                                                    expected=expected)
    if dictionary == None:
        return dictionary
    if type(mx) == type(None) and type(my) == type(None):
        return dictionary
    if type(mx) == type(None) or type(my) == type(None):
        logger.error("Requesting upper limits for mx = %s and my = %s", mx, my)
        return None
    if not (type(mx) == type(1.) or type(mx) == type(GeV)):
        return dictionary
    if not getInterpolatedUpperLimitDelaunay(dictionary, mx, my): return None
    return getInterpolatedUpperLimitDelaunay(dictionary, mx, my) * pb


def getInterpolatedUpperLimitDelaunay(dictionary, inmx, inmy):
    """
    Get interpolated upper limit from dictionary at point (inmx, inmy).

    :param dictionary: dictionary (sms.py), contains upper limits of one
                       analysis and one topology
    :param inmx: mass point on x-axis
    :param inmy: mass point on y-axis
    :returns: interpolated upper limit at point (inmx, inmy)

    """
    import numpy as np
    import scipy.interpolate as ip

    try:
        mx = inmx / GeV
        my = inmy / GeV
        if not inConvexHull(dictionary, mx, my):
            logger.debug("Cannot interpolate for (%s, %s), point is not in "
                         "convex hull.", str(inmx), str(inmy))
            return None
        n = 0
        for k in dictionary:
            n += len(dictionary[k])
        points = np.zeros((n, 2))
        values = np.zeros((n))
        i = 0
        for x in dictionary:
            for y in dictionary[x]:
                points[i] = [x, y]
                values[i] = dictionary[x][y]
                i += 1
        gridX = np.zeros((1, 1))
        gridY = np.zeros((1, 1))
        gridX = mx
        gridY = my
        return float(ip.griddata(points, values, (gridX, gridY),
                                 method='linear'))
    except Exception as e:
        logger.error("Cannot interpolate %s. Using closest value instead.", e)
        if not inConvexHull(dictionary, inmx, inmy):
            return False
        return getClosestValue(dictionary, inmx, inmy)


def inConvexHull(dictionary, mx, my):
    """
    Check if (mx,my) point is in the data dictionary.

    """
    import numpy
    from scipy.spatial import Delaunay

    pointlist = []
    for k in dictionary.keys():
        for ki in dictionary[k].keys():
            pointlist.append([k, ki])
    p = numpy.array(pointlist)
    dela = Delaunay(p)
    return dela.find_simplex((mx, my)) >= 0


def getClosestValue(dictionary, mx, my):
    """
    Get the upper limit of the point in dictionary that is closest to mx and
    my, assuming that dictionary is a dictionary of mx, my, ul.

    """
    closest = 9999999
    retul = None
    for (dmx, dmv) in dictionary.items():
        for (dmy, ul) in dmv.items():
            dist = (mx - dmx) ** 2 + (my - dmy) ** 2
            if dist < closest:
                closest = dist
                retul = ul
    if closest > 20.**2:
        # Return False, if distance > 20 GeV from closest point
        return False
    return retul

