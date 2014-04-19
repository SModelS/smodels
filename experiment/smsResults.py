"""
.. module:: smsResults
   :synopsis: Centralized facility to access the SMS results.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

from . import smsHelpers
from tools.physicsUnits import addunit, rmvunit
import logging
from experiment.experimentExceptions import MetaInfoError

logger = logging.getLogger(__name__)

allresults = {}
constraints = {}
conditions = {}


def getSmartUpperLimit(analysis, topo, masses):
    """
    Return the upper limit for analysis-topology pair, given an ordered
    sequence of the mass (mother, intermediate, LSP).
    
    """
    import smsInterpolation
    return smsInterpolation.upperLimit(analysis, topo, masses)


def getAllResults(run=None):
    """
    Return all analyses and the topologies they have results for.
    
    """
    import os
    key = str(run)
    if key in allresults:
        return allresults[key]
    runs = smsHelpers.runs
    if run: runs = [ run ]
    ret = {}
    for r in runs:
        dirs = os.listdir ("%s/%s/" % (smsHelpers.base, r))
        for ana in dirs:
            if os.path.exists ("%s/%s/%s/info.txt" % (smsHelpers.base, r, ana)):
                topos = getTopologies (ana, run )
                ret[ana] = topos
    allresults[key] = ret
    return ret


def getTopologies(analysis, run=None):
    """
    Return all topologies that this analysis has constraints for.
    
    """
    run = smsHelpers.getRun(analysis, run)
    x = getConstraints (analysis, run=run)
    return x.keys()


def getConstraints(analysis, topology="all", run=None):
    """
    Get constraints of an analysis. If topology is "all", return a dictionary,
    else return the constraint only for the given topology, None if
    non-existent.
    
    """
    key = analysis + topology + str(run)
    if key in constraints: return constraints[key]
    run = smsHelpers.getRun(analysis, run)
    ret = smsHelpers.constraints(analysis, run)
    if topology == "all":
        constraints[key] = ret
        return ret
    if not topology in ret:
        constraints[key] = None
        return None
    constraints[key] = ret[topology]
    return ret[topology]


def getSqrts(analysis, run=None):
    """ TODO: what is s_hat?
    Get s_hat for this analysis.
    
    """
    sqrts = smsHelpers.getMetaInfoField(analysis, "sqrts", run)
    try:
        return addunit(float(sqrts), "TeV")
    except: # TODO: except what?
        pass
    return sqrts


def getConditions (analysis, topology="all", fuzzy=True, run=None):
    """
    Get conditions of an analysis. if topology is "all", return a dictionary,
    else return the condition only for the given topology, None if
    non-existent.
    
    """
    key = analysis + topology + str(fuzzy) + str(run)
    if key in conditions: return conditions[key]
    run = smsHelpers.getRun (analysis, run)
    if fuzzy:
        ret = smsHelpers.getLines(analysis, run, "fuzzycondition")
    else:
        ret = smsHelpers.getLines(analysis, run, "condition")
    if topology == "all":
        conditions[key] = ret
        return ret
    if not topology in ret:
        conditions[key] = None
        return None
    conditions[key] = ret[topology]
    return ret[topology]


def getaxes(analysis, topology=None, run=None):
    """ TODO: improve docstring
    Get information about the histogram axes for an analysis.
    
    For each topology list of dictionary, each dictionary corresponds to one
    histogram. The key axes gives string (mx-my), the key mz gives information
    on other masses, if you supply a topology, returns list for this topology
    only.
    
    """
    if not exists(analysis, topology=None):    # analysis is not known, we return None
        return None
    try:
        st = smsHelpers.getMetaInfoField (analysis, "axes", run)
    except MetaInfoError:
        logger.error("Meta info field 'axes' does not exist in %s." %analysis)
        st = None
    if not st:
        if not topology: return None    # cannot return default without info on topology
# if there is no information about the axes, return the default
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
            else: d[nm].append({'axes': n.pop(0) + '-' + n.pop(0), 'mz': n})

    if topology:
        topology = topology.replace(" ", "")
        if not d or not topology in d:
    # # topology does not exist, we return None
            return None
        else: return d[topology]

    return d


def exists(analysis, topology, run=None):
    """
    Check if the dictioanary 'limit_topo' in run/analysis/sms.py exists.
    
    For topologies with intermediate masses, check if all dictionaries listed
    in the axes-information exist. If topology==None, simply check if
    run/analysis/sms.py exists.
    
    """
    run2 = smsHelpers.getRun(analysis, run)
    if not topology:
        import os
        base = smsHelpers.base
        pydict = "%s/%s/%s/sms.py" % (base, run2, analysis)
        if os.path.exists(pydict):
            return True
        else:
            return False
    axes = getaxes(analysis, topology)
    if not axes:
        return False
    hasDict = smsHelpers.hasDictionary (analysis, run2)
    if not hasDict: return False
    return smsHelpers.hasDictionary (analysis, run2, topology )


def getUpperLimit (analysis, topo, mx=None, my=None, run=None, interpolate=False, expected=False):
    """
    Get the upper limit for run/analysis/topo.
    
    :returns: None, if it does not exist. If mx and my are None, return the
    entire dictionary, if mx and my are floats, return the upper limit at this
    point.
    
    """
    run = smsHelpers.getRun(analysis, run)
    if smsHelpers.hasDictionary(analysis, run):
        return getUpperLimitFromDictionary (analysis, topo, mx, my, run, interpolate=interpolate, expected=expected)
    logger.warning ("no upper limits found for %s" %analysis)
    return None


def getUpperLimitFromDictionary (analysis, topo, mx=None, my=None, run=None, png=None, interpolate=False, expected=False):
    """
    Get an upper limit from the python dictionary.
    
    """
    Dict = smsHelpers.getUpperLimitDictionary(analysis, topo, run, expected=expected)
    if Dict == None:
        return Dict
    if mx == None and my == None:
        return Dict
    if mx == None or my == None:
        logger.error("Requesting upper limits for mx = %s and my = %s" % (mx, my))
        return None
    if rmvunit(mx, 'GeV') == None: return Dict
    return addunit(getInterpolatedUpperLimitDelaunay(Dict, mx, my), "pb")


def getInterpolatedUpperLimitDelaunay(Dict, inmx, inmy):
    """
    Get interpolated upper limit from dictionary at point (inmx, inmy).
    
    :param Dict: dictionary (sms.py), contains upper limits of one analysis and
    one topology
    :param inmx: mass point on x-axis
    :param inmy: mass point on y-axis
    :returns: interpolated upper limit at point (inmx, inmy)
    
    """
    try:
        import numpy as np
        import scipy.interpolate as ip
        mx = rmvunit(inmx, 'GeV')
        my = rmvunit(inmy, 'GeV')
        if not inConvexHull(Dict, mx, my):
            logger.debug("Cannot interpolate for (%f, %f), point is not in \
                          convex hull." % (inmx, inmy))
            return None
        n = 0
        for k in Dict:
            n += len(Dict[k])
        points = np.zeros((n, 2))
        values = np.zeros((n))
        i = 0
        for x in Dict:
            for y in Dict[x]:
                points[i] = [x, y]
                values[i] = Dict[x][y]
                i += 1
        grid_x = np.zeros((1, 1))
        grid_y = np.zeros((1, 1))
        grid_x = mx
        grid_y = my
        return float(ip.griddata(points, values, (grid_x, grid_y), method='linear'))
    except Exception as e: # TODO: which exception?
        logger.error("Cannot interpolate %s. Using closest value instead." % str(e))
        if not inConvexHull(Dict, inmx, inmy):
            return False
        return getClosestValue(Dict, inmx, inmy)
    

def inConvexHull(Dict, mx, my):
    pointlist = []
    for k in Dict.keys():
        for ki in Dict[k].keys():
            pointlist.append([k, ki])
#    try:
    import numpy
    p = numpy.array(pointlist)
    from scipy.spatial import Delaunay
    dela = Delaunay(p)
    return dela.find_simplex((mx, my)) >= 0


def getClosestValue (Dict, mx, my):
    """
    Get the upper limit of the point in Dict that is closest to mx and my,
    assuming that Dict is a dictionary of mx,my,ul.
    
    """
    closest = 9999999
    retul = None
    for (dmx, dmv) in Dict.items():
        for (dmy, ul) in dmv.items():
            dist = (mx - dmx) ** 2 + (my - dmy) ** 2
            if dist < closest:
                closest = dist
                retul = ul
    if closest > 20.**2:    # # if we're more than 20 gev from the closest point, we return False
        return False
    return retul

