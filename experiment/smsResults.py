"""
.. module:: experiment.smsResults
   :synopsis: Centralized facility to access the SMS results.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

from smodels.tools.physicsUnits import addunit, rmvunit
from smodels.experiment import smsHelpers
from smodels.experiment.experimentExceptions import MetaInfoError
from smodels.experiment.smsHelpers import getRun
import logging


logger = logging.getLogger(__name__)

allresults = {}
constraints = {}
conditions = {}


def getAllResults(run=None):
    """
    Get all analyses and topologies that have results.

    """
    import os
    key = str(run)
    if key in allresults:
        return allresults[key]
    runs = smsHelpers.runs
    if run:
        runs = [run]
    ret = {}
    for r in runs:
        dirs = os.listdir("%s/%s/" % (smsHelpers.base, r))
        for ana in dirs:
            if os.path.exists("%s/%s/%s/info.txt" % (smsHelpers.base, r, ana)):
                topos = getTopologies(ana, run)
                ret[ana] = topos
    allresults[key] = ret
    return ret


def getTopologies(analysis, run=None):
    """
    Get all topologies of an analysis with constraints.

    """
    run = smsHelpers.getRun(analysis, run)
    x = getConstraints(analysis, run=run)
    return x.keys()


def getConstraints(analysis, topology="all", run=None):
    """
    Get constraints of an analysis.

    :returns: dictionary of constraints, if topology == "all"; single
    constraint for the passed topology, if only one topology is passed; None if
    non-existent;

    """
    key = analysis + topology + str(run)
    if key in constraints:
        return constraints[key]
    run = smsHelpers.getRun(analysis, run)
    ret = smsHelpers.getLines(analysis, run, "constraint")
    if topology == "all":
        constraints[key] = ret
        return ret
    if not topology in ret:
        constraints[key] = None
        return None
    constraints[key] = ret[topology]
    return ret[topology]


def getSqrts(analysis, run=None):
    """ get the center-of-mass energy of the analysis.
    """
    sqrts = smsHelpers.getMetaInfoField(analysis, "sqrts", run)
    try:
        return addunit(float(sqrts), "TeV")
    except:  # TODO: except what?
        pass
    return sqrts


def getConditions(analysis, topology="all", fuzzy=True, run=None):
    """
    Get conditions of an analysis.

    :returns: dictionary of conditions, if topology == "all"; single condition
    for the passed topology, if only one topology is passed; None if
    non-existent

    """
    key = analysis + topology + str(fuzzy) + str(run)
    if key in conditions:
        return conditions[key]
    run = smsHelpers.getRun(analysis, run)
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
    if not _exists(analysis, topology=None):
        return None
    try:
        st = smsHelpers.getMetaInfoField(analysis, "axes", run)
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


def getURL(analysis, run=None):
    """
    Get the URL of an analysis.
    
    """
    logger.warning("getURL is deprecated")
    return smsHelpers.getMetaInfoField(analysis, "url", run)


def hasURL(analysis, run=None):
    """
    Check if URL of an analysis exists."""
    logger.warning("hasURL is deprecated")
    return smsHelpers.getMetaInfoField(analysis, "url", run)


def getPAS(analysis, run=None):
    """
    Get the PAS of an analysis.
    
    """
    logger.warning("getPAS is deprecated")
    return smsHelpers.getMetaInfoField(analysis, "pas", run)


def getJournal(analysis, run=None):
    """
    Get the journal of an analysis.
    
    """
    logger.warning("getJournal is deprecated")
    return smsHelpers.getMetaInfoField(analysis, "journal", run)


def getLumi(analysis, run=None):
    """
    Get the integrated luminosity for an analysis.
    
    """
    logger.warning("getLumi is deprecated")
    lumifb = float(smsHelpers.getMetaInfoField(analysis, "lumi", run))
    return addunit(lumifb, "fb-1")


def getExperiment(analysis, run=None):
    """
    Check if run is ATLAS8TeV, else return CMS.
    
    """
    logger.warning("getExperiment is deprecated")
    run1 = getRun(analysis, run)
    if run1.find("ATLAS") > -1:
        return "ATLAS"
    return "CMS"


def getComment(analysis, run=None):
    """
    Get the comment of an analysis.
    
    """
    logger.warning("getComment is deprecated")
    return smsHelpers.getMetaInfoField(analysis, "comment", run)


def considerRuns(runs):
    """
    TODO: write docstring
    
    """
    smsHelpers.runs = runs


def _exists(analysis, topology, run=None):
    """
    Check if the dictionary 'limit_topo' in <run>/<analysis>/sms.py exists.

    For topologies with intermediate masses, check if all dictionaries listed
    in the axes-information exist. If topology == None, check if
    <run>/<analysis>/sms.py exists.

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
    hasDict = smsHelpers.hasDictionary(analysis, run2)
    if not hasDict:
        return False
    return smsHelpers.hasDictionary(analysis, run2, topology)


def getUpperLimit(analysis, topology, mx=None, my=None, run=None,
                  interpolate=False, expected=False):
    """
    Get the upper limit for run/analysis/topology.

    :returns: None, if it does not exist; entire dictionary, if mx and my are
    None; upper limit at mx/my, if mx and my are floats;

    """
    run = smsHelpers.getRun(analysis, run)
    if smsHelpers.hasDictionary(analysis, run):
        return getUpperLimitFromDictionary(analysis, topology, mx, my, run,
                                           interpolate=interpolate,
                                           expected=expected)
    logger.warning("No upper limits found for %s", analysis)
    return None


def getUpperLimitFromDictionary(analysis, topology, mx=None, my=None,
                                run=None, png=None, interpolate=False,
                                expected=False):
    """ TODO: unused arguments png and interpolate
    Get an upper limit from the python dictionary.

    """
    dictionary = smsHelpers.getUpperLimitDictionary(analysis, topology, run,
                                                    expected=expected)
    if dictionary == None:
        return dictionary
    if mx == None and my == None:
        return dictionary
    if mx == None or my == None:
        logger.error("Requesting upper limits for mx = %s and my = %s", mx, my)
        return None
    if rmvunit(mx, 'GeV') == None:
        return dictionary
    return addunit(getInterpolatedUpperLimitDelaunay(dictionary, mx, my), "pb")


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
        mx = rmvunit(inmx, 'GeV')
        my = rmvunit(inmy, 'GeV')
        if not inConvexHull(dictionary, mx, my):
            logger.debug("Cannot interpolate for (%f, %f), point is not in "
                         "convex hull.", inmx, inmy)
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
    except Exception as e:  # TODO: which exception?
        logger.error("Cannot interpolate %s. Using closest value instead.", e)
        if not inConvexHull(dictionary, inmx, inmy):
            return False
        return getClosestValue(dictionary, inmx, inmy)


def inConvexHull(dictionary, mx, my):
    """
    TODO: write docstring

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

