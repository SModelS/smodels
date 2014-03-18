#!/usr/bin/env python

"""
.. module:: SMSResults
   :synopsis: Centralized facility to access the SMS results.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

from experiment import smsHelpers
from tools.PhysicsUnits import addunit, rmvunit
from tools import PhysicsUnits, RCFile
from experiment import logger
from experiment.experimentExceptions import MetaInfoError
from experiment.smsHelpers import databaseVersion, getRun

def setBase (base):
    """ just sets the base directory of the database """
    smsHelpers.Base = base

def useUnits (b=True):
    PhysicsUnits.useUnits = b

alldirectories = ['8TeV', 'ATLAS8TeV', '2012', '2011']

def considerRuns(run=None):
    """ 
        defines what runs are to be considered when asking for results.

        :param run: a list of runs to be considered, e.g. [ '2012', '8TeV' ]). If None, all runs are taken into account.
        :type run: list or NoneType
    """
    allruns = ["8TeV", "ATLAS8TeV", "RPV8", "2012", "RPV7", "2011"]
    runsort = []
    if run:
        for r in allruns:
            if r in run:
                runsort.append(r)
        smsHelpers.runs = runsort
        alldirectories = runsort
        for r in run:
            if not r in allruns:
                print "%s is not a run!!!" % r
    else:
        smsHelpers.runs = allruns
        alldirectories = allruns

def getTopologies (analysis, run=None ):
    """ return all topologies that this analysis has constraints for """
    run = smsHelpers.getRun (analysis, run)
    x = getConstraints (analysis, run=run)
    return x.keys()

def getExperiment (analysis, run=None):
    """ return experiment name for given analysis
            for now: check if run is ATLAS8TeV, else return CMS """
    run1 = getRun(analysis, run)
    if run1 == "ATLAS8TeV": return "ATLAS"
    else: return "CMS"


def getPrettyName (analysis, run=None, latex=False):
    value = smsHelpers.getMetaInfoField (analysis, "prettyname", run)
    if value == None or value == "": return analysis
    if not latex:
        value = value.replace("\\", "#")
    return value

def getAnalyses (topo, run=None):
    """return all analyses that have results for topo
    """
    import os
    runs = smsHelpers.runs
    if run: runs = [ run ]
    analyses = {}
    for r in runs:
# # so thats the runs I really have to think about
        dirs = os.listdir ("%s/%s/" % (smsHelpers.Base, r))
        for ana in dirs:
            if os.path.exists ("%s/%s/%s/info.txt" % (smsHelpers.Base, r, ana)):
#                e=getExclusion ( ana, topo, r )
#                if e: analyses[ana]=True
                if exists(ana, topo, r): analyses[ana] = True

    return analyses.keys()

allresults = {}

def getAllResults (run=None ):
    """returns all analyses and the topologies they have results for
    """
    import os
    key = str(run)
    if allresults.has_key (key):
        return allresults[key]
    runs = smsHelpers.runs
    if run: runs = [ run ]
    ret = {}
    for r in runs:
# # so thats the runs I really have to think about
        dirs = os.listdir ("%s/%s/" % (smsHelpers.Base, r))
        for ana in dirs:
            if os.path.exists ("%s/%s/%s/info.txt" % (smsHelpers.Base, r, ana)):
                topos = getTopologies (ana, run )
                ret[ana] = topos
    allresults[key] = ret
    return ret

def getClosestValue (Dict, mx, my):
    """assuming that Dict is a dictionary of mx,my,ul, get the upper limit
       of the point in Dict that is closest to mx and my.
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

def getInterpolatedUpperLimitDelaunay (Dict, inmx, inmy):
    """ get interpolated upper limit from dictionary at point (inmx, inmy)
            :param Dict: dictionray (sms.py), contains upper limits of one analysis and one topology
            :param inmx: mass point on x-axis
            :param inmy: mass point on y-axis
            :returns: interpolated upper limit at point (inmx, inmy) """
    try:
        import numpy as np
        import scipy.interpolate as ip
        mx = rmvunit(inmx, 'GeV')
        my = rmvunit(inmy, 'GeV')
        if not inConvexHull(Dict, mx, my):
            logger.debug ("Can\'t interpolate for (%f,%f), point is not in convex hull." % (inmx, inmy))
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
    except Exception, e:
        logger.error ("cannot interpolate: %s. use closest value." % str(e))
        if not inConvexHull (Dict, inmx, inmy): return False
        return getClosestValue (Dict, inmx, inmy)

def getInterpolatedUpperLimit (Dict, inmx, inmy):
    """ get interpolated upper limit from dictionary at point (inmx, inmy)
            :param Dict: dictionray (sms.py), contains upper limits of one analysis and one topology
            :param inmx: mass point on x-axis
            :param inmy: mass point on y-axis
            :returns: interpolated upper limit at point (inmx, inmy) """
    import scipy.interpolate as ip
    mx = rmvunit(inmx, 'GeV')
    my = rmvunit(inmy, 'GeV')
    cv = getClosestValue (Dict, inmx, inmy)
    xV, yV, zV = [], [], []
    keys = Dict.keys()
    keys.sort()
    for xvalue in keys:
        yvalues = Dict[xvalue]
        for (yvalue, zvalue) in yvalues.items():
            xV.append (xvalue)
            yV.append (yvalue)
            zV.append (zvalue)
    # for i in range(len(xV)):
    #    print xV[i],yV[i],zV[i]
    try:
        ip2d = ip.interpolate.interp2d (xV, yV, zV, kind='cubic', fill_value=999999.)
        tmp = ip2d(mx, my)
        if len(tmp) != 1:
            logger.warning ("return value of interpolation is not a single-element array?")
            return None
        if abs (tmp[0] - cv) / cv < 0.3: return tmp[0]
        logger.warning ("when interpolating, would have had to interpolate too much. returning None.")
        return None
    except Exception, e:
        logger.warning ("when interpolating, caught exception: " + str(e))
        return None

    # print "mx=",mx,"my=",my,"res=",ip2d( mx, my )

def getUpperLimitFromDictionary (analysis, topo, mx=None, my=None, run=None, png=None, interpolate=False, expected=False):
    """ shouldnt have to call this directly. It's obtaining an upper limit from the python dictionary """
#    if interpolate:
#     print "[SMSResults.py] error: need to implement interpolation function for getUpperLimitFromDictionary"
#        import sys
#        sys.exit(0)
    Dict = smsHelpers.getUpperLimitDictionary (analysis, topo, run, expected=expected)
    if Dict == None: return Dict
    # # Dict=addunit ( Dict, "pb" )
    # print "[SMSResults.py] mx=",mx
    if rmvunit(mx, 'GeV') == None: return Dict
    # # return getClosestValue ( Dict, mx, my )
    return addunit (getInterpolatedUpperLimitDelaunay (Dict, mx, my), "pb")

def getSmartUpperLimit (analysis, topo, masses, massesbranch2=None, debug=False):
    """ returns the upper limit for analysis/topo, given an ordered sequence of
            the mass (mother, intermediate, LSP) """
    import smsInterpolation
    return smsInterpolation.upperLimit(analysis, topo, masses, debug)
#    return getUpperLimit ( analysis, topo, mx=masses[0], my=masses[-1], interpolate=True )

def getUpperLimit (analysis, topo, mx=None, my=None, run=None, png=None, interpolate=False, expected=False):
    """ get the upper limit for run/analysis/topo.
            return none if it doesnt exist.
            if mx and my are none, return the entire histogram,
            if mx and my are floats, return the upper limit at this
            point
            if png==True, return path of pngfile containing the histogram"""
    run = smsHelpers.getRun (analysis, run)
    if smsHelpers.hasDictionary (analysis, run):
        return  getUpperLimitFromDictionary (analysis, topo, mx, my, run, interpolate=interpolate, expected=expected)
    if png == True:
        pngfile = smsHelpers.getUpperLimitPng(analysis, topo, run)
        return pngfile
    if smsHelpers.hasHistogram (analysis, run):
        histo = smsHelpers.getUpperLimitFromHisto (analysis, topo, run, expected=expected)
        if rmvunit(mx, 'GeV') == None:
            return histo
        value = smsHelpers.getUpperLimitAtPoint (histo, mx, my, interpolate=interpolate)
        if value == 0.0: value = None    # 0.0 usually means out of bounds
        return addunit (value, "pb")
    logger.warning ("no upper limits found for %s" %analysis)
    return None 

def getEfficiency (analysis, topo, mx=None, my=None, run=None):
    """ get the efficiency for run/analysis/topo.
            return none if it doesnt exist.
            if mx and my are none, return the entire histogram,
            if mx and my are floats, return the upper limit at this
            point """
    run = smsHelpers.getRun (analysis, run)
    histo = smsHelpers.getEfficiencyHisto (analysis, topo, run)
    if mx == None:
        return histo
    value = smsHelpers.getEfficiencyAtPoint (histo, mx, my)
    return value

def getExplanationForLackOfUpperLimit (analysis, topo, mx=None, my=None, run=None, number=False):
    """if there's no upper limit, we want to know what's wrong.
       If number is false, return a text, if number is true
       return the error code
    """
    value = getUpperLimit (analysis, topo, run=run)
    msg = smsHelpers.getErrorMessage (value, mx, my)
    if number: return msg[0]
    return msg[1]

def getLumi (analysis, run=None):
    """ get the integrated luminosity for this analysis """
    lumifb = float(smsHelpers.getMetaInfoField (analysis, "lumi", run))
    return addunit (lumifb, "fb-1")

def getSqrts (analysis, run=None):
    """ get s_hat for this analysis """
    sqrts = smsHelpers.getMetaInfoField (analysis, "sqrts", run)
    try:
        return addunit (float(sqrts), "TeV")
    except:
        pass
    return sqrts

def getPAS (analysis, run=None):
    """ get the PAS for this analysis """
    return smsHelpers.getMetaInfoField (analysis, "pas", run)

def hasDictionary (analysis, run=None):
    """ are the upper limits available in dictionary format? """
    return smsHelpers.hasDictionary (analysis, run)

def getx (analysis, topo=None, run=None):
    """ get the description of the x-values for this analysis, if you supply a
            topo, then the return value is the x-values only for this topo """

    st = smsHelpers.getMetaInfoField (analysis, "x", run)
    if not st:
        return None
    st = st.split(',')
    d = {}
    for i in range(len(st)):
        l = st[i].split(':')
        x = l[1].split()
        d[l[0].replace(" ", "")] = x

    if topo:
        topo = topo.replace(" ", "")
        if not d or not d.has_key (topo): return None
        else: return d[topo]

    return d

def getFigures (analysis, run=None):
    """ get the figure number for this analysis """
    return smsHelpers.getMetaInfoField (analysis, "figures", run)

def getComment (analysis, run=None):
    """ an option comment? """
    return smsHelpers.getMetaInfoField (analysis, "comment", run)

conditions = {}

def getConditions (analysis, topo="all", fuzzy=True, run=None):
    """ get the conditions. if topo is "all",
            returns a dictionary, else it returns the condition
            only for the given topo, None if non-existent. """
    key = analysis + topo + str(fuzzy) + str(run)
    if conditions.has_key (key): return conditions[key]
    run = smsHelpers.getRun (analysis, run)
    if fuzzy: ret = smsHelpers.fuzzyconditions (analysis, run)
    else: ret = smsHelpers.conditions (analysis, run)
    if topo == "all":
        conditions[key] = ret
        return ret
    if not ret.has_key (topo):
        conditions[key] = None
        return None
    conditions[key] = ret[topo]
    return ret[topo]

constraints = {}

def getConstraints (analysis, topo="all", run=None):
    """ get the constraints. if topo is "all", 
            returns a dictionary, else it returns the constraint
            only for the given topo, None if non-existent. """
    key = analysis + topo + str(run)
    if constraints.has_key (key): return constraints[key]
    run = smsHelpers.getRun (analysis, run)
    ret = smsHelpers.constraints (analysis, run)
    if topo == "all":
        constraints[key] = ret
        return ret
    if not ret.has_key (topo):
        constraints[key] = None
        return None
    constraints[key] = ret[topo]
    return ret[topo]

def getJournal (analysis, run=None):
    """ get the journal of this analysis """
    return smsHelpers.getMetaInfoField (analysis, "journal", run)

def getBibtex (analysis, run=None):
    """ get the inspire page with the bibtex entry for this analysis """
    return smsHelpers.getMetaInfoField (analysis, "bibtex", run)

def getURL (analysis, run=None):
    """ get the URL for this analysis """
    return smsHelpers.getMetaInfoField (analysis, "url", run)

def hasURL (analysis, run=None):
    """ see if an URL is known """
    return smsHelpers.hasMetaInfoField (analysis, "url", run)

def exists(analysis, topo, run=None):
    """ check if the histogram ``limit_topo'' in 
     run/analysis/sms.py|root exists.
     For topologies with intermediate masses, check if all histograms (dictionaries)
     listed in the axes-information exist.
     If topo==None, simply check
     if run/analysis/sms.py|root exists. """
    import smsInterpolation
    run2 = smsHelpers.getRun(analysis, run)
    if not topo:
        import os
        Base = smsHelpers.Base
        rootfile = "%s/%s/%s/sms.root" % (Base, run2, analysis)
        pydict = "%s/%s/%s/sms.py" % (Base, run2, analysis)
        if os.path.exists(rootfile) or os.path.exists(pydict):
            return True
        else:
            return False
    axes = getaxes(analysis, topo)
    if not axes: return False
    hasDict = smsHelpers.hasDictionary (analysis, run2)
    for a in axes:
# print "a=",a
        mzname = None
        if a['mz'] and len(a['mz']): mzname = a['mz'][0]
        toponame = smsInterpolation.getHistName(topo, mzname)
        if hasDict:
# print "BBB 1 ana=%s run=%s run2=%s" % ( analysis,run,run2 )
            Dict = smsHelpers.getUpperLimitDictionary (analysis, toponame, run2)
            if not Dict or len(Dict) == 0: return False
            continue
        histo = smsHelpers.getUpperLimitFromHisto(analysis, toponame, run2)
        if not histo: return False

    return True


def getaxes (analysis, topo=None, run=None):
    """ get information about the histogram axes for this analysis: for each topo
        list of dictionary, each dictionary corresponds to one histogram, the key
        axes gives string (mx-my), the key mz gives information on other masses, if
        you supply a topo, returns list for this topo only. """
    if not exists(analysis, topo=None):    # # analysis is not known, we return None
        return None
    try:
        st = smsHelpers.getMetaInfoField (analysis, "axes", run)
    except MetaInfoError:
        logger.error("Meta info field 'axes' does not exist in %s." %analysis)
        st = None
    if not st:
        if not topo: return None    # cannot return default without info on topology
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

    if topo:
        topo = topo.replace(" ", "")
        if not d or not d.has_key (topo):
    # # topology does not exist, we return None
            return None
        else: return d[topo]

    return d

