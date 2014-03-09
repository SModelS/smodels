#!/usr/bin/env python

"""
.. module:: SMSResults
   :synopsis: Centralized facility to access the SMS results.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

from experiment import SMSHelpers
from tools.PhysicsUnits import addunit, rmvunit
from experiment.SMSResultsCollector import description
from tools import PhysicsUnits, VariousHelpers, RCFile
from experiment.experimentExceptions import MetaInfoError
logger = VariousHelpers.logging.getLogger(__name__)

def describeTx (topo, short=True):
    """ describe a Tx name, e.g. T2tt -> "#tilde{t} -> t #tilde{#chi}^{0} """
    ret = topo + ": " + description (topo, plot='ROOT', kerning=False, short=short).replace("pp #rightarrow ", "")
    ret = ret[:ret.find(";")]
    return ret

def getAllHistNames (ana, topo, run=None):
    """ for a given analysis, topology, return list of all available histograms"""
    from experiment import SMSInterpolation
    ret = []
    dic_list = getaxes(ana, topo, run)
    if not dic_list: return None    # topology not found
    if len(dic_list) == 1 and not dic_list[0]['mz']: return [SMSInterpolation.gethistname(topo, None)]
    for dic in dic_list:
        ret.append(SMSInterpolation.gethistname(topo, dic["mz"][0]))
    return ret

def setLogLevel (l=[ "error" ]):
    """ defines what is written out, and what isnt """
    SMSHelpers.logLevel = l

def setBase (base):
    """ just sets the base directory of the database """
    SMSHelpers.Base = base

def useUnits (b=True):
    PhysicsUnits.useUnits = b

def considerRuns(run=None):
    """ 
        defines what runs are to be considered when asking for results.

        :param run: a list of runs to be considered, e.g. [ '2012', '8TeV' ]). If None, all runs are taken into account.
        :type run: list or NoneType
    """
    from experiment import SMSResultsCollector
    allruns = ["8TeV", "ATLAS8TeV", "RPV8", "2012", "RPV7", "2011"]
    runsort = []
    if run:
        for r in allruns:
            if r in run:
                runsort.append(r)
        SMSHelpers.runs = runsort
        SMSResultsCollector.alldirectories = runsort
        for r in run:
            if not r in allruns:
                print "%s is not a run!!!" % r
    else:
        SMSHelpers.runs = allruns
        SMSResultsCollector.alldirectories = allruns

def verbosity (level="error"):
    if level == "error":
        SMSHelpers.verbose = False
    else:
        SMSHelpers.verbose = True

def getExclusion (analysis, topo, run=None):
    """ get the exclusions on the mother particles,
            for the summary plots """
    run = SMSHelpers.getRun (analysis, run)
    ex = SMSHelpers.motherParticleExclusions (analysis, run)
    for t in SMSHelpers.getPotentialNames (topo):
        if ex.has_key (t): return ex[t]
    return None

def getBinWidthX (analysis, topo, run=None):
    """ get the bin width in X """
    run = SMSHelpers.getRun (analysis, run)
    if SMSHelpers.hasDictionary (analysis, run):
        return addunit(25, "GeV")
    histo = SMSHelpers.getUpperLimitFromHisto (analysis, topo, run)
    if not histo: return None
    w = histo.GetXaxis().GetBinWidth(1)
    return addunit (w, "GeV")

def getLowX (analysis, topo, run=None):
    """ get the lower edge of the x axis """
    run = SMSHelpers.getRun (analysis, run)
    if SMSHelpers.hasDictionary (analysis, run):
        dic = SMSHelpers.getUpperLimitDictionary(analysis, topo, run)
        if not dic: return None
        minx = min(dic.keys()) - 12.5
        return addunit (minx, "GeV")
    histo = SMSHelpers.getUpperLimitFromHisto (analysis, topo, run)
    if not histo: return None
    w = histo.GetXaxis().GetBinLowEdge(histo.GetXaxis().GetFirst())
    return addunit (w, "GeV")

def getUpX(analysis, topo, run=None):
    """ get the upper edge of the x axis """
    run = SMSHelpers.getRun (analysis, run)
    if SMSHelpers.hasDictionary (analysis, run):
        dic = SMSHelpers.getUpperLimitDictionary(analysis, topo, run)
        if not dic: return None
        maxx = max(dic.keys()) + 12.5
        return addunit (maxx, "GeV")
    histo = SMSHelpers.getUpperLimitFromHisto (analysis, topo, run)
    if not histo: return None
    w = histo.GetXaxis().GetBinUpEdge(histo.GetXaxis().GetLast())
    return addunit (w, "GeV")

def getLowY (analysis, topo, run=None):
    """ get the lower edge of the y axis """
    run = SMSHelpers.getRun (analysis, run)
    if SMSHelpers.hasDictionary (analysis, run):
        dic = SMSHelpers.getUpperLimitDictionary(analysis, topo, run)
        if not dic: return None
        ylist = []
        for ydic in dic.values():
            ylist.extend(ydic.keys())
        miny = min(ylist) - 12.5
        return addunit(miny, "GeV")
    histo = SMSHelpers.getUpperLimitFromHisto (analysis, topo, run)
    if not histo: return None
    w = histo.GetYaxis().GetBinLowEdge(histo.GetYaxis().GetFirst())
    return addunit (w, "GeV")

def getUpY (analysis, topo, run=None):
    """ get the upper edge of the y axis """
    run = SMSHelpers.getRun (analysis, run)
    if SMSHelpers.hasDictionary (analysis, run):
        dic = SMSHelpers.getUpperLimitDictionary(analysis, topo, run)
        if not dic: return None
        ylist = []
        for ydic in dic.values():
            ylist.extend(ydic.keys())
        maxy = max(ylist) + 12.5
        return addunit(maxy, "GeV")
    histo = SMSHelpers.getUpperLimitFromHisto (analysis, topo, run)
    if not histo: return None
    w = histo.GetYaxis().GetBinUpEdge(histo.GetYaxis().GetLast())
    return addunit (w, "GeV")

def getBinWidthY (analysis, topo, run=None):
    """ get the bin width in Y """
    run = SMSHelpers.getRun (analysis, run)
    if SMSHelpers.hasDictionary (analysis, run):
        return addunit (25, "GeV")
    histo = SMSHelpers.getUpperLimitFromHisto (analysis, topo, run)
    if not histo: return None
    w = histo.GetYaxis().GetBinWidth(1)
    return addunit (w, "GeV")

def getExclusionLine(topo, ana, expected=False, plusminussigma=0, extendedinfo=False, xvalue=None, factor=1.0):
    """ get the exclusion line, as a TGraph """
    if xvalue == None: xvalue = ''
    from experiment import SMSResultsCollector
    ex = SMSResultsCollector.exclusionline(topo, ana, xvalue, factor, extendedinfo=extendedinfo, expected=expected, plusminussigma=plusminussigma)
    return ex

def getTopologies (analysis, run=None, allHistos=False):
    """ return all topologies that this analysis has results for
            if allHistos=True is selected: return all available histograms """
    run = SMSHelpers.getRun (analysis, run)
    if allHistos or not getaxes(analysis):
# we used the exclusion info to get the list
        x = SMSHelpers.motherParticleExclusions (analysis, run)
        return x.keys()
    topolist = []
    axes = getaxes(analysis)
    for topo in axes.keys():
        if exists(analysis, topo): topolist.append(topo)
    return topolist

def getRun (analysis, run=None):
    """ tell us, which run the results will be fetched for.
            None if the request cannot be met """
    run = SMSHelpers.getRun (analysis, run)
    return run

def getExperiment (analysis, run=None):
    """ return experiment name for given analysis
            for now: check if run is ATLAS8TeV, else return CMS """
    run1 = getRun(analysis, run)
    if run1 == "ATLAS8TeV": return "ATLAS"
    else: return "CMS"


def getPrettyName (analysis, run=None, latex=False):
    value = SMSHelpers.getMetaInfoField (analysis, "prettyname", run)
    if value == None or value == "": return analysis
    if not latex:
        value = value.replace("\\", "#")
    return value

def getAnalyses (topo, run=None):
    """return all analyses that have results for topo
    """
    import os
    runs = SMSHelpers.runs
    if run: runs = [ run ]
    analyses = {}
    for r in runs:
# # so thats the runs I really have to think about
        dirs = os.listdir ("%s/%s/" % (SMSHelpers.Base, r))
        for ana in dirs:
            if os.path.exists ("%s/%s/%s/info.txt" % (SMSHelpers.Base, r, ana)):
#                e=getExclusion ( ana, topo, r )
#                if e: analyses[ana]=True
                if exists(ana, topo, r): analyses[ana] = True

    return analyses.keys()

allresults = {}

def getAllResults (run=None, allHistos=False):
    """returns all analyses and the topologies they have results for
    """
    import os
    key = str(run)
    if allresults.has_key (key):
        return allresults[key]
    runs = SMSHelpers.runs
    if run: runs = [ run ]
    ret = {}
    for r in runs:
# # so thats the runs I really have to think about
        dirs = os.listdir ("%s/%s/" % (SMSHelpers.Base, r))
        for ana in dirs:
            if os.path.exists ("%s/%s/%s/info.txt" % (SMSHelpers.Base, r, ana)):
                topos = getTopologies (ana, run, allHistos=allHistos)
                ret[ana] = topos
    allresults[key] = ret
    return ret


dbresults = {}

def getDatabaseResults (run=None, category=None):
    """ returns all public analyses and the topologies we have constraints for """
    import os, copy
    key = str(run)+str(category)
    if dbresults.has_key (key):
        return dbresults[key]
    runs = SMSHelpers.runs
    if run: runs = [ run ]
    print "runs=",runs
    ret = {}
    for r in runs:
# # so thats the runs I really have to think about
        dirs = os.listdir ("%s/%s/" % (SMSHelpers.Base, r))
        for ana in dirs:
            if os.path.exists ("%s/%s/%s/info.txt" % (SMSHelpers.Base, r, ana)):
                topos = getConstraints (ana, run=r).keys()
                if not topos: continue
                if category:
                    if isPrivate(ana): continue
                    newTopos=copy.deepcopy(topos)
                    for t in newTopos:
                        if category not in getCategories(ana, topo=t):
                            topos.remove(t)
                if not topos: continue
                ret[ana] = topos
    dbresults[key] = ret
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
    Dict = SMSHelpers.getUpperLimitDictionary (analysis, topo, run, expected=expected)
    if Dict == None: return Dict
    # # Dict=addunit ( Dict, "pb" )
    # print "[SMSResults.py] mx=",mx
    if rmvunit(mx, 'GeV') == None: return Dict
    # # return getClosestValue ( Dict, mx, my )
    return addunit (getInterpolatedUpperLimitDelaunay (Dict, mx, my), "pb")

def getSmartUpperLimit (analysis, topo, masses, massesbranch2=None, debug=False):
    """ returns the upper limit for analysis/topo, given an ordered sequence of
            the mass (mother, intermediate, LSP) """
    import SMSInterpolation
    return SMSInterpolation.UpperLimit(analysis, topo, masses, debug)
#    return getUpperLimit ( analysis, topo, mx=masses[0], my=masses[-1], interpolate=True )

def getUpperLimit (analysis, topo, mx=None, my=None, run=None, png=None, interpolate=False, expected=False):
    """ get the upper limit for run/analysis/topo.
            return none if it doesnt exist.
            if mx and my are none, return the entire histogram,
            if mx and my are floats, return the upper limit at this
            point
            if png==True, return path of pngfile containing the histogram"""
    run = SMSHelpers.getRun (analysis, run)
    if SMSHelpers.hasDictionary (analysis, run):
        return  getUpperLimitFromDictionary (analysis, topo, mx, my, run, interpolate=interpolate, expected=expected)
    if png == True:
        pngfile = SMSHelpers.getUpperLimitPng(analysis, topo, run)
        return pngfile
    if SMSHelpers.hasHistogram (analysis, run):
        histo = SMSHelpers.getUpperLimitFromHisto (analysis, topo, run, expected=expected)
        if rmvunit(mx, 'GeV') == None:
            return histo
        value = SMSHelpers.getUpperLimitAtPoint (histo, mx, my, interpolate=interpolate)
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
    run = SMSHelpers.getRun (analysis, run)
    histo = SMSHelpers.getEfficiencyHisto (analysis, topo, run)
    if mx == None:
        return histo
    value = SMSHelpers.getEfficiencyAtPoint (histo, mx, my)
    return value

def getExplanationForLackOfUpperLimit (analysis, topo, mx=None, my=None, run=None, number=False):
    """if there's no upper limit, we want to know what's wrong.
       If number is false, return a text, if number is true
       return the error code
    """
    value = getUpperLimit (analysis, topo, run=run)
    msg = SMSHelpers.getErrorMessage (value, mx, my)
    if number: return msg[0]
    return msg[1]

def isPrivate(analysis, run=None):
    """Check if analysis is flagged as private.
         
       :param analysis: analysis to be checked for its public/private state
       :param run: run to be selected (default: None)
       :returns: True, if analysis is flagged as private, else False
    """
    try:
        return bool(int(SMSHelpers.getMetaInfoField(analysis, "private", run)))
    except MetaInfoError:
        logger.exception("Could not parse field 'private'.")


def isPrivateTopology(analysis, topology, run=None):
    """Check if topology of analysis is flagged as private.
    
       :param analysis: analysis to be selected
       :param topology: topology to be checked for its public/private state
       :param run: run to be selected (default: None)
       :returns: True, if analysis or topology is selected, else False
    """
    if isPrivate(analysis, run):
        return True
    else:
        try:
            return topology in SMSHelpers.getMetaInfoField(analysis, "private_topologies", run).split()
        except MetaInfoError:
            logger.exception("Could not parse field 'private_topologies'.")

def isPublic(analysis, run=None):
    """DEPRECATED - Check if analysis is NOT flagged as private.
         
       :param analysis: analysis to be checked for its public/private state
       :param run: specific run of the analysis to be checked (default: None)
       :returns: True, if analysis is NOT flagged as private, else False
    """
    return not isPrivate(analysis, run)

def hasDataPublished (analysis, run=None):
    """ has the analysis published their data in digitized form? """
    value = SMSHelpers.getMetaInfoField (analysis, "publisheddata", run)
    if type(value) == type("string"):
        if value == "None": return None
        if value == "True" or value == "1": return True
        if value == "False" or value == "0": return False
    if not value:
        return None
    try:
        return bool(value)
    except Exception:
        logger.exception("Could not parse field 'publisheddata'.")
    return None

def getLumi (analysis, run=None):
    """ get the integrated luminosity for this analysis """
    lumifb = float(SMSHelpers.getMetaInfoField (analysis, "lumi", run))
    return addunit (lumifb, "fb-1")

def getSqrts (analysis, run=None):
    """ get s_hat for this analysis """
    sqrts = SMSHelpers.getMetaInfoField (analysis, "sqrts", run)
    try:
        return addunit (float(sqrts), "TeV")
    except:
        pass
    return sqrts

def getPAS (analysis, run=None):
    """ get the PAS for this analysis """
    return SMSHelpers.getMetaInfoField (analysis, "pas", run)

def getOrder (analysis, run=None):
    """ get the order in perturbation theory that the exclusion lines correspond
            with """
    return SMSHelpers.getMetaInfoField (analysis, "order", run)

def hasDictionary (analysis, run=None):
    """ are the upper limits available in dictionary format? """
    return SMSHelpers.hasDictionary (analysis, run)

def getx (analysis, topo=None, run=None):
    """ get the description of the x-values for this analysis, if you supply a
            topo, then the return value is the x-values only for this topo """

    st = SMSHelpers.getMetaInfoField (analysis, "x", run)
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
    return SMSHelpers.getMetaInfoField (analysis, "figures", run)

def getComment (analysis, run=None):
    """ an option comment? """
    return SMSHelpers.getMetaInfoField (analysis, "comment", run)

conditions = {}

def getConditions (analysis, topo="all", fuzzy=True, run=None):
    """ get the conditions. if topo is "all",
            returns a dictionary, else it returns the condition
            only for the given topo, None if non-existent. """
    key = analysis + topo + str(fuzzy) + str(run)
    if conditions.has_key (key): return conditions[key]
    run = SMSHelpers.getRun (analysis, run)
    if fuzzy: ret = SMSHelpers.fuzzyconditions (analysis, run)
    else: ret = SMSHelpers.conditions (analysis, run)
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
    run = SMSHelpers.getRun (analysis, run)
    ret = SMSHelpers.constraints (analysis, run)
    if topo == "all":
        constraints[key] = ret
        return ret
    if not ret.has_key (topo):
        constraints[key] = None
        return None
    constraints[key] = ret[topo]
    return ret[topo]

categories = {}

def getCategories (analysis, topo="all", run=None):
    """ get the categories. if topo is "all", 
            returns a dictionary, else it returns the constraint
            only for the given topo, None if non-existent. """
    key = analysis + topo + str(run)
    if categories.has_key (key): return categories[key]
    run = SMSHelpers.getRun (analysis, run)
    ret = SMSHelpers.categories (analysis, run)
    if not ret: return None
    if topo == "all":
        categories[key] = ret
        return ret
    if not ret.has_key (topo):
        categories[key] = None
        return None
    categories[key] = ret[topo]
    return ret[topo]

def getRequirement (analysis, run=None):
    """ any requirements that come with this analysis? (e.g. onshellness) """
    return SMSHelpers.getMetaInfoField (analysis, "requires", run)

def getCheckedBy (analysis, run=None):
    """ has the result been validated? """
    return SMSHelpers.getMetaInfoField (analysis, "checked", run)

def getJournal (analysis, run=None):
    """ get the journal of this analysis """
    return SMSHelpers.getMetaInfoField (analysis, "journal", run)

def getBibtex (analysis, run=None):
    """ get the inspire page with the bibtex entry for this analysis """
    return SMSHelpers.getMetaInfoField (analysis, "bibtex", run)

def getURL (analysis, run=None):
    """ get the URL for this analysis """
    return SMSHelpers.getMetaInfoField (analysis, "url", run)

def hasURL (analysis, run=None):
    """ see if an URL is known """
    return SMSHelpers.hasMetaInfoField (analysis, "url", run)

def getContact (analysis, run=None):
    """ get the contact for this analysis """
    return SMSHelpers.getMetaInfoField (analysis, "contact", run)

def getPerturbationOrder (analysis, run=None):
    """ get the contact for this analysis """
    return SMSHelpers.getMetaInfoField (analysis, "order", run)

def particles(topo, plot='ROOT'):
    """return the production mode for a given topology:
         latex code either compatible to ROOT or Python"""
    if SMSHelpers.dicparticle.has_key(topo):
        part = SMSHelpers.dicparticle[topo]
        if plot == 'ROOT':
    # print "[SMSResults.py] debug",part
            return part
        if plot == 'python':
            return part.replace('#', '\\')
    else:
        return None

def particleName(topo):
    """return the production mode for a given topology:
         write out the name in plain letters, no latex """
    if topo[:2] == "TGQ": return "associate"
    if topo == "TChiSlep" or topo == "TChiNuSlep": return "chargino/neutralino"
    if topo == "TChiSlepSlep": return "chargino/neutralino"
    if topo == "TChiWZ": return "chargino/neutralino"
    if topo[:4] == "TChi": return "chargino/neutralino"
    if not SMSHelpers.dicparticle.has_key(topo):
        return "???"
    part = SMSHelpers.dicparticle[topo].replace("#tilde", "").replace("{", "").replace("}", "")
    if part == "g": part = "gluino"
    if part == "b": part = "sbottom"
    if part == "t": part = "stop"
    if part == "q": part = "squark"
    return part

def massDecoupling_ (topo):
    if topo == "T2tt":
        return "m(#tilde{g},#tilde{q})>>m(#tilde{t})"
    if topo == "T2FVttcc":
        return "m(#tilde{g},#tilde{q})>>m(#tilde{t})"
    if topo == "T2bb":
        return "m(#tilde{g},#tilde{q})>>m(#tilde{b})"
    if topo == "TChiSlep":
# return "m(#tilde{g}),m(#tilde{q})>>m(#tilde{#chi}^{0}_{2}),m(#tilde{#chi}^{#pm})"
        return "m(#tilde{g},m(#tilde{q})>>m(#tilde{#chi}^{0}_{2},#tilde{#chi}^{#pm})"
    if topo == "TChiSlepSlep":
# return "m(#tilde{g}),m(#tilde{q})>>m(#tilde{#chi}^{0}_{2})"
        return "m(#tilde{g},#tilde{q})>>m(#tilde{#chi}^{0}_{2},#tilde{#chi}^{#pm})"
    if topo == "TChiNuSlep":
        return "m(#tilde{g},#tilde{q})>>m(#tilde{#chi}^{0}_{2},#tilde{#chi}^{#pm})"
# return "m(#tilde{g}),m(#tilde{q})>>m(#tilde{#chi}^{0}_{2}),m(#tilde{#chi}^{#pm})"
    if topo == "TChiwz":
        return "m(#tilde{g},#tilde{q})>>m(#tilde{#chi}^{0}_{2},#tilde{#chi}^{#pm})"
    T2 = topo[:2]
    if T2 == "T1" or T2 == "T3" or T2 == "T5":
        return "m(#tilde{q})>>m(#tilde{g})"
    if T2 == "T2" or T2 == "T4" or T2 == "T6":
        return "m(#tilde{g})>>m(#tilde{q})"
    return ""

def massDecoupling (topo, plot='ROOT', kerning=True):
    """ describe the assumed mass decoupling """
    md = massDecoupling_ (topo)
    if kerning:
        md = md.replace(">>", ">#kern[-.2]{>}")
    if plot != 'ROOT':
        md = '$' + md.replace('#', '\\') + '$'
    return md



def exists(analysis, topo, run=None):
    """ check if the histogram run/analysis/sms.root(limit_topo) exists.
                for topologies with intermediate masses, check if all histograms (dictionaries)
                listed in the axes-information exist."""
    """ or sms.py, if it's stored in dictionaries. If topo==None, simply check
            if run/analysis/sms.root or run/analysis/sms.py exists. """
    import SMSInterpolation
    run2 = SMSHelpers.getRun(analysis, run)
    if not topo:
        import os
        Base = SMSHelpers.Base
        rootfile = "%s/%s/%s/sms.root" % (Base, run2, analysis)
        pydict = "%s/%s/%s/sms.py" % (Base, run2, analysis)
        if os.path.exists(rootfile) or os.path.exists(pydict):
            return True
        else:
            return False
    axes = getaxes(analysis, topo)
    if not axes: return False
    hasDict = SMSHelpers.hasDictionary (analysis, run2)
    for a in axes:
# print "a=",a
        mzname = None
        if a['mz'] and len(a['mz']): mzname = a['mz'][0]
        toponame = SMSInterpolation.gethistname(topo, mzname)
        if hasDict:
# print "BBB 1 ana=%s run=%s run2=%s" % ( analysis,run,run2 )
            Dict = SMSHelpers.getUpperLimitDictionary (analysis, toponame, run2)
            if not Dict or len(Dict) == 0: return False
            continue
        histo = SMSHelpers.getUpperLimitFromHisto(analysis, toponame, run2)
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
        st = SMSHelpers.getMetaInfoField (analysis, "axes", run)
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
