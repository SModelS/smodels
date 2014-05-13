"""
.. module:: theory.theoryPrediction
   :synopsis: Provides a class to encapsulate the results of the computation of
              reference cross sections and related functions.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

import copy
from smodels.theory import clusterTools, crossSection, element
from smodels.theory.particleNames import elementsInStr
from smodels.theory.auxiliaryFunctions import cSim, cGtr  # pylint: disable=W0611
from smodels.theory.analysis import SRanalysis
from smodels.theory.analysis import ULanalysis
from smodels.theory.printer import Printer
import logging

logger = logging.getLogger(__name__)


class TheoryPrediction(object):
    """
    An instance of this class represents the results of the theory prediction
    for an analysis.
    
    """
    def __init__(self):
        self.analysis = None
        self.value = None
        self.conditions = None
        self.mass = None


class TheoryPredictionList(Printer):
    """
    An instance of this class represents the a collection of theory prediction
    objects.
    
    """
    def __init__(self, theoryPredictions):
        self._theoryPredictions = theoryPredictions


    def __iter__(self):
        for theoryPrediction in self._theoryPredictions:
            yield theoryPrediction


    def prepareData(self):
        """
        Select data preparation method through dynamic binding.
        
        """
        return Printer.prepareTheoryPredictionData(self)


def theoryPredictionFor(analysis, smsTopList, maxMassDist=0.2):
    """
    Compute theory predictions.
    
    Collect elements and efficiencies, combine the masses (if needed) and
    compute the conditions (if existing).
    
    :returns: list of TheoryPrediction objects
    
    """
    # Select elements constrained by analysis and apply efficiencies
    elements = _getElementsFrom(smsTopList, analysis)
    if len(elements) == 0:
        return None
    # Combine masses
    clusters = _combineElements(elements, analysis, maxDist=maxMassDist)
    # Collect results and evaluate conditions
    predictions = []
    for cluster in clusters:
        theoryPrediction = TheoryPrediction()
        theoryPrediction.analysis = analysis
        theoryPrediction.value = cluster.getTotalXSec()
        theoryPrediction.conditions = _evalConditions(cluster, analysis)
        theoryPrediction.mass = cluster.getAvgMass()
        predictions.append(theoryPrediction)

    if len(predictions) == 0:
        return None
    return TheoryPredictionList(predictions)


def _getElementsFrom(smsTopList, analysis):
    """
    Get elements, that are constrained by the analysis.
    
    Loop over all elements in smsTopList and returns the elements which are
    constrained by the analysis (have efficiency != 0). The elements weights
    are multiplied by their respective efficiency and the cross-sections not
    matching the analysis sqrts are removed.
    
    """
    elements = []
    for el in smsTopList.getElements():
        eff = analysis.getEfficiencyFor(el)
        if eff == 0.:
            continue
        element = el.copy()
        element.weight = crossSection.XSectionList()
        for xsec in el.weight:
            if xsec.info.sqrts == analysis.sqrts:
                element.weight.add(copy.deepcopy(xsec * eff))
        if len(element.weight) > 0:
            elements.append(element)

    return elements


def _combineElements(elements, analysis, maxDist):
    """
    Combine elements according to the analysis type.
    
    If analysis == upper limit type, group elements into mass clusters. If
    analysis == signal region type, group all elements into a single cluster.
    
    """
    if type(analysis) == type(SRanalysis()):
        clusters = [clusterTools.groupAll(elements)]
    elif type(analysis) == type(ULanalysis()):
        clusters = clusterTools.clusterElements(elements, analysis, maxDist)
    return clusters


def _evalConditions(cluster, analysis):
    """
    Evaluate the analysis conditions inside an element cluster.
    
    If analysis type == upper limit, evaluates the analysis conditions inside
    an element cluster.
    
    :retunrs: None, if analysis type == signal region
    
    """
    if type(analysis) == type(SRanalysis()):
        return None
    elif type(analysis) == type(ULanalysis()):
        if not analysis.conditions:
            return analysis.conditions
        conditions = {}
        # Loop over conditions
        for cond in analysis.conditions:
            # Get elements appearing in conditions
            condElements = [element.Element(elStr) \
                            for elStr in elementsInStr(cond)]
            newcond = cond
            for iel, el in enumerate(condElements):
                newcond = newcond.replace(str(el), "condElements[" + str(iel) +
                                          "].weight")
                for el1 in cluster.elements:
                    if el1.particlesMatch(el):
                        el.weight.combineWith(el1.weight)


            if newcond.find("Cgtr") >= 0:
                newcond = newcond.replace("Cgtr", "cGtr")
                logger.warning(analysis.label + " using deprecated function "
                               "'Cgtr'. Auto-replacing with 'cGtr'.")

            if newcond.find("Csim") >= 0:
                newcond = newcond.replace("Csim", "cSim")
                logger.warning(analysis.label + " using deprecated function "
                               "'Csim'. Auto-replacing with 'cSim'.")

            conditions[cond] = eval(newcond)

        return conditions
