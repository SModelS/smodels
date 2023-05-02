#!/usr/bin/env python3

"""
.. module:: theoryPredictionsCombiner
   :synopsis: a module that deals with combining signal regions from
              different analyses, offering the same API as TheoryPrediction.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Jamie Yellen <j.yellen.1@research.gla.ac.uk>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from smodels.theory.theoryPrediction import TheoryPrediction
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
from smodels.tools.statsTools import StatsComputer
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError


class TheoryPredictionsCombiner(TheoryPrediction):
    """
    Facility used to combine theory predictions from different analyes.
    If a list with a single TheoryPrediction is given, return the TheoryPrediction
    object.
    """

    def __new__(cls,theoryPredictions: list, slhafile=None, deltas_rel=None):
        """
        If called with a list containing a single TheoryPrediction, return the TheoryPrediction object.
        Otherwise, create a TheoryPredictionsCombiner object.
        """

        if len(theoryPredictions) == 1:
            return theoryPredictions[0]
        else:
            tpCombiner = super(TheoryPredictionsCombiner, cls).__new__(cls)
            return tpCombiner

    def __init__(self, theoryPredictions: list, slhafile=None, deltas_rel=None):
        """
        Constructor.
        :param theoryPredictions: the List of theory predictions
        :param slhafile: optionally, the slhafile can be given, for debugging
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        """

        if len(theoryPredictions) == 0:
            raise SModelSError("asking to combine zero theory predictions")
        
        self.theoryPredictions = theoryPredictions
        self.slhafile = slhafile
        if deltas_rel is None:
            from smodels.tools.runtime import _deltas_rel_default

            deltas_rel = _deltas_rel_default
        self.deltas_rel = deltas_rel
        self.cachedObjs = {False: {}, True: {}, "posteriori": {}}
        self.cachedNlls = {False: {}, True: {}, "posteriori": {}}
        self._statsComputer = None

    @classmethod
    def selectResultsFrom(cls, theoryPredictions, anaIDs):
        """
        Select the results from theoryPrediction list which match one
        of the IDs in anaIDs. If there are multiple predictions for the
        same ID for which a likelihood is available, it gives priority
        to the ones with largest expected r-values.

        :param theoryPredictions: list of TheoryPrediction objects
        :param anaIDs: list with the analyses IDs (in string format) to be combined
        :return: a TheoryPredictionsCombiner object for the selected predictions.
                 If no theory prediction was selected, return None.
        """

        # First select the theory predictions which correspond to the analyses to be combined
        filteredTPs = [tp for tp in theoryPredictions if tp.analysisId() in anaIDs]
        filteredIDs = set([tp.analysisId() for tp in filteredTPs])
        # Now remove results with no likelihood available
        selectedTPs = [tp for tp in filteredTPs if tp.likelihood() is not None]
        selectedIDs = set([tp.analysisId() for tp in selectedTPs])
        # Warn the user concerning results with no likelihoods:
        for anaID in filteredIDs.difference(selectedIDs):
            logger.info(
                "No likelihood available for %s. This analysis will not be used in analysis combination."
                % anaID
            )
        # If no results are available, return None
        if len(selectedTPs) == 0:
            return None

        # Define a hierarchy for the results:
        priority = {"combined": 2, "efficiencyMap": 1, "upperLimit": 0}
        # Now sort by highest priority and then by highest expected r-value:
        selectedTPs = sorted(
            selectedTPs, key=lambda tp: (priority[tp.dataType()], tp.getRValue(expected=True))
        )
        # Now get a single TP for each result
        # (the highest ranking analyses come last and are kept in the dict)
        uniqueTPs = {tp.analysisId(): tp for tp in selectedTPs}
        uniqueTPs = list(uniqueTPs.values())

        combiner = cls(uniqueTPs)
        return combiner

    def dataId(self):
        """
        Return a string with the IDs of all the datasets used in the combination.
        """
        ids = [str(tp.dataset.getID()) for tp in self.theoryPredictions]
        ret = ",".join(ids)

        return ret

    def analysisId(self):
        """
        Return a string with the IDs of all the experimental results used in the combination.
        """

        ret = ",".join(sorted([tp.analysisId() for tp in self.theoryPredictions]))

        return ret

    def dataType(self, short=False):
        """
        Return its type (combined)
        :param: short, if True, return abbreviation (anacomb)
        """
        if short:
            return "comb"
        else:
            return "combined"

    def totalXsection(self):
        ret = 0.0 * fb
        if self.theoryPredictions is not None:
            for tp in self.theoryPredictions:
                ret += tp.xsection.value
        return ret

    def getmaxCondition(self):
        """
        Returns the maximum xsection from the list conditions

        :returns: maximum condition xsection (float)
        """
        conditions = [tp.getmaxCondition() for tp in self.theoryPredictions]
        return max(conditions)

    def setStatsComputer(self):
        """
        Creates and instance of StatsComputer depending on the
        type of TheoryPrediction/dataset. In case it is not possible
        to define a statistical computer (upper limit result or no expected
        upper limits), set the computer to 'N/A'.
        """

        # First make sure all theory predictions in the combiner
        # have well-defined stats models
        if any(tp.statsComputer == 'N/A' for tp in self.theoryPredictions):
            computer = 'N/A'
        else:
            computer = StatsComputer.forAnalysesComb(self.theoryPredictions, self.deltas_rel)     
            
        self._statsComputer = computer

    def getLlhds(self,muvals,expected=False,normalize=True):
        """
        Compute the likelihoods for the individual analyses and the combined
        likelihood.
        Returns a dictionary with the analysis IDs as keys and the likelihood values as values.

        :param muvals: List with values for the signal strenth for which the likelihoods must
                       be evaluated.
        :param expected: If True returns the expected likelihood values.
        :param normalize: If True normalizes the likelihood by its integral over muvals.
        """

        return self.statsComputer.likelihoodComputer.getLlhds(muvals,expected,normalize)

    def describe(self):
        """returns a string containing a list of all analysisId and dataIds"""
        ids = []
        for tp in self.theoryPredictions:
            ids.append(f"{tp.analysisId()}:{tp.dataId()}")
        return f"SRs: {', '.join(ids)}"


