# !/usr/bin/env python3

"""
.. module:: coverage
   :synopsis: Definitions of classes used to find, group and format missing topologies

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from smodels.tools.physicsUnits import fb
from smodels.experiment.reweighting import reweightFactorFor
from bisect import bisect
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.theory.element import Element
from smodels.experiment.defaultFinalStates import (WList, lList, tList, taList, nuList, jetList,
                                                   MET, HSCP, RHadronG, RHadronQ)


#  Default definitions for the uncovered topology categories/groups:
#  Element filters for each group:
#  (it should be a function which takes an Element object as input
#  and returns True if the element belongs to the group and False otherwise)
filtersDefault = {'missing (prompt)': lambda el: not ('prompt' in el.coveredBy),
                  'missing (displaced)': lambda el: not ('displaced' in el.coveredBy),
                  # 'missing (long cascade)' : lambda el: (not el.coveredBy) and el._getLength() > 3,
                  'missing (all)': lambda el: (not el.coveredBy),
                  'outsideGrid (all)': lambda el: (el.coveredBy and not el.testedBy)}

# Description for each group (optional and only used for printing):
# (if not defined, the group label will be used instead)
descriptionDefault = {'missing (prompt)': 'missing topologies with prompt decays',
                      'missing (displaced)': 'missing topologies with displaced decays',
                      # 'missing (long cascade)' : 'missing topologies with long cascade decays',
                      'missing (all)': 'missing topologies',
                      'outsideGrid (all)': 'topologies outside the grid'}

# Default BSM inclusive particles used for grouping final states:
bsmDefault = [MET, HSCP, RHadronG, RHadronQ]

# Defaul SM inclusive particles used for grouping final states
smDefault = [WList, lList, tList, taList, nuList, jetList]

# Weight factors for each group:
# (it should be a function which takes an Element object as input
# and returns the reweighting factor to be applied to the element weight. It is relevant if only
# the fraction of the weight going into prompt or displaced decays is required)
factorsDefault = {}
for key in filtersDefault:
    if 'prompt' in key.lower():
        factorsDefault[key] = lambda el: reweightFactorFor(el, 'prompt')
    elif 'displaced' in key.lower():
        factorsDefault[key] = lambda el: reweightFactorFor(el, 'displaced')
    else:
        # If not specified assumed all fractions
        # (note that we can not include any long-lived fraction since this is already included in
        # the topologies where the meta-stable particle appears as a final state,
        # so the total is = (fraction of all decays being prompt)
        # + (fraction of at least one displaced decay and no long-lived decays)
        factorsDefault[key] = lambda el: reweightFactorFor(el, 'prompt') + reweightFactorFor(el, 'displaced')


class Uncovered(object):
    """
    Wrapper object for defining and holding a list of coverage groups  (UncoveredGroup objects).

    The class builds a series of UncoveredGroup objects and stores them.
    """

    def __init__(self, topDict, sqrts=None, sigmacut=0*fb,
                 groupFilters=filtersDefault,
                 groupFactors=factorsDefault,
                 groupdDescriptions=descriptionDefault,
                 smFinalStates=smDefault, bsmFinalStates=bsmDefault):
        """
        Inititalize the object.

        :param topDict: TopologyDict object used to select elements from.
        :param sqrts: Value (with units) for the center of mass energy used to compute the missing cross sections.
                     If not specified the largest value available will be used.
        :param sigmacut: Minimum cross-section/weight value (after applying the reweight factor)
                       for an element to be included. The value should in fb (unitless)
        :param groupFilters: Dictionary containing the groups' labels and the method for selecting
                            elements.
        :param groupFactors: Dictionary containing the groups' labels and the method for reweighting
                            cross sections.
        :param groupdDescriptions: Dictionary containing the groups' labels and strings describing the group
                                  (used for printout)
        :param smFinalStates: List of (inclusive) Particle or MultiParticle objects used for grouping SM
                              particles when creating FinalStateElements.
        :param bsmFinalSates: List of (inclusive) Particle or MultiParticle objects used for grouping
                              BSM particles when creating FinalStateElements.
        """

        # Sanity checks:
        if not isinstance(groupFilters, dict):
            raise SModelSError("groupFilters input should be a Dictionary and not %s" % type(groupFilters))
        if not isinstance(groupFactors, dict):
            raise SModelSError("groupFactors input should be a Dictionary and not %s" % type(groupFactors))
        if sorted(groupFilters.keys()) != sorted(groupFactors.keys()):
            raise SModelSError("Keys in groupFilters and groupFactors do not match")
        if any(not hasattr(gFilter, '__call__') for gFilter in groupFilters.values()):
            raise SModelSError("Group filters must be callable methods")

        if sqrts is None:
            sqrts = max([xsec.info.sqrts for xsec in topDict.getTotalWeight()])
        else:
            sqrts = sqrts

        # Store the relevant element cross-sections to improve performance:
        for el in topDict.getElements():
            xsec = el.weightList.getXsecsFor(sqrts)
            if xsec:
                el._totalXsec = xsec[0].value.asNumber(fb)
            else:
                el._totalXsec = 0.

        self.groups = []
        # Create each uncovered group and get the topologies from topDict
        for gLabel, gFilter in groupFilters.items():
            # Initialize the uncovered topology list:
            uncoveredTopos = UncoveredGroup(label=gLabel, elementFilter=gFilter,
                                            reweightFactor=groupFactors[gLabel],
                                            smFinalStates=smFinalStates,
                                            bsmFinalStates=bsmFinalStates,
                                            sqrts=sqrts, sigmacut=sigmacut.asNumber(fb))
            if groupdDescriptions and gLabel in groupdDescriptions:
                uncoveredTopos.description = groupdDescriptions[gLabel]
            else:
                uncoveredTopos.description = gLabel
            # Fill the list with the elements in topDict:
            uncoveredTopos.getElementsFrom(topDict)
            self.groups.append(uncoveredTopos)

    def getGroup(self, label):
        """
        Returns the group with the required label. If not found, returns None.

        :param label: String corresponding to the specific group label

        :return: UncoveredGroup object which matches the label
        """

        for group in self.groups:
            if group.label == label:
                return group

        return None


class UncoveredGroup(object):
    """
    Holds information about a single coverage group: criteria for selecting and grouping elements,
    function for reweighting cross sections, etc.
    """

    def __init__(self, label, elementFilter, reweightFactor,
                 smFinalStates, bsmFinalStates, sqrts, sigmacut=0.):
        """
        :param label: Group label
        :param elementFilter: Function which takes an element as argument and returns True (False) if
                             the element should (not) be selected.
        :param reweightFactor: Function which takes an element as argument and returns the reweighting
                              factor to be applied to the element weight.
        :param smFinalStates: List of Particle/MultiParticle objects used to group Z2-even particles appearing
                            in the final state
        :param bsmFinalStates: List of Particle/MultiParticle objects used to group Z2-odd particles appearing
                            in the final state
        :param sqrts: Value (with units) for the center of mass energy used to compute the missing cross sections.
                     If not specified the largest value available will be used.
        :param sigmacut: Minimum cross-section/weight value (after applying the reweight factor)
                       for an element to be included. The value should in fb (unitless)
        """

        self.finalStateElements = []
        self.smFinalStates = smFinalStates
        self.bsmFinalStates = bsmFinalStates
        self.sqrts = sqrts
        self.sigmacut = sigmacut
        self.label = label
        self.elementFilter = elementFilter
        self.reweightFactor = reweightFactor

    def __str__(self):
        return self.label

    def __repr__(self):
        return str(self)

    def getElementsFrom(self, topDict):
        """
        Select the elements from topDict according to self.elementFilter
        and build FinalStateElements from the selected elements.
        The FinalStateElements weights corresponds to the missing cross-section
        with double counting from compressed elements already accounted for.
        """

        # First select all elements according to the filter (type of uncovered/missing topology):
        elementList = [el for el in topDict.getElements() if self.elementFilter(el)]
        # Get missing xsections including the reweight factor:
        missingXandEls = [(self.getMissingX(el)*self.reweightFactor(el), el) for el in elementList]
        # Only keep the elements above sigmacut:
        missingXandEls = [x for x in missingXandEls[:] if x[0] > self.sigmacut]
        # Sort according to largest missingX, smallest size and largest ID
        missingXandEls = sorted(missingXandEls, key=lambda pt: (pt[0], -pt[1].canonName, -pt[1].elID), reverse=True)

        # Split lists of elements and missingX:
        missingXsecs = [pt[0] for pt in missingXandEls]
        elementList = [pt[1] for pt in missingXandEls]

        # Remove all elements which are related to each other in order to avoid double counting
        # (keep always the first appearance in the list, so we always keep the ones with largest missing xsec)
        elementListUnique = []
        missingXsecsUnique = []
        ancestors = set()  # Keep track of all the ancestor ids of the elements in the unique list
        for i, element in enumerate(elementList):
            # Get ancestor object ids
            # (safer than element.elID since some elements can share elID = 0 if they
            # were never directly inserted in the TopologyList)
            ancestorsIDs = set([id(element)] + [id(el) for el in element.getAncestors()])
            # If the element has any common ancestor with any of the previous elements
            # or if it is an ancestor of any of the previous elements
            # skip it to avoid double counting
            if ancestors.intersection(ancestorsIDs):
                continue
            elementListUnique.append(element)
            missingXsecsUnique.append(missingXsecs[i])
            ancestors = ancestors.union(ancestorsIDs)

        # Now that we only have unique elements with their effective missing cross-sections
        # we create General Elements out of them
        for iel, el in enumerate(elementListUnique):
            missingX = missingXsecsUnique[iel]
            if not missingX:
                continue
            self.addToFSElements(el, missingX)

        # Finally sort general elements by their missing cross-section:
        self.finalStateElements = sorted(self.finalStateElements, key=lambda fsEl: fsEl.missingX, reverse=True)

    def getMissingX(self, element):
        """
        Calculate total missing cross section of an element, by recursively checking if its
        mothers already appear in the list.
        :param element: Element object

        :returns: missing cross section without units (in fb)
        """

        ancestorList = element.getAncestors()
        alreadyChecked = []  # keep track of which elements have already been checked
        # Get the (pre-loaded) total element weight in fb:
        missingX = element._totalXsec
        if not missingX:
            return 0.
        overlapXsec = 0.
        for ancestor in ancestorList:  # check all ancestors (ancestorList is sorted by generation)
            # Skip entries which correspond to the element itself
            if ancestor is element:
                continue
            # Make sure we do not subtract the same mother twice
            if any(ancestor is el for el in alreadyChecked):
                continue
            alreadyChecked.append(ancestor)
            # Check if ancestor passes the group filter (if it has been covered/tested or not):
            if self.elementFilter(ancestor):
                continue
            # Subtract the weight of the ancestor and skip checking for all the older family tree of the ancestor
            # (avoid subtracting the same weightList twice from mother and grandmother for instance)
            # (since the ancestorList is sorted by generation, the mother always
            # appears before the grandmother in the list)
            alreadyChecked += ancestor.getAncestors()
            if not hasattr(ancestor, '_totalXsec'):
                xsec = ancestor.weightList.getXsecsFor(self.sqrts)
                ancestor._totalXsec = xsec[0].value.asNumber(fb)
            overlapXsec += ancestor._totalXsec

        return missingX-overlapXsec

    def getTotalXSec(self, sqrts=None):
        """
        Calculate total missing topology cross section at sqrts. If no sqrts is given use self.sqrts
        :ivar sqrts: sqrts
        """
        xsec = 0.
        if not sqrts:
            sqrts = self.sqrts
        for fsEl in self.finalStateElements:
            xsec += fsEl.missingX
        return xsec

    def addToFSElements(self, el, missingX):
        """
        Adds an element to the list of missing topologies = final state elements.
        If the element contributes to an element that is already in the list, add element and weight to
        the element.
        :parameter el: element to be added
        :parameter missingX: missing cross-section for the element (in fb)
        """

        newGenEl = FinalStateElement(el, missingX, self.smFinalStates, self.bsmFinalStates)

        index = bisect(self.finalStateElements, newGenEl)
        if index != len(self.finalStateElements) and self.finalStateElements[index] == newGenEl:
            self.finalStateElements[index]._contributingElements.append(el)
            self.finalStateElements[index].missingX += missingX
        else:
            self.finalStateElements.insert(index, newGenEl)


class FinalStateElement(Element):
    """
    This class represents a simplified element which does
    only holds information about the final states. It holds a simple
    tree with one root (PV), having the final state nodes as its
    daughters.
    """

    def __new__(self, el, missingX=None,
                smFinalStates=smDefault,
                bsmFinalStates=bsmDefault):

        if missingX is None:
            missingX = el.weightList.getMaxXsec()

        # Get an element holding only the final states
        finalStatesEl = el.compressToFinalStates()
        finalStatesEl.missingX = missingX
        # Replace particles by inclusive particles, if possible:
        for node in finalStatesEl.tree.nodes:
            if node == finalStatesEl.tree.root:
                continue
            ptc = node.particle
            if ptc.isSM:
                for smFS in smFinalStates:
                    if ptc == smFS:
                        node.particle = smFS
                        break
            else:
                for bsmFS in bsmFinalStates:
                    if ptc == bsmFS:
                        node.particle = bsmFS
                        break

        finalStatesEl.tree.sort(force=True)
        finalStatesEl.tree.setGlobalProperties()
        finalStatesEl._contributingElements = [el]

        return finalStatesEl
