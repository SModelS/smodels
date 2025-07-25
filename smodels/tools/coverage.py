# !/usr/bin/env python3

"""
.. module:: coverage
   :synopsis: Definitions of classes used to find, group and format missing topologies

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from bisect import bisect
from itertools import product
from smodels.base.physicsUnits import fb
from smodels.base.exceptions import SModelSBaseError as SModelSError
from smodels.experiment.reweighting import reweightFactorFor
from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.defaultFinalStates import (WList, lList, tList, taList, nuList, jetList,
                                                   MET, HSCP, RHadronG, RHadronQ,anyBSM)


#  Default definitions for the uncovered topology categories/groups:
#  SMS filters for each group:
#  (it should be a function which takes an SMS object as input
#  and returns True if the SMS belongs to the group and False otherwise)
filtersDefault = {'missing (prompt)': lambda sms: not ('prompt' in sms.coveredBy),
                  'missing (displaced)': lambda sms: not ('displaced' in sms.coveredBy),
                  # 'missing (long cascade)' : lambda el: (not el.coveredBy) and el._getLength() > 3,
                  'missing (all)': lambda sms: (not sms.coveredBy),
                  'outsideGrid (all)': lambda sms: (sms.coveredBy and not sms.testedBy)}

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
# (it should be a function which takes an SMS object as input
# and returns the reweighting factor to be applied to the SMS weight. It is relevant if only
# the fraction of the weight going into prompt or displaced decays is required)
factorsDefault = {}
for key in filtersDefault:
    if 'prompt' in key.lower():
        factorsDefault[key] = lambda sms: reweightFactorFor(sms, 'prompt')
    elif 'displaced' in key.lower():
        factorsDefault[key] = lambda sms: reweightFactorFor(sms, 'displaced')
    else:
        # If not specified assumed all fractions
        # (note that we can not include any long-lived fraction since this is already included in
        # the topologies where the meta-stable particle appears as a final state,
        # so the total is = (fraction of all decays being prompt)
        # + (fraction of at least one displaced decay and no long-lived decays)
        factorsDefault[key] = lambda sms: reweightFactorFor(sms, 'prompt') + reweightFactorFor(sms, 'displaced')


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

        :param topDict: TopologyDict object used to select SMS from.
        :param sqrts: Value (with units) for the center of mass energy used to compute the missing cross sections.
                     If not specified the largest value available will be used.
        :param sigmacut: Minimum cross-section/weight value (after applying the reweight factor)
                       for an SMS to be included. The value should in fb (unitless)
        :param groupFilters: Dictionary containing the groups' labels and the method for selecting
                            SMS.
        :param groupFactors: Dictionary containing the groups' labels and the method for reweighting
                            cross sections.
        :param groupdDescriptions: Dictionary containing the groups' labels and strings describing the group
                                  (used for printout)
        :param smFinalStates: List of (inclusive) Particle or MultiParticle objects used for grouping SM
                              particles when creating FinalStateSMS.
        :param bsmFinalSates: List of (inclusive) Particle or MultiParticle objects used for grouping
                              BSM particles when creating FinalStateSMS.
        """

        # Sanity checks:
        if not isinstance(groupFilters, dict):
            raise SModelSError(f"groupFilters input should be a Dictionary and not {type(groupFilters)}")
        if not isinstance(groupFactors, dict):
            raise SModelSError(f"groupFactors input should be a Dictionary and not {type(groupFactors)}")
        if sorted(groupFilters.keys()) != sorted(groupFactors.keys()):
            raise SModelSError("Keys in groupFilters and groupFactors do not match")
        if any(not hasattr(gFilter, '__call__') for gFilter in groupFilters.values()):
            raise SModelSError("Group filters must be callable methods")

        if sqrts is None:
            sqrts = max([xsec.info.sqrts for xsec in topDict.getTotalWeight()])
        else:
            sqrts = sqrts

        # Store the relevant SMS cross-sections to improve performance:
        for sms in topDict.getSMSList():
            xsec = sms.weightList.getXsecsFor(sqrts)
            if xsec:
                sms._totalXsec = xsec[0].value.asNumber(fb)
            else:
                sms._totalXsec = 0.

        self.groups = []
        # Create each uncovered group and get the topologies from topDict
        for gLabel, gFilter in groupFilters.items():
            # Initialize the uncovered topology list:
            uncoveredTopos = UncoveredGroup(label=gLabel, smsFilter=gFilter,
                                            reweightFactor=groupFactors[gLabel],
                                            smFinalStates=smFinalStates,
                                            bsmFinalStates=bsmFinalStates,
                                            sqrts=sqrts, sigmacut=sigmacut.asNumber(fb))
            if groupdDescriptions and gLabel in groupdDescriptions:
                uncoveredTopos.description = groupdDescriptions[gLabel]
            else:
                uncoveredTopos.description = gLabel
            # Fill the list with the SMS in topDict:
            uncoveredTopos.getSMSFrom(topDict)
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
    Holds information about a single coverage group: criteria for selecting and grouping SMS,
    function for reweighting cross sections, etc.
    """

    def __init__(self, label, smsFilter, reweightFactor,
                 smFinalStates, bsmFinalStates, sqrts, sigmacut=0.):
        """
        :param label: Group label
        :param smsFilter: Function which takes an SMS as argument and returns True (False) if
                             the SMS should (not) be selected.
        :param reweightFactor: Function which takes an SMS as argument and returns the reweighting
                              factor to be applied to the SMS weight.
        :param smFinalStates: List of Particle/MultiParticle objects used to group Z2-even particles appearing
                            in the final state
        :param bsmFinalStates: List of Particle/MultiParticle objects used to group Z2-odd particles appearing
                            in the final state
        :param sqrts: Value (with units) for the center of mass energy used to compute the missing cross sections.
                     If not specified the largest value available will be used.
        :param sigmacut: Minimum cross-section/weight value (after applying the reweight factor)
                       for an SMS to be included. The value should in fb (unitless)
        """

        self.finalStateSMS = []
        self.smFinalStates = smFinalStates
        self.bsmFinalStates = bsmFinalStates
        self.sqrts = sqrts
        self.sigmacut = sigmacut
        self.label = label
        self.smsFilter = smsFilter
        self.reweightFactor = reweightFactor

    def __str__(self):
        return self.label

    def __repr__(self):
        return str(self)

    def getSMSFrom(self, topDict):
        """
        Select the SMS from topDict according to self.smsFilter
        and build FinalStateSMS from the selected SMS.
        The FinalStateSMS weights corresponds to the missing cross-section
        with double counting from compressed SMS already accounted for.
        """

        # First select all SMS according to the filter (type of uncovered/missing topology):
        smsList = [sms for sms in topDict.getSMSList() if self.smsFilter(sms)]
        # Get missing xsections including the reweight factor:
        missingXandSMS = [(self.getMissingX(sms)*self.reweightFactor(sms), sms) for sms in smsList]
        # Only keep the SMS above sigmacut:
        missingXandSMS = [x for x in missingXandSMS[:] if x[0] > self.sigmacut]
        # Sort according to largest missingX, smallest size and largest ID
        missingXandSMS = sorted(missingXandSMS, key=lambda pt: (pt[0], -pt[1].canonName, -pt[1].smsID), reverse=True)

        # Split lists of SMS and missingX:
        missingXsecs = [pt[0] for pt in missingXandSMS]
        smsList = [pt[1] for pt in missingXandSMS]

        # Remove all SMS which are related to each other in order to avoid double counting
        # (keep always the first appearance in the list, so we always keep the ones with largest missing xsec)
        smsListUnique = []
        missingXsecsUnique = []
        ancestors = set()  # Keep track of all the ancestor ids of the SMS in the unique list
        for isms, sms in enumerate(smsList):
            # Get ancestor object ids
            # (safer than sms.smsID since some SMS can share smsID = 0 if they
            # were never directly inserted in the TopologyList)
            ancestorsIDs = set([id(sms)] + [id(sms) for sms in sms.getAncestors()])
            # If the SMS has any common ancestor with any of the previous SMS
            # or if it is an ancestor of any of the previous SMS
            # skip it to avoid double counting
            if ancestors.intersection(ancestorsIDs):
                continue
            smsListUnique.append(sms)
            missingXsecsUnique.append(missingXsecs[isms])
            ancestors = ancestors.union(ancestorsIDs)

        # Now that we only have unique SMS with their effective missing cross-sections
        # we create General SMS out of them
        for isms, sms in enumerate(smsListUnique):
            missingX = missingXsecsUnique[isms]
            if not missingX:
                continue
            self.addToFinalStateSMS(sms, missingX)

        # Finally sort general SMS by their missing cross-section:
        self.finalStateSMS = sorted(self.finalStateSMS, key=lambda fsEl: fsEl.missingX, reverse=True)

    def getMissingX(self, sms):
        """
        Calculate total missing cross section of an SMS, by recursively checking if its
        mothers already appear in the list.
        :param sms: SMS object

        :returns: missing cross section without units (in fb)
        """

        ancestorList = sms.getAncestors()
        alreadyChecked = []  # keep track of which SMS have already been checked
        # Get the (pre-loaded) total SMS weight in fb:
        missingX = sms._totalXsec
        if not missingX:
            return 0.
        overlapXsec = 0.
        for ancestor in ancestorList:  # check all ancestors (ancestorList is sorted by generation)
            # Skip entries which correspond to the SMS itself
            if ancestor is sms:
                continue
            # Make sure we do not subtract the same mother twice
            if any(ancestor is smsB for smsB in alreadyChecked):
                continue
            alreadyChecked.append(ancestor)
            # Check if ancestor passes the group filter (if it has been covered/tested or not):
            if self.smsFilter(ancestor):
                continue
            # Subtract the weight of the ancestor and skip checking for all the older family tree of the ancestor
            # (avoid subtracting the same weightList twice from mother and grandmother for instance)
            # (since the ancestorList is sorted by generation, the mother always
            # appears before the grandmother in the list)
            alreadyChecked += ancestor.getAncestors()
            if not hasattr(ancestor, '_totalXsec'):
                xsec = ancestor.weightList.getXsecsFor(self.sqrts)
                if not xsec:
                    ancestor._totalXsec = 0.0
                else:
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
        for fsEl in self.finalStateSMS:
            xsec += fsEl.missingX
        return xsec

    def addToFinalStateSMS(self, sms, missingX):
        """
        Adds an SMS to the list of missing topologies = final state SMS.
        If the SMS contributes to an SMS that is already in the list, add SMS and weight to
        the SMS.
        :parameter sms: SMS object to be added
        :parameter missingX: missing cross-section for the SMS (in fb)
        """

        newGenSMS = FinalStateSMS(sms, missingX, self.smFinalStates, self.bsmFinalStates)

        # Get index where to insert the new SMS
        # (if the SMS to the left of index is equal, add missing xsec)
        index = bisect(self.finalStateSMS, newGenSMS)
        if index != 0 and self.finalStateSMS[index-1] == newGenSMS:
            self.finalStateSMS[index-1]._contributingSMS.append(sms)
            self.finalStateSMS[index-1].missingX += missingX
        else:
            self.finalStateSMS.insert(index, newGenSMS)


class FinalStateSMS(TheorySMS):
    """
    This class represents a simplified SMS which
    only holds information about the final states. It holds a simple
    tree with one root (PV), having the final state nodes as its
    daughters.
    """

    def __init__(self, sms, missingX=None,
                smFinalStates=smDefault,
                bsmFinalStates=bsmDefault):

        TheorySMS.__init__(self)

        if missingX is None:
            missingX = sms.weightList.getMaxXsec()
        self.missingX = missingX

        # Store the primary mothers in dict, so
        # they will be replace by generic (anyBSM) particles
        nodeDict = {}
        for d in sms.daughterIndices(sms.rootIndex):
            # If mother is final state, skip it
            if sms.out_degree(d) == 0:
                continue
            daughter = sms.indexToNode(d)
            newNode = daughter.copy()
            newNode.particle = anyBSM
            nodeDict[d] = newNode

        compressedSMS = self.compressToFinalStates(sms)

        # Replace particles by inclusive particles, if possible:
        for nodeIndex in compressedSMS.nodeIndices:
            # Skip the primary mothers
            if nodeIndex in nodeDict:
                continue
            # Add root
            if nodeIndex == compressedSMS.rootIndex:
                nodeDict[nodeIndex] = compressedSMS.root
                continue
            node = compressedSMS.indexToNode(nodeIndex)
            ptc = node.particle
            newParticle = None
            if ptc.isSM:
                for smFS in smFinalStates:
                    if ptc == smFS:
                        newParticle = smFS
                        break
            else:
                for bsmFS in bsmFinalStates:
                    if ptc == bsmFS:
                        newParticle = bsmFS
                        break
            if newParticle is not None:
                newNode = node.copy()
                newNode.particle = newParticle
                nodeDict[nodeIndex] = newNode
            else:
                nodeDict[nodeIndex] = node

        self.copyTreeFrom(compressedSMS,nodeDict)

        self.sort(force=True)
        self.setGlobalProperties(weight=False)
        self._contributingSMS = [sms]


    def compressToFinalStates(self,sms):
        """
        Compress the SMS subtrees generated by the primary mothers.
        After the compression the SMS will contain the primary mothers with
        direct edges to its final states. Returns a compressed copy of sms.

        :param sms: SMS object
        """

        smsComp = sms.copy()
        for d in sms.daughterIndices(sms.rootIndex):
            # Skip final state mothers
            if sms.out_degree(d) == 0:
                continue
            # Remove all nodes for subtree
            # and collect the final state nodes
            fsNodes = []
            allDaughters = list(sms.dfsIndexIterator(d))
            fsNodes = [sms.indexToNode(n) for n in allDaughters
                       if sms.out_degree(n) == 0]
            # If the root daughters are all final states,
            # do nothing
            if len(fsNodes) == len(allDaughters):
                continue

            # Remove daughter nodes:
            smsComp.remove_nodes_from(allDaughters)
            # Add final state nodes with edges to root
            nIndices = smsComp.add_nodes_from(fsNodes)
            smsComp.add_edges_from(product([d],nIndices))

        # Re-compute canonical name and sort the tree
        smsComp._sorted = False
        smsComp._canonName = smsComp.computeCanonName()
        return smsComp

    def __str__(self):
        """
        Defines a slightly simplified version of the SMS string
        """

        fsStrs = []
        for d in self.daughterIndices(self.rootIndex):
            if self.out_degree(d) == 0:
                daughter = self.indexToNode(d)
                finalStates = f'({str(daughter)})'
            else:
                finalStates = str(tuple(self.daughters(d)[::-1]))
            fsStrs.append(finalStates.replace(' ',''))
        smsStr = 'PV > '
        smsStr += ', '.join(fsStrs)

        return smsStr


    def oldStr(self):
        """
        Generates a string using the old format (bracket notaion),
        if possible. For non Z2-like SMS, return the process string.

        :returns: string representation of the SMS (in bracket notation)
        """

        try:
            evenParticles, finalBSM, _ = self.treeToBrackets()
            smsStr = str(evenParticles).replace("'","").replace(" ","")
            smsStr += ' '+str(tuple(finalBSM)).replace("'","").replace(" ","")
            smsStr = smsStr.replace('~','')
        except:
            smsStr = str(self)

        return smsStr

