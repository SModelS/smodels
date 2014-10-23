#!/usr/bin/env python

"""
.. module:: theory.slhaDecomposer
   :synopsis: Decomposition of SLHA events and creation of TopologyLists.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import copy
import time
from smodels.theory import element, topology, crossSection
from smodels.theory.branch import Branch, decayBranches
from smodels.tools import modpyslha as pyslha
from smodels.tools.physicsUnits import fb, GeV
import smodels.particles
import logging

logger = logging.getLogger(__name__)


def decompose(slhafile, sigcut=0.1, doCompress=False, doInvisible=False,
              minmassgap=-1.*GeV, useXSecs=None):
    """
    Perform SLHA-based decomposition.
    :param sigcut: minimum sigma*BR to be generated, by default sigcut = 0.1 fb
    :param doCompress: turn mass compressed topologies on/off
    :param doInvisible: turn invisibly compressed topologies on/off
    :param minmassgap: maximum value for considering two R-odd particles
                       degenerate (only revelant for doCompress=True )
    :param useXSecs: optionally a dictionary with cross-sections for pair
                 production, by default reading the cross sections
                 from the SLHA file.
    :returns: TopologyList

    """
    t1 = time.time()

    if doCompress and minmassgap / GeV < 0.:
        logger.error("Asked for compression without specifying minmassgap. Please set minmassgap.")
        import sys
        sys.exit()

    if type(sigcut) == type(1.):
        sigcut = sigcut * fb

    # Get cross-section from file
    xSectionList = crossSection.getXsecFromSLHAFile(slhafile, useXSecs)
    # Get BRs and masses from file
    brDic, massDic = _getDictionariesFromSLHA(slhafile)
    # Only use the highest order cross-sections for each process
    xSectionList.removeLowerOrder()


    # Get maximum cross-sections (weights) for single particles (irrespective
    # of sqrtS)
    maxWeight = {}
    for pid in xSectionList.getPIDs():
        maxWeight[pid] = xSectionList.getXsecsFor(pid).getMaxXsec()

    # Create 1-particle branches with all possible mothers
    branchList = []
    for pid in maxWeight:
        branchList.append(Branch())
        branchList[-1].momID = pid
        branchList[-1].daughterID = pid
        branchList[-1].masses = [massDic[pid]]
        branchList[-1].maxWeight = maxWeight[pid]

    # Generate final branches (after all R-odd particles have decayed)
    finalBranchList = decayBranches(branchList, brDic, massDic, sigcut)

    smsTopList = topology.TopologyList()
    # Combine pairs of branches into elements according to production
    # cross-section list
    for pids in xSectionList.getPIDpairs():
        for branch1 in finalBranchList:
            for branch2 in finalBranchList:
                if branch1.momID == pids[0] and branch2.momID == pids[1]:
                    finalBR = branch1.maxWeight * branch2.maxWeight / \
                            (maxWeight[pids[0]] * maxWeight[pids[1]])
                    if type(finalBR) == type( 1. * fb):
                        finalBR = finalBR.asNumber()
                    weightList = xSectionList.getXsecsFor(pids) * finalBR

                    newElement = element.Element([branch1, branch2])
                    newElement.weight = weightList
                    allElements = [newElement]
                    # Perform compression
                    if doCompress or doInvisible:
                        allElements += newElement.compressElement(doCompress,
                                                                  doInvisible,
                                                                  minmassgap)

                    for el in allElements:
                        if el.weight.getMaxXsec() < sigcut:
                            # Skip elements with xsec below sigcut
                            continue
                        top = topology.Topology(el)
                        smsTopList.addList([top])
    logger.debug("slhaDecomposer done in " + str(time.time() - t1) + " s.")

    return smsTopList


def _getDictionariesFromSLHA(slhafile):
    """
    Create mass and BR dictionaries from an SLHA file.
    Ignore decay blocks with R-parity violating or unknown decays

    """

    res = pyslha.readSLHAFile(slhafile)

    rOdd = smodels.particles.rOdd.keys()
    rEven = smodels.particles.rEven.keys()
    knownParticles = rOdd + rEven

    # Get mass and branching ratios for all particles
    brDic = {}
    for pid in res.decays.keys():
        if pid in rEven:
            logger.warning("Ignoring %s decays",smodels.particles.rEven[pid])
            continue
        brs = []
        for decay in res.decays[pid].decays:
            nEven = nOdd = 0.
            for id in decay.ids:
                if id in rOdd: nOdd += 1
                elif id in rEven: nEven += 1
            if nOdd + nEven == len(decay.ids) and nOdd == 1:
                brs.append(decay)
            else:
                logger.warning("Ignoring decay: %i -> [%s]",pid,decay.ids)

        brsConj = copy.deepcopy(brs)
        for br in brsConj:
            br.ids = [-x for x in br.ids]
        brDic[pid] = brs
        brDic[-pid] = brsConj
    # Get mass list for all particles
    massDic = {}
    for pid in res.decays.keys():
        if pid and res.decays[pid].mass != None:
            massDic[pid] = abs(res.decays[pid].mass)* GeV
            massDic[-pid] = abs(res.decays[pid].mass) * GeV

    return brDic, massDic

if __name__ == "__main__":
    """
    Decomposes a given SLHA file.

    """
    import argparse
    import types
    argparser = argparse.ArgumentParser(description="Decomposes a given SLHA file.")
    argparser.add_argument('file', type=types.StringType, nargs=1,
                           help="SLHA file to decompose." )
    argparser.add_argument('-c', '--compress', action='store_true',
                           help="turn mass compressed topologies on." )
    argparser.add_argument('-i', '--invisible', action='store_true',
                           help="turn invisibly compressed topologies on." )
    argparser.add_argument('-x', '--crosssections', action='store_true',
                           help="print the production cross sections." )
    argparser.add_argument('-M', '--masses', action='store_true',
                           help="print the masses." )
    argparser.add_argument('-s', '--sigmacut', type=float, default=0.1,
                           help="minimum sigma*BR to be generated, in fb" )       
    argparser.add_argument('-m', '--minmassgap', type=float, default=-1.,
                           help="minimum sigma*BR to be generated, in fb" )       

    args = argparser.parse_args()
    File=args.file[0]
    def boolean(b):
        if b: return "true"
        return "false"
    print "Now trying to decompose:",File
    print "                       - compression=%s" % boolean(args.compress)
    print "                       - invisible=%s" % boolean(args.invisible)
    print "                       - sigmacut=%.2f" % args.sigmacut
    print "                       - minmassgap=%.2f" % args.minmassgap
    topolist=decompose ( File, args.sigmacut, args.compress, args.invisible,
            args.minmassgap * GeV )
    if len(topolist)>0:
        print "Found the following topologies: "
        for t in topolist:
            print
            print t
            for el in t.getElements():
                print el, ## ,el.getMasses()
                if args.masses:
                    print el.getMasses(),
                if args.crosssections:
                    print el.weight,
                print
