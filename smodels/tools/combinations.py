#!/usr/bin/env python3

"""
.. module:: combinations
   :synopsis: a module to contain the logic around combinations, be they
              SL-based or pyhf-based.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from smodels.tools.simplifiedLikelihoods import LikelihoodComputer, Data, UpperLimitComputer
from smodels.tools.smodelsLogging import logger
import numpy as np

def getCombinedUpperLimitFor(dataset, nsig, expected=False, deltas_rel=0.2):
    """
    Get combined upper limit. If covariances are given in globalInfo then
    simplified likelihood is used, else if json files are given pyhf
    cimbination is performed.

    :param nsig: list of signal events in each signal region/dataset. The list
                    should obey the ordering in globalInfo.datasetOrder.
    :param expected: return expected, not observed value
    :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.

    :returns: upper limit on sigma*eff
    """

    if dataset.type == "simplified":
        cov = dataset.globalInfo.covariance
        if type(cov) != list:
            raise SModelSError( "covariance field has wrong type: %s" % type(cov))
        if len(cov) < 1:
            raise SModelSError( "covariance matrix has length %d." % len(cov))

        computer = UpperLimitComputer(ntoys=10000)

        nobs = [x.dataInfo.observedN for x in dataset._datasets]
        bg = [x.dataInfo.expectedBG for x in dataset._datasets]
        no = nobs

        ret = computer.ulSigma(Data(observed=no, backgrounds=bg, covariance=cov,
                                    third_moment=None, nsignal=nsig, deltas_rel=deltas_rel),
                                    marginalize=dataset._marginalize,
                                    expected=expected)

        if ret != None:
            #Convert limit on total number of signal events to a limit on sigma*eff
            ret = ret/dataset.getLumi()
        logger.debug("SL upper limit : {}".format(ret))
        return ret
    elif dataset.type == "pyhf":
        logger.debug("Using pyhf")
        if all([s == 0 for s in nsig]):
            logger.warning("All signals are empty")
            return None
        ulcomputer = _getPyhfComputer( dataset, nsig )
        ret = ulcomputer.ulSigma(expected=expected)
        if ret == None:
            return None
        ret = ret/dataset.getLumi()
        logger.debug("pyhf upper limit : {}".format(ret))
        return ret
    else:
        logger.error ( "no covariance matrix or json file given in globalInfo.txt for %s" % dataset.globalInfo.id )
        raise SModelSError( "no covariance matrix or json file given in globalInfo.txt for %s" % dataset.globalInfo.id )

def computeCombinedLikelihood ( dataset, nsig, marginalize=False, deltas_rel=0.2 ):
    """ compute only lBSM
    :param nsig: predicted signal (list)
    :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    """
    if dataset.type == "pyhf":
        # Getting the path to the json files
        # Loading the jsonFiles
        ulcomputer = _getPyhfComputer( dataset, nsig, False )
        index = ulcomputer.getBestCombinationIndex()
        lbsm = ulcomputer.likelihood( index )
        return lbsm
    lbsm = _combinedLikelihood( dataset, nsig, marginalize, deltas_rel )
    return lbsm

def computeCombinedStatistics ( dataset, nsig, marginalize=False, deltas_rel=0.2 ):
    """ compute lBSM, lmax, and LSM in a single run 
    :param nsig: predicted signal (list)
    :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    """
    if dataset.type == "pyhf":
        # Getting the path to the json files
        # Loading the jsonFiles
        ulcomputer = _getPyhfComputer( dataset, nsig, False )
        index = ulcomputer.getBestCombinationIndex()
        lbsm = ulcomputer.likelihood( index )
        lmax = ulcomputer.lmax ( index )
        ulcomputer = _getPyhfComputer( dataset, [0.]*len(nsig), False )
        lsm = ulcomputer.likelihood ( index )
        return lbsm, lmax, lsm
    lbsm = _combinedLikelihood( dataset, nsig, marginalize, deltas_rel )
    lmax = _combinedLmax ( dataset, nsig, marginalize, deltas_rel )
    lsm =  _combinedLikelihood ( dataset, [0.]*len(nsig), marginalize, deltas_rel )
    return lbsm, lmax, lsm

def _getPyhfComputer ( dataset, nsig, normalize = True ):
    """ create the pyhf ul computer object
    :param normalize: if true, normalize nsig
    :returns: pyhf upper limit computer, and combinations of signal regions
    """
    # Getting the path to the json files
    jsonFiles = [js for js in dataset.globalInfo.jsonFiles]
    jsons = dataset.globalInfo.jsons.copy()
    datasets = [ds.getID() for ds in dataset._datasets]
    total = sum(nsig)
    if total == 0.: # all signals zero? can divide by anything!
        total = 1.
    if normalize:
        nsig = [s/total for s in nsig] # Normalising signals to get an upper limit on the events count
    # Filtering the json files by looking at the available datasets
    for jsName in dataset.globalInfo.jsonFiles:
        if all([ds not in dataset.globalInfo.jsonFiles[jsName] for ds in datasets]):
            # No datasets found for this json combination
            jsIndex = jsonFiles.index(jsName)
            jsonFiles.pop(jsIndex)
            jsons.pop(jsIndex)
            continue
        if not all([ds in datasets for ds in dataset.globalInfo.jsonFiles[jsName]]):
            # Some SRs are missing for this json combination
            logger.error("Wrong json definition in globalInfo.jsonFiles for json : %s" % jsName)
    logger.debug("list of datasets: {}".format(datasets))
    logger.debug("jsonFiles after filtering: {}".format(jsonFiles))
    # Constructing the list of signals with subsignals matching each json
    nsignals = list()
    for jsName in jsonFiles:
        subSig = list()
        for srName in dataset.globalInfo.jsonFiles[jsName]:
            try:
                index = datasets.index(srName)
            except ValueError:
                logger.error("%s signal region provided in globalInfo is not in the list of datasets" % srName)
            sig = nsig[index]
            subSig.append(sig)
        nsignals.append(subSig)
    # Loading the jsonFiles, unless we already have them (because we pickled)
    from smodels.tools.pyhfInterface import PyhfData, PyhfUpperLimitComputer
    data = PyhfData(nsignals, jsons, jsonFiles )
    if data.errorFlag: return None
    ulcomputer = PyhfUpperLimitComputer(data)
    return ulcomputer

def _combinedLikelihood( dataset, nsig, marginalize=False, deltas_rel=0.2 ):
    """
    Computes the (combined) likelihood to observe nobs events, given a
    predicted signal "nsig", with nsig being a vector with one entry per
    dataset.  nsig has to obey the datasetOrder. Deltas is the error on
    the signal.
    :param nsig: predicted signal (list)
    :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    :returns: likelihood to observe nobs events (float)
    """

    if dataset.type == "simplified":
        if len(dataset._datasets) == 1:
            if isinstance(nsig,list):
                nsig = nsig[0]
            return dataset._datasets[0].likelihood(nsig,marginalize=marginalize)
        nobs = [ x.dataInfo.observedN for x in dataset._datasets]
        bg = [ x.dataInfo.expectedBG for x in dataset._datasets]
        cov = dataset.globalInfo.covariance
        computer = LikelihoodComputer(Data(nobs, bg, cov, None, nsig,
                                           deltas_rel=deltas_rel))
        return computer.likelihood(nsig, marginalize=marginalize )
    elif dataset.type == "pyhf":
        # Getting the path to the json files
        # Loading the jsonFiles
        ulcomputer = _getPyhfComputer( dataset, nsig, False )
        return ulcomputer.likelihood()
    else:
        logger.error("Asked for combined likelihood, but no covariance or json file given: %s" % dataset.type )
        return None

def _combinedLmax ( dataset, nsig, marginalize, deltas_rel, nll=False, expected=False,
           allowNegativeSignals = False ):
    """ compute likelihood at maximum """
    if dataset.type == "simplified":
        nobs = [ x.dataInfo.observedN for x in dataset._datasets]
        if expected:
            # nobs = [ x.dataInfo.expectedBG for x in dataset._datasets]
            nobs = [ int(np.round(x.dataInfo.expectedBG)) for x in dataset._datasets]
        bg = [ x.dataInfo.expectedBG for x in dataset._datasets]
        cov = dataset.globalInfo.covariance
        if type(nsig) in [ list, tuple ]:
            nsig = np.array(nsig)
        computer = LikelihoodComputer(Data(nobs, bg, cov, None, nsig, deltas_rel=deltas_rel))
        mu_hat = computer.findMuHat ( nsig, allowNegativeSignals=allowNegativeSignals )
        musig = nsig * mu_hat
        return computer.likelihood ( musig, marginalize=marginalize, nll=nll )
    if dataset.type == "pyhf":
        ulcomputer = _getPyhfComputer( dataset, nsig, False )
        return ulcomputer.lmax ( nll=nll )
    return -1.

