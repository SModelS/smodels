"""
.. module:: datasetObj
   :synopsis: Holds the classes and methods used to read and store the information in the
              data folders.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import os
import glob
import numpy as np
from smodels.experiment import txnameObj, infoObj
from smodels.tools.physicsUnits import fb
from smodels.tools.simplifiedLikelihoods import LikelihoodComputer, Data, UpperLimitComputer
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.theory.auxiliaryFunctions import getAttributesFrom, getValuesForObj
from smodels.tools.smodelsLogging import logger
from smodels.theory.auxiliaryFunctions import elementsInStr
from smodels.theory.element import Element

from typing import Text, Union, List, Dict, Optional
import itertools
from spey import get_uncorrelated_region_statistical_model, get_multi_region_statistical_model, ExpectationType

# if on, will check for overlapping constraints
_complainAboutOverlappingConstraints = True


class DataSet(object):
    """
    Holds the information to a data set folder (TxName objects, dataInfo,...)
    """

    def __init__(self, path=None, info=None, createInfo=True,
                 discard_zeroes=True, databaseParticles=None):
        """
        :param discard_zeroes: discard txnames with zero-only results
        """

        self.path = path
        self.globalInfo = info
        self.txnameList = []

        if path and createInfo:
            logger.debug('Creating object based on data folder : %s' % self.path)

            # Get data folder info:
            if not os.path.isfile(os.path.join(path, "dataInfo.txt")):
                logger.error("dataInfo.txt file not found in " + path)
                raise TypeError
            self.dataInfo = infoObj.Info(os.path.join(path, "dataInfo.txt"))

            # Get list of TxName objects:
            for txtfile in glob.iglob(os.path.join(path, "*.txt")):
                try:
                    txname = txnameObj.TxName(txtfile, self.globalInfo,
                                              self.dataInfo, databaseParticles)
                    if discard_zeroes and txname.hasOnlyZeroes():
                        logger.debug("%s, %s has only zeroes. discard it." %
                                     (self.path, txname.txName))
                        continue
                    self.txnameList.append(txname)
                except TypeError:
                    continue

            self.txnameList.sort()
            self.checkForRedundancy(databaseParticles)

    def isCombinableWith(self, other):
        """
        Function that reports if two datasets are mutually uncorrelated = combinable.

        :param other: datasetObj to compare self with
        """
        id1, id2 = self.globalInfo.id, other.globalInfo.id
        if id1 == id2:  # we are always correlated with ourselves
            return False
        from smodels.tools.physicsUnits import TeV
        ds = abs(self.globalInfo.sqrts.asNumber(TeV) - other.globalInfo.sqrts.asNumber(TeV))
        if ds > 1e-5:  # not the same
            return True

        def getCollaboration(ds):
            return "CMS" if "CMS" in ds.globalInfo.id else "ATLAS"
        coll1, coll2 = getCollaboration(self), getCollaboration(other)
        if coll1 != coll2:
            return True

        if self.isGlobalFieldCombinableWith_(other):
            return True
        if other.isGlobalFieldCombinableWith_(self):
            return True
        if self.isLocalFieldCombinableWith_(other):
            return True
        if other.isLocalFieldCombinableWith_(self):
            return True
        if self.isCombMatrixCombinableWith_(other):
            return True
        return False

    def isCombMatrixCombinableWith_(self, other):
        """ check for combinability via the combinations matrix """
        if not hasattr(self.globalInfo, "_combinationsmatrix"):
            return False
        if self.globalInfo._combinationsmatrix is None:
            return False
        idSelf = self.globalInfo.id
        didSelf = self.dataInfo.dataId
        selflabel = f"{idSelf}:{didSelf}"
        idOther = other.globalInfo.id
        didOther = other.dataInfo.dataId
        otherlabel = f"{idOther}:{didOther}"
        for label, combs in self.globalInfo._combinationsmatrix.items():
            if label in [idSelf, didSelf]:
                # match! with self! is "other" in combs?
                if idOther in combs or otherlabel in combs:
                    return True
            if label in [idOther, didOther]:
                # match! with other! is "self" in combs?
                if idSelf in combs or selflabel in combs:
                    return True
        return False

    def isGlobalFieldCombinableWith_(self, other):
        """ check for 'combinableWith' fields in globalInfo, check if <other> matches.
        this check is at analysis level (not at dataset level).

        :params other: a dataset to check against
        :returns: true, if pair is marked as combinable, else false
        """
        if not hasattr(self.globalInfo, "combinableWith"):
            return False
        tokens = self.globalInfo.combinableWith.split(",")
        idOther = other.globalInfo.id
        for t in tokens:
            if ":" in t:
                logger.error("combinableWith field in globalInfo is at the analysis level. You specified a dataset-level combination %s." % t)
                raise SModelSError()
        if idOther in tokens:
            return True
        return False

    def isLocalFieldCombinableWith_(self, other):
        """ check for 'combinableWith' fields in globalInfo, check if <other> matches.
        this check is at dataset level (not at dataset level).

        :params other: a dataset to check against
        :returns: true, if pair is marked as combinable, else false
        """
        if not hasattr(self.dataInfo, "combinableWith"):
            return False
        tokens = self.dataInfo.combinableWith.split(",")
        for t in tokens:
            if ":" not in t:
                logger.error("combinableWith field in dataInfo is at the dataset level. You specified an analysis-level combination %s." % t)
                raise SModelSError()
        idOther = other.globalInfo.id
        didOther = other.dataInfo.dataId
        label = f"{idOther}:{didOther}"
        if label in tokens:
            return True
        return False

    def checkForRedundancy(self, databaseParticles):
        """ In case of efficiency maps, check if any txnames have overlapping
            constraints. This would result in double counting, so we dont
            allow it. """
        if self.getType() == "upperLimit":
            return False
        logger.debug("checking for redundancy")
        datasetElements = []
        for tx in self.txnameList:
            if hasattr(tx, 'finalState'):
                finalState = tx.finalState
            else:
                finalState = ['MET', 'MET']
            if hasattr(tx, 'intermediateState'):
                intermediateState = tx.intermediateState
            else:
                intermediateState = None
            for el in elementsInStr(str(tx.constraint)):
                newEl = Element(el, finalState, intermediateState,
                                model=databaseParticles)
                datasetElements.append(newEl)
        combos = itertools.combinations(datasetElements, 2)
        for x, y in combos:
            if x == y and _complainAboutOverlappingConstraints:
                errmsg = "Constraints (%s) and (%s) appearing in dataset %s:%s overlap "\
                        "(may result in double counting)." % \
                        (x, y, self.getID(), self.globalInfo.id)
                logger.error(errmsg)
                raise SModelSError(errmsg)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        if self.dataInfo.dataId:
            ret = "Dataset %s: %s" % (self.dataInfo.dataId, ", ".join(map(str, self.txnameList)))
        else:
            ret = "Dataset: %s" % (", ".join(map(str, self.txnameList)))
        return ret

    def __repr__(self):
        if self.dataInfo.dataId:
            return self.dataInfo.dataId
        else:
            return 'Dataset'

    def __eq__(self, other):
        if type(other) != type(self):
            return False
        if self.dataInfo != other.dataInfo:
            return False
        if len(self.txnameList) != len(other.txnameList):
            return False
        return True

    def getType(self):
        """
        Return the dataset type (EM/UL)
        """

        return self.dataInfo.dataType

    def getID(self):
        """
        Return the dataset ID
        """

        return self.dataInfo.dataId

    def getLumi(self):
        """
        Return the dataset luminosity. If not defined for the dataset, use
        the value defined in globalInfo.lumi.
        """

        if hasattr(self, 'lumi'):
            return self.lumi
        else:
            return self.globalInfo.lumi

    def getTxName(self, txname):
        """
        get one specific txName object.
        """
        for tn in self.txnameList:
            if tn.txName == txname:
                return tn
        return None

    def getEfficiencyFor(self, txname, mass):
        """
        Convenience function. Get efficiency for mass
        assuming no lifetime rescaling. Same as self.getTxName(txname).getEfficiencyFor(m)
        """
        txname = self.getTxName(txname)
        if txname:
            return txname.getEfficiencyFor(mass)
        else:
            return None

    def getValuesFor(self, attribute):
        """
        Returns a list for the possible values appearing in the ExpResult
        for the required attribute (sqrts,id,constraint,...).
        If there is a single value, returns the value itself.

        :param attribute: name of a field in the database (string).
        :return: list of unique values for the attribute
        """

        return getValuesForObj(self, attribute)

    def likelihood(self, nsig, deltas_rel=0.2, marginalize=False, expected=False, backend="simplified_likelihoods"):
        """
        Computes the likelihood to observe nobs events,
        given a predicted signal "nsig", assuming "deltas_rel"
        error on the signal efficiency.
        The values observedN, expectedBG, and bgError are part of dataInfo.

        :param nsig: predicted signal (float)

        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param expected: if true, compute prior expected, if false compute observed
                         if "posteriori" compute posterior expected
        :param backend: the backend used to compute the likelihood. Can be "pyhf" or "SL".
        :returns: likelihood to observe nobs events (float)
        """
        if deltas_rel != 0.2:
            logger.warning("Relative uncertainty on signal not supported by spey for a single region.")

        if backend == "pyhf":
            if marginalize == True:
                logger.error("pyhf backend cannot marginalize.")
            args={}
        elif backend == "simplified_likelihoods":
            args={"marginalize":marginalize}
        else:
            logger.error('%s is not a valid backend. Possible backends are "pyhf" and "simplified_likelihoods".' %backend)
            return None

        statModel = get_uncorrelated_region_statistical_model(observations=float(self.dataInfo.observedN),
                                                                backgrounds=float(self.dataInfo.expectedBG),
                                                                background_uncertainty=float(self.dataInfo.bgError),
                                                                signal_yields=float(nsig),
                                                                xsection=nsig/self.getLumi(),
                                                                analysis=self.globalInfo.id,
                                                                backend=backend
                                                            )
        expectedDict = {False:ExpectationType.observed,
                        True:ExpectationType.apriori,
                        "posteriori":ExpectationType.aposteriori}
        if expected not in expectedDict.keys():
            logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
            return None

        ret = statModel.likelihood(poi_test=1., expected=expectedDict[expected], return_nll=False, **args)
        return ret

    # def likelihood(self, nsig, deltas_rel=0.2, marginalize=False, expected=False):
    #     """
    #     Computes the likelihood to observe nobs events,
    #     given a predicted signal "nsig", assuming "deltas_rel"
    #     error on the signal efficiency.
    #     The values observedN, expectedBG, and bgError are part of dataInfo.
    #
    #     :param nsig: predicted signal (float)
    #
    #     :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    #     :param marginalize: if true, marginalize nuisances. Else, profile them.
    #     :param expected: if true, compute prior expected, if false compute observed
    #                      if "posteriori" compute posterior expected
    #     :returns: likelihood to observe nobs events (float)
    #     """
    #     obs = self.dataInfo.observedN
    #     if expected:  # this step is done for prior and posterior expected
    #         obs = self.dataInfo.expectedBG
    #
    #     m = Data(obs, self.dataInfo.expectedBG, self.dataInfo.bgError**2,
    #              deltas_rel=deltas_rel)
    #     computer = LikelihoodComputer(m)
    #     if expected == "posteriori":
    #         thetahat = computer.findThetaHat(0.)
    #         if type(self.dataInfo.expectedBG) in [float, np.float64, int]:
    #             thetahat = float(thetahat[0])
    #
    #         obs = self.dataInfo.expectedBG + thetahat
    #         m = Data(obs, self.dataInfo.expectedBG, self.dataInfo.bgError**2,
    #                  deltas_rel=deltas_rel)
    #         # if abs ( nsig[0]-1 ) < 1e-5:
    #         #    print ( f"COMB ebg={self.dataInfo.expectedBG:.3f} obs={obs:.3f} nsig {nsig[0]:.3f}" )
    #         computer = LikelihoodComputer(m)
    #     ret = computer.likelihood(nsig, marginalize=marginalize)
    #     if hasattr ( computer, "theta_hat" ):
    #         ## seems like someone wants to debug them
    #         self.theta_hat = computer.theta_hat
    #
    #     return ret

    def lmax(self, nsig, deltas_rel=0.2, marginalize=False, expected=False, allowNegativeSignals=False, backend="simplified_likelihoods"):
        """
        Convenience function, computes the likelihood at nsig = observedN - expectedBG,
        assuming "deltas_rel" error on the signal efficiency.
        The values observedN, expectedBG, and bgError are part of dataInfo.

        :param nsig: predicted signal (float)
        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param expected: Compute expected instead of observed likelihood
        :param allowNegativeSignals: if False, then negative nsigs are replaced with 0.
        :param backend: the backend used to compute the likelihood. Can be "pyhf" or "SL".

        :returns: likelihood to observe nobs events (float)
        """
        if deltas_rel != 0.2:
            logger.warning("Relative uncertainty on signal not supported by spey for a single region.")

        if backend == "pyhf":
            if marginalize == True:
                logger.error("pyhf backend cannot marginalize.")
            args={"iteration_threshold":3} #default in spey
        elif backend == "simplified_likelihoods":
            args={}
        else:
            logger.error('%s is not a valid backend. Possible backends are "pyhf" and "simplified_likelihoods".' %backend)
            return None

        statModel = get_uncorrelated_region_statistical_model(observations=float(self.dataInfo.observedN),
                                                                backgrounds=float(self.dataInfo.expectedBG),
                                                                background_uncertainty=float(self.dataInfo.bgError),
                                                                signal_yields=float(nsig),
                                                                xsection=nsig/self.getLumi(),
                                                                analysis=self.globalInfo.id,
                                                                backend=backend
                                                            )

        expectedDict = {False:ExpectationType.observed,
                        True:ExpectationType.apriori,
                        "posteriori":ExpectationType.aposteriori}
        if expected not in expectedDict.keys():
            logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
            return None

        muhat, maxLL = statModel.maximize_likelihood(return_nll=False, expected=expectedDict[expected], allow_negative_signal=allowNegativeSignals , **args)

        self.muhat = muhat
        self.sigma_mu = np.sqrt(self.dataInfo.observedN + self.dataInfo.bgError**2) / nsig

        return maxLL

    # def lmax(self, deltas_rel=0.2, marginalize=False, expected=False,allowNegativeSignals=False):
    #     """
    #     Convenience function, computes the likelihood at nsig = observedN - expectedBG,
    #     assuming "deltas_rel" error on the signal efficiency.
    #     The values observedN, expectedBG, and bgError are part of dataInfo.
    #
    #     :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    #     :param marginalize: if true, marginalize nuisances. Else, profile them.
    #     :param expected: Compute expected instead of observed likelihood
    #     :param allowNegativeSignals: if False, then negative nsigs are replaced with 0.
    #
    #     :returns: likelihood to observe nobs events (float)
    #     """
    #     obs = self.dataInfo.observedN
    #     if expected:
    #         obs = self.dataInfo.expectedBG
    #         if expected == "posteriori":
    #             m = Data(obs, self.dataInfo.expectedBG, self.dataInfo.bgError**2,
    #                      deltas_rel=deltas_rel)
    #             computer = LikelihoodComputer(m)
    #             thetahat = computer.findThetaHat(0.)
    #             if type(self.dataInfo.expectedBG) in [float, np.float64,
    #                                             np.float32, int, np.int64, np.int32]:
    #                 thetahat = float(thetahat[0])
    #             obs = self.dataInfo.expectedBG + thetahat
    #
    #     m = Data(obs, self.dataInfo.expectedBG, self.dataInfo.bgError**2,
    #              deltas_rel=deltas_rel)
    #     computer = LikelihoodComputer(m)
    #     ret = computer.lmax ( marginalize=marginalize, nll=False,
    #                           allowNegativeSignals=allowNegativeSignals )
    #     if hasattr ( computer, "theta_hat" ):
    #         ## seems like someone wants to debug them
    #         self.theta_hat = computer.theta_hat
    #     if hasattr ( computer, "muhat" ):
    #         ## seems like someone wants to debug them
    #         self.muhat = computer.muhat
    #     if hasattr ( computer, "sigma_mu" ):
    #         ## seems like someone wants to debug them
    #         self.sigma_mu = computer.sigma_mu
    #     return ret

    def chi2(self, nsig, deltas_rel=0.2, marginalize=False, backend="simplified_likelihoods"):
        """
        Computes the chi2 for a given number of observed events "nobs",
        given number of signal events "nsig", and error on signal "deltas".

        nobs, expectedBG and bgError are part of dataInfo.
        :param nsig: predicted signal (float)
        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param backend: the backend used to compute the likelihood. Can be "pyhf" or "SL".
        :return: chi2 (float)
        """

        if deltas_rel != 0.2:
            logger.warning("Relative uncertainty on signal not supported by spey for a single region.")

        if backend == "pyhf":
            if marginalize == True:
                logger.error("pyhf backend cannot marginalize.")
            args={"iteration_threshold":3} #default in spey
        elif backend == "simplified_likelihoods":
            args={"marginalize":marginalize}
        else:
            logger.error('%s is not a valid backend. Possible backends are "pyhf" and "simplified_likelihoods".' %backend)
            return None

        statModel = get_uncorrelated_region_statistical_model(observations=float(self.dataInfo.observedN),
                                                                backgrounds=float(self.dataInfo.expectedBG),
                                                                background_uncertainty=float(self.dataInfo.bgError),
                                                                signal_yields=float(nsig),
                                                                xsection=nsig/self.getLumi(),
                                                                analysis=self.globalInfo.id,
                                                                backend=backend
                                                            )
        ret = statModel.chi2(poi_test=1., expected = ExpectationType.observed, allow_negative_signal=True)

        return ret

    # def chi2(self, nsig, deltas_rel=0.2, marginalize=False):
    #     """
    #     Computes the chi2 for a given number of observed events "nobs",
    #     given number of signal events "nsig", and error on signal "deltas".
    #
    #     nobs, expectedBG and bgError are part of dataInfo.
    #     :param nsig: predicted signal (float)
    #     :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    #     :param marginalize: if true, marginalize nuisances. Else, profile them.
    #     :return: chi2 (float)
    #     """
    #
    #     m = Data(self.dataInfo.observedN, self.dataInfo.expectedBG,
    #              self.dataInfo.bgError**2, deltas_rel=deltas_rel)
    #     computer = LikelihoodComputer(m)
    #     ret = computer.chi2(nsig, marginalize=marginalize)
    #
    #     return ret

    def folderName(self):
        """
        Name of the folder in text database.
        """
        return os.path.basename(self.path)

    def getAttributes(self, showPrivate=False):
        """
        Checks for all the fields/attributes it contains as well as the
        attributes of its objects if they belong to smodels.experiment.

        :param showPrivate: if True, also returns the protected fields (_field)
        :return: list of field names (strings)

        """

        attributes = getAttributesFrom(self)

        if not showPrivate:
            attributes = list(filter(lambda a: a[0] != '_', attributes))

        return attributes

    def getUpperLimitFor(self, element=None, expected=False, txnames=None, compute=False, alpha=0.05, deltas_rel=0.2):
        """
        Returns the upper limit for a given element (or mass) and txname. If
        the dataset hold an EM map result the upper limit is independent of
        the input txname or mass.
        For UL results if an Element object is given the corresponding upper limit
        will be rescaled according to the lifetimes of the element intermediate particles.
        On the other hand, if a mass is given, no rescaling will be applied.

        :param txname: TxName object or txname string (only for UL-type results)
        :param element: Element object or mass array with units (only for UL-type results)
        :param alpha: Can be used to change the C.L. value. The default value is 0.05
                      (= 95% C.L.) (only for  efficiency-map results)
        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
        :param expected: Compute expected limit, i.e. Nobserved = NexpectedBG
                         (only for efficiency-map results)
        :param compute: If True, the upper limit will be computed
                        from expected and observed number of events.
                        If False, the value listed in the database will be used
                        instead.
        :return: upper limit (Unum object)
        """

        if self.getType() == 'efficiencyMap':
            upperLimit = self.getSRUpperLimit(expected=expected)
            if type(upperLimit) == type(None):
                return None
            if (upperLimit/fb).normalize()._unit:
                logger.error("Upper limit defined with wrong units for %s and %s"
                             % (self.globalInfo.id, self.getID()))
                return False
            else:
                return upperLimit

        elif self.getType() == 'upperLimit':
            if not txnames or not element:
                logger.error("A TxName and mass array must be defined when \
                             computing ULs for upper-limit results.")
                return False
            elif isinstance(txnames, list):
                if len(txnames) != 1:
                    logger.error("txnames must be a TxName object, a string or a list with a single Txname object")
                    return False
                else:
                    txname = txnames[0]
            else:
                txname = txnames

            if not isinstance(txname, txnameObj.TxName) and \
                    not isinstance(txname, str):
                logger.error("txname must be a TxName object or a string")
                return False

            if not isinstance(element, list) and not isinstance(element, Element):
                logger.error("Element must be an element object or a mass array")
                return False

            for tx in self.txnameList:
                if tx == txname or tx.txName == txname:
                    upperLimit = tx.getULFor(element, expected)

            return upperLimit

        else:
            logger.warning("Unkown data type: %s. Data will be ignored.",
                           self.getType())
            return None

    def getSRUpperLimit(self,expected=False):
        """
        Returns the 95% upper limit on the signal*efficiency for a given dataset (signal region).
        Only to be used for efficiency map type results.

        :param expected: If True, return the expected limit ( i.e. Nobserved = NexpectedBG )

        :return: upper limit value
        """

        if not self.getType() == 'efficiencyMap':
            logger.error("getSRUpperLimit can only be used for efficiency map results!")
            raise SModelSError()

        if expected:
            if hasattr(self.dataInfo, "upperLimit") and not hasattr(self.dataInfo, "expectedUpperLimit"):
                logger.info("expectedUpperLimit field not found. Returning None instead.")
                return None

            if hasattr(self.dataInfo, "expectedUpperLimit"):
                return self.dataInfo.expectedUpperLimit
        else:
            if hasattr(self.dataInfo, "upperLimit"):
                return self.dataInfo.upperLimit

    # !TP
    def getStatModel(self,
        signal: Union[float, np.ndarray],
        backend: Text = "simplified_likelihoods"
    ):
        """
        Create statistical model from a single bin or multiple uncorrelated regions

        :param signal: signal yields
        :param backend: "pyhf" or "simplified_likelihoods"

        :return: spey StatisticalModel object

        :raises NotImplementedError: If requested backend has not been recognised.
        """

        from spey import get_uncorrelated_region_statistical_model

        if hasattr(self, "statModel"):
            return self.statModel
        else:
            self.statModel =  get_uncorrelated_region_statistical_model(observations = float(self.dataInfo.observedN),
                                                                        backgrounds = float(self.dataInfo.expectedBG),
                                                                        background_uncertainty = float(self.dataInfo.bgError),
                                                                        signal_yields = signal,
                                                                        xsection = nsig/self.getLumi(),
                                                                        analysis = self.globalInfo.id,
                                                                        backend = backend
                                                                        )
        return self.statModel



class CombinedDataSet(object):
    """
    Holds the information for a combined dataset (used for combining multiple datasets).
    """

    def __init__(self, expResult):

        self.path = expResult.path
        self.globalInfo = expResult.globalInfo
        self._datasets = expResult.datasets[:]
        self.origdatasets = expResult.origdatasets[:]
        self._marginalize = False
        self.sortDataSets()
        self.findType()

    def findType(self):
        """ find the type of the combined dataset """
        self.type = "bestSR"  # type of combined dataset, e.g. pyhf, or simplified
        if hasattr(self.globalInfo, "covariance"):
            self.type = "simplified"
        if hasattr(self.globalInfo, "jsonFiles"):
            self.type = "pyhf"

    def __str__(self):
        ret = "Combined Dataset (%i datasets)" % len(self._datasets)
        return ret

    def getIndex(self, dId, datasetOrder):
        """ get the index of dataset within the datasetOrder,
            but allow for abbreviated names
        :param dId: id of dataset to search for, may be abbreviated
        :param datasetOrder: the ordered list of datasetIds, long form
        :returns: index, or -1 if not found
        """
        if dId in datasetOrder:
            # easy peasy, we found the dId
            return datasetOrder.index(dId)
        foundIndex = -1
        for i, ds in enumerate(datasetOrder):
            if ds.endswith(":" + dId):
                # ok, dId is the abbreviated form
                if foundIndex == -1:
                    foundIndex = i
                else:
                    line = f"abbreviation {dId} matches more than one id in datasetOrder"
                    logger.error(line)
                    raise SModelSError(line)
        return foundIndex

    def sortDataSets(self):
        """
        Sort datasets according to globalInfo.datasetOrder.
        """
        if hasattr(self.globalInfo, "covariance"):
            datasets = self.origdatasets[:]
            if not hasattr(self.globalInfo, "datasetOrder"):
                raise SModelSError("datasetOrder not given in globalInfo.txt for %s" % self.globalInfo.id)
            datasetOrder = self.globalInfo.datasetOrder
            if isinstance(datasetOrder, str):
                datasetOrder = [datasetOrder]

            if len(datasetOrder) != len(datasets):
                raise SModelSError( f"Number of datasets in the datasetOrder field {len(datasetOrder)} does not match the number of datasets {len(datasets)}/{len(self.origdatasets)} for {self.globalInfo.id}" )
            ## need to reinitialise, we might have lost some datasets when filtering
            self._datasets = [ None ] * len(datasets)
            for dataset in datasets:
                idx = self.getIndex(dataset.getID(), datasetOrder)
                if idx == -1:
                    raise SModelSError("Dataset ID %s not found in datasetOrder" % dataset.getID())
                self._datasets[idx] = dataset
                # dsIndex = datasetOrder.index(dataset.getID())
                # self._datasets[dsIndex] = dataset

    def getType(self):
        """
        Return the dataset type (combined)
        """

        return 'combined'

    def getID(self):
        """
        Return the ID for the combined dataset
        """

        return '(combined)'

    def getLumi(self):
        """
        Return the dataset luminosity. For CombinedDataSet always return
        the value defined in globalInfo.lumi.
        """

        return self.globalInfo.lumi

    def getDataSet(self, datasetID):
        """
        Returns the dataset with the corresponding dataset ID.
        If the dataset is not found, returns None.

        :param datasetID: dataset ID (string)

        :return: DataSet object if found, otherwise None.
        """

        for dataset in self._datasets:
            if datasetID == dataset.getID():
                return dataset

        return None

    # !TP
    def getWSInfo(self, jsons):
        """
        Getting informations from the json files

        :param jsons: list of json instances.

        :return: wsInfo list of dictionaries (one dictionary for each json file) containing useful information about the json files.
            - :key signalRegions: list of dictonaries with 'json path' and 'size' (number of bins) of the 'signal regions' channels in the json files
            - :key otherRegions: list of strings indicating the path to the control and validation region channels
        """
        # Identifying the path to the SR and VR channels in the main workspace files
        wsInfo = []  # workspace specifications
        if not isinstance(jsons, list):
            logger.error("The `jsons` parameter must be of type list")
            return
        for ws in jsons:
            wsChannelsInfo = {}
            wsChannelsInfo["signalRegions"] = []
            wsChannelsInfo["otherRegions"] = []
            if not "channels" in ws.keys():
                logger.error(
                    "Json file number {} is corrupted (channels are missing)".format(
                        jsons.index(ws)
                    )
                )
                wsInfo = None
                return
            for i_ch, ch in enumerate(ws["channels"]):
                if ch["name"][:2] == "SR":  # if channel name starts with 'SR'
                    wsChannelsInfo["signalRegions"].append(
                        {
                            "path": "/channels/"
                            + str(i_ch)
                            + "/samples/0",  # Path of the new sample to add (signal prediction)
                            "size": len(ch["samples"][0]["data"]),
                        }
                    )  # Number of bins
                else:
                    wsChannelsInfo["otherRegions"].append("/channels/" + str(i_ch))
            wsChannelsInfo["otherRegions"].sort(
                key=lambda path: path.split("/")[-1], reverse=True
            )  # Need to sort correctly the paths to the channels to be removed
            wsInfo.append(wsChannelsInfo)
        return wsInfo

    # !TP
    def patchMaker(self, jsons, wsInfo, nsignals, includeCRs):
        """
        Method that creates the list of patches to be applied to the `jsons` workspaces, one for each region given the `nsignals` and the informations available in `wsInfo` and the content of the `jsons`

        :param jsons: list of json instances.
        :param wsInfo: list of dictionaries (one dictionary for each json file) containing useful information about the json files
        :param nsignals: list of list of signal yields (one list for each json file).
        :param includeCRs: if True leaves the pacth unchanged,
                           if False adds to the patch an operation that removes the CRs from the json files.

        NB: It seems we need to include the change of the "modifiers" in the patches as well

        :return: the list of patches, one for each workspace
        """
        if wsInfo == None:
            return None
        # Constructing the patches to be applied on the main workspace files
        patches = []
        for ws, info, subSig in zip(jsons, wsInfo, nsignals):
            patch = []
            for srInfo in info["signalRegions"]:
                nBins = srInfo["size"]
                operator = {}
                operator["op"] = "add"
                operator["path"] = srInfo["path"]
                value = {}
                value["data"] = subSig[:nBins]
                subSig = subSig[nBins:]
                value["modifiers"] = []
                value["modifiers"].append({"data": None, "type": "normfactor", "name": "mu_SIG"})
                value["modifiers"].append({"data": None, "type": "lumi", "name": "lumi"})
                value["name"] = "bsm"
                operator["value"] = value
                patch.append(operator)
            if includeCRs:
                logger.debug("keeping the CRs")
            else:
                for path in info["otherRegions"]:
                    patch.append({"op": "remove", "path": path})
            patches.append(patch)
        return patches

    # !TP
    def _getPatches(self, nsig):
        """
        :param nsig: list of signal yields (not relative).
        :param normalize: if true, normalize nsig

        :returns: the list of patches, one for each wkspace and the list of list of signals (one list for each json file).
        """
        # Getting the path to the json files
        datasets = [ds.getID() for ds in self.origdatasets]
        jsonFiles = [js for js in self.globalInfo.jsonFiles] #List of json files names (list of str)
        jsons = self.globalInfo.jsons.copy() #List of json files (list of dict)
        if not isinstance(jsons, list):
            logger.error("The `jsons` parameter must be of type list.")
            return None
        # Filtering the json files by looking at the available datasets
        listOfSRInJson = []
        for jsName in self.globalInfo.jsonFiles:
            if all([ds not in self.globalInfo.jsonFiles[jsName] for ds in datasets]):
                # No datasets found for this json combination
                jsIndex = jsonFiles.index(jsName)
                jsonFiles.pop(jsIndex)
                jsons.pop(jsIndex)
                continue
            if not all([ds in datasets for ds in self.globalInfo.jsonFiles[jsName]]):
                # Some SRs are missing for this json combination
                logger.error( "Wrong json definition in globalInfo.jsonFiles for json: %s" % jsName)
            listOfSRInJson += self.globalInfo.jsonFiles[jsName]
        logger.debug("list of datasets: {}".format(datasets))
        logger.debug("jsonFiles after filtering: {}".format(jsonFiles))
        # Constructing the list of signals with subsignals matching each json
        listOfSignals = list()
        for jsName in jsonFiles:
            subSig = list()
            for srName in self.globalInfo.jsonFiles[jsName]:
                try:
                    index = datasets.index(srName)
                except ValueError:
                    line = (
                        f"{srName} signal region provided in globalInfo is not in the list of datasets, {jsName}:{','.join(datasets)}"
                    )
                    raise ValueError(line)
                sig = nsig[index]
                subSig.append(sig)
            listOfSignals.append(subSig)

        if hasattr(self.globalInfo, "includeCRs"):
            includeCRs = self.globalInfo.includeCRs
        else:
            includeCRs = False

        wsInfo = getWSInfo(jsons)
        self.patches = patchMaker(jsons, wsInfo, listOfSignals, includeCRs)
        return self.patches, listOfSignals

    # !TP
    def _getBestStatModel(self, nsig, allow_negative_signal=False):
            """
            find the index of the best expected combination.

            :param dataset: the CombinedDataSet object.
            :param nsig: list of signal yields.
            :param allow_negative_signal: if True, the expected upper limit on mu, used to find the best statistical model, can be negative.

            :return: the minimal apriori expected poi upper limit, and the spey StatisticalModel object that produced that result.
            """
            # Get the list of the names of the signal regions used in the json files
            listOfSRInJson=[]
            for SRnames in self.globalInfo.jsonFiles.values():
                listOfSRInJson += SRnames

            patches, listOfSignals = _getPatches(self, nsig)
            mu_ul_exp_min = np.inf
            # Find best combination of signal regions
            for index, (patch, json) in enumerate(zip(patches,self.globalInfo.jsons)):
                # If the expected signal is 0 for each SR in the combined set of SRs, skip
                if all([sig==0. for sig in listOfSignals[index]]):
                    continue
                # The x-section is at the level of the TheoryPrediction
                # if there are multiple sets of SRs, set a xsec_UL for the whole analysis, i.e. that uses all the SRs,
                # so that the resulting R value Is for the whole analysis
                xsec = sum(nsig)/self.getLumi()
                # It is possible to do differently and to set a xsec_UL on each set of SRs but that is not how it done in SModelS so far
                # xsec = sum(listOfSignals[index])/self.getLumi()
                statModel = get_multi_region_statistical_model(analysis=self.globalInfo.id,
                                                                signal=patch,
                                                                observed=json,
                                                                xsection=xsec
                                                                )
                # If all the SRs are used in the json files and there is only one json files, there is only one statModel.
                # No need to compute mu_ul_exp if not needed.
                if all([ds.dataInfo.dataId in listOfSRInJson for ds in self._datasets]) and len(self.globalInfo.jsons) == 1 and return_mu_ul_exp_min == False:
                    return statModel
                config = statModel.backend.model.config()
                bounds = [(suggested[0]-200,suggested[1]+200) for suggested in config.suggested_bounds]
                if allow_negative_signal:
                    bounds[config.poi_index] = (config.minimum_poi, 100)
                else:
                    bounds[config.poi_index] = (0, 100)
                mu_ul_exp = statModel.poi_upper_limit(expected=ExpectationType.apriori,allow_negative_signal=allow_negative_signal,par_bounds=bounds)
                if mu_ul_exp == None:
                    continue
                elif mu_ul_exp < mu_ul_exp_min:
                    mu_ul_exp_min = mu_ul_exp
                    bestStatModel = statModel

            # Check if a non-combined (uncorrelated) signal region is more contraining than the best combination obtained above
            # Check if a signal region is not in the list of SR names used in the json files
            for sig,ds in zip(nsig,self._datasets):
                ds = ds.dataInfo
                if ds.dataId not in listOfSRInJson:
                    xsec = sig/self.getLumi()
                    # Don't bother to compute eUL again (one could do it again if needed)
                    mu_ul_exp = ds.expectedUpperLimit/xsec
                    if mu_ul_exp < muull_exp_min:
                        logger.info("Best constraining model is a single uncorrelated model.")
                        mu_ul_exp_min = mu_ul_exp
                        bestStatModel = get_uncorrelated_region_statistical_model(observations=float(ds.observedN),
                                                                                backgrounds=float(ds.expectedBG),
                                                                                background_uncertainty=float(ds.bgError),
                                                                                signal_yields=float(sig),
                                                                                xsection=xsec,
                                                                                analysis=self.globalInfo.id,
                                                                                backend="simplified_likelihoods" # simplified likelhood backend by default
                                                                                )
            if mu_ul_exp_min == np.inf:
                logger.error(f'No minimal upper limit on POI found for {self.globalInfo.id}')
                return None
            if return_mu_ul_exp_min:
                return bestStatModel, mu_ul_exp_min
            else:
                return bestStatModel

    # !TP
    def getStatModel(self,
        nsig: Union[np.ndarray, List[Dict[Text, List]], List[float]],
        delta_sys: float = 0.0,
        allow_negative_signal = False
    ):
        """
        Create a statistical model from multibin data.

        :param nsig: number of signal events. For simplified likelihood backend this input can
                       contain `np.array` or `List[float]` which contains signal yields per region.
                       For `pyhf` backend this input expected to be a JSON-patch i.e. `List[Dict]`,
                       see `pyhf` documentation for details on JSON-patch format.
        :param delta_sys: systematic uncertainty on signal. Only used for simplified likelihood backend.
        :param allow_negative_signal: if True, the expected upper limit on mu, used to find the best statistical model, can be negative.

        :return: spey StatisticalModel object

        :raises NotImplementedError: if input patter does not match to any backend specific input option
        """

        if hasattr(self, "statModel"):
            return self.statModel
        else:
            if self.type == "simplified":
                nobs = [x.dataInfo.observedN for x in self.origdatasets]
                cov = self.globalInfo.covariance
                bg = [x.dataInfo.expectedBG for x in self.origdatasets]
                third_moment = self.globalInfo.third_moment if hasattr(self.globalInfo, "third_moment") else None
                xsec = sum(nsig)/self.getLumi()

                self.statModel = get_multi_region_statistical_model(analysis = self.globalInfo.id,
                                                                    signal = nsig,
                                                                    observed = nobs,
                                                                    covariance = cov,
                                                                    nb = bg,
                                                                    third_moment = third_moment,
                                                                    delta_sys = delta_sys,
                                                                    xsection = xsec
                                                                    )
            elif self.type == "pyhf":
                self.statModel = self._getBestStatModel(nsig, allow_negative_signal=allow_negative_signal)
            else:
                logger.error(f'Dataset of type "{self.type}" for analysis {self.globalInfo.id} is not of type "simplified" or "pyhf".')
