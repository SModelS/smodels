#!/usr/bin/env python3

"""
.. module:: statsTools
   :synopsis: a module that contains the class responsible for
              all statistical computations. Designed to
              eventually become simply a frontend for spey

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "StatsComputer" ]

from typing import Union, Text, Tuple, Dict, List
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
import numpy as np

class StatsComputer:
    __slots__ = [ "nsig", "dataset", "data", "marginalize", "likelihoodComputer",
                  "upperLimitComputer", "type" ]
    def __init__ ( self, dataset, nsig : Union[float,list],
                   deltas_rel : Union[None,float] = None,
                   marginalize : bool = False, **kwargs ):
        """ initialise with dataset.
        :param dataset: a smodels (combined)dataset
        :param nsig: signal yield, either as float or as list
        :param deltas_rel: relative error on signal. currently unused
        :param marginalize: profile (False) or marginalize (True)
        """
        #if deltas_rel != None:
        #    logger.warning("Relative uncertainty on signal not supported for a single region.")
        if dataset.getType() not in [ "efficiencyMap", "combined" ]:
            logger.error ( f"I do not recognize the dataset type {dataset.getType()}" )

        self.dataset = dataset
        self.nsig = nsig
        self.getComputer ( nsig, deltas_rel, marginalize, **kwargs )

    def getComputer( self, nsig : Union [ float, np.ndarray ],
                     deltas_rel : Union [ float, None ],
                     marginalize : bool = False,
                     **kwargs ):
        """ retrieve the statistical model """
        self.marginalize = marginalize
        if self.dataset.getType() == "efficiencyMap":
            return self.getComputerSingleBin ( nsig, deltas_rel )
        # dataset.getType() is "combined"
        assert self.dataset.type in ("simplified", "pyhf" ), \
            f"I do not recognize the datatype {self.dataset.type}. it is not one of simplified, pyhf"
        if self.dataset.type == "simplified":
            return self.getComputerMultiBinSL ( nsig, deltas_rel )
        self.getComputerPyhf ( nsig, deltas_rel, **kwargs )

    def getComputerSingleBin(self, nsig: Union[float, np.ndarray],
            delta_sys : Union[None,float] = None,
#            allow_negative_signal : bool = False
    ):
        """
        Create computer from a single bin
        :param nsig: signal yields.
        """
        self.type = "1bin"
        if delta_sys == None:
            delta_sys = 0.
        dataset = self.dataset
        from smodels.tools.simplifiedLikelihoods import LikelihoodComputer, UpperLimitComputer, Data
        data = Data( dataset.dataInfo.observedN, dataset.dataInfo.expectedBG,
                     dataset.dataInfo.bgError**2, deltas_rel = delta_sys )
        self.data = data
        self.likelihoodComputer = LikelihoodComputer ( data )
        self.upperLimitComputer = UpperLimitComputer ( )

    def getComputerMultiBinSL(self, nsig: Union[float, np.ndarray],
            delta_sys : Union[None,float] = None
    ):
        """
        Create computer from a multi bin SL result
        :param nsig: signal yields.
        """
        self.type = "SL"
        if delta_sys == None:
            delta_sys = 0.
        dataset = self.dataset
        cov = dataset.globalInfo.covariance
        nobs = [ x.dataInfo.observedN for x in dataset._datasets ]
        bg = [ x.dataInfo.expectedBG for x in dataset._datasets ]
        third_momenta = [ getattr ( x.dataInfo, "thirdMoment", None ) for x in dataset._datasets ]
        c = third_momenta.count ( None )
        if c > 0:
            if c < len(third_momenta):
                logger.warning ( f"third momenta given for some but not all signal regions in {dataset.globalInfo.id}" )
            third_momenta = None
        from smodels.tools.simplifiedLikelihoods import LikelihoodComputer, UpperLimitComputer, Data
        data = Data( nobs, bg, cov, third_moment=third_momenta, nsignal = nsig,
                     deltas_rel = delta_sys, lumi=dataset.getLumi() )
        self.data = data
        self.likelihoodComputer = LikelihoodComputer ( data )
        self.upperLimitComputer = UpperLimitComputer ( ntoys = 10000 )

    def get_five_values ( self, expected : Union [ bool, Text ],
                      return_nll : bool = False, allowNegativeSignals : bool =False )-> Dict:
        """ return the Five Values: l(bsm), l(sm), muhat, l(muhat), sigma(mu_hat) """
        ret = self.maximize_likelihood ( expected = expected, allowNegativeSignals =allowNegativeSignals, return_nll = return_nll  )
        ret["lmax"] = ret.pop("llhd")
        lbsm = self.likelihood ( poi_test = 1., expected=expected, return_nll = return_nll )
        ret["lbsm"] = lbsm
        lsm = self.likelihood ( poi_test = 0., expected=expected, return_nll = return_nll )
        ret["lsm"] = lsm
        return ret

    def getComputerPyhf(self, nsig: Union[float, np.ndarray],
            delta_sys : Union[None,float] = None, **kwargs
    ):
        """
        Create computer for a pyhf result
        :param nsig: signal yields.
        """
        self.type = "pyhf"
        normalize = kwargs["normalize"] if "normalize" in kwargs else False
        globalInfo = self.dataset.globalInfo
        jsonFiles = [js for js in globalInfo.jsonFiles]
        jsons = globalInfo.jsons.copy()
        # datasets = [ds.getID() for ds in dataset._datasets]
        datasets = [ds.getID() for ds in self.dataset.origdatasets]
        total = sum(nsig)
        # if total == 0.0:  # all signals zero? can divide by anything!
        #     total = 1.0
        if normalize and total != 0.:
            nsig = [
                s / total for s in nsig
            ]  # Normalising signals to get an upper limit on the events count
        # Filtering the json files by looking at the available datasets
        for jsName in globalInfo.jsonFiles:
            if all([ds not in globalInfo.jsonFiles[jsName] for ds in datasets]):
                # No datasets found for this json combination
                jsIndex = jsonFiles.index(jsName)
                jsonFiles.pop(jsIndex)
                jsons.pop(jsIndex)
                continue
            if not all([ds in datasets for ds in globalInfo.jsonFiles[jsName]]):
                # Some SRs are missing for this json combination
                logger.error( "Wrong json definition in globalInfo.jsonFiles for json : %s" % jsName)
        logger.debug("list of datasets: {}".format(datasets))
        logger.debug("jsonFiles after filtering: {}".format(jsonFiles))
        # Constructing the list of signals with subsignals matching each json
        nsignals = list()
        for jsName in jsonFiles:
            subSig = list()
            for srName in globalInfo.jsonFiles[jsName]:
                try:
                    index = datasets.index(srName)
                except ValueError:
                    line = (
                        f"{srName} signal region provided in globalInfo is not in the list of datasets, {jsName}:{','.join(datasets)}"
                    )
                    raise ValueError(line)
                sig = nsig[index]
                subSig.append(sig)
            nsignals.append(subSig)
        # Loading the jsonFiles, unless we already have them (because we pickled)
        from smodels.tools.pyhfInterface import PyhfData, PyhfUpperLimitComputer

        data = PyhfData(nsignals, jsons, total, jsonFiles)
        if data.errorFlag:
            return None
        if hasattr(globalInfo, "includeCRs"):
            includeCRs = globalInfo.includeCRs
        else:
            includeCRs = False
        self.upperLimitComputer = PyhfUpperLimitComputer(data, includeCRs=includeCRs,
                                            lumi=self.dataset.getLumi() )
        self.likelihoodComputer = self.upperLimitComputer # for pyhf its the same

    def likelihood ( self, poi_test : float, expected : Union[bool,Text],
                            return_nll : bool, **kwargs ) -> float:
        """ simple frontend to individual computers """
        self.transform ( expected )
        if self.type == "pyhf":
            if not "workspace_index" in kwargs:
                index = self.likelihoodComputer.getBestCombinationIndex()
                kwargs["workspace_index"] = index
            return self.likelihoodComputer.likelihood (
                    poi_test, nll = return_nll,
                    expected = expected, **kwargs )
        return self.likelihoodComputer.likelihood ( poi_test * np.array (self.nsig),
                nll = return_nll, marginalize = self.marginalize, **kwargs )

    def transform ( self, expected ):
        """ SL only. transform the data to expected or observed """
        if self.type == "pyhf":
            return
        self.likelihoodComputer.transform ( expected )

    def maximize_likelihood ( self, expected : Union[bool,Text],
           allowNegativeSignals : bool = True,
           return_nll : bool = False, **kwargs  ) -> dict:
        """ simple frontend to the individual computers, later spey
        :param return_nll: if True, return negative log likelihood
        :param allowNegativeSignals: allow also negative muhats
        :returns: Dictionary of llhd (llhd at mu_hat),
                  muhat, sigma_mu (sigma of mu_hat),
                  optionally also theta_hat
        """
        self.transform ( expected )
        ret = self.likelihoodComputer.lmax ( nll = return_nll,
               allowNegativeSignals = allowNegativeSignals, **kwargs )
        #self.sigma_mu = self.likelihoodComputer.sigma_mu
        #self.muhat = self.likelihoodComputer.muhat
        #ret = { "llhd": lmax, "muhat": self.likelihoodComputer.muhat,
        #        "sigma_mu": self.likelihoodComputer.sigma_mu }
        #if hasattr ( self.likelihoodComputer, "theta_hat" ):
        #    ret["theta_hat"] = self.likelihoodComputer.theta_hat
        return ret

    def poi_upper_limit ( self, expected : Union [ bool, Text ],
           limit_on_xsec : bool = False ) -> float:
        """ simple frontend to the upperlimit computers, later
            to spey::poi_upper_limit
        :param limit_on_xsec: if True, then return the limit on the
                              cross section
        """
        if self.type == "pyhf":
            index = self.likelihoodComputer.getBestCombinationIndex()
            if limit_on_xsec:
                ret = self.upperLimitComputer.getUpperLimitOnSigmaTimesEff(
                       expected = expected, workspace_index = index )
            else:
                ret = self.upperLimitComputer.getUpperLimitOnMu(
                       expected = expected, workspace_index = index )
        else:
            if limit_on_xsec:
                ret = self.upperLimitComputer.getUpperLimitOnSigmaTimesEff( self.data,
                       expected = expected )
            else:
                ret = self.upperLimitComputer.getUpperLimitOnMu( self.data,
                       expected = expected )
        return ret

class SimpleStatsDataSet:
    """ a very simple data class that can replace a smodels.dataset,
    for 1d SL data only. used for testing and in dataPreparation """
    class SimpleInfo:
        def __init__ ( self, observedN : float, expectedBG : float,
                       bgError : float ):
            self.observedN = observedN
            self.expectedBG = expectedBG
            self.bgError = bgError

    class GlobalInfo:
        def __init__ ( self, lumi ):
            self.id = "SimpleStatsDataSet"
            self.lumi = lumi

    def __init__ ( self, observedN : float, expectedBG : float,
                   bgError : float, lumi : fb ):
        """ initialise the dataset with all relevant stats """
        self.dataInfo = self.SimpleInfo ( observedN, expectedBG, bgError )
        self.globalInfo = self.GlobalInfo( lumi )

    def getLumi ( self ):
        return self.globalInfo.lumi

    def getType ( self ):
        return "efficiencyMap"

if __name__ == "__main__":
    # nobs,bg,bgerr,lumi = 3., 4.1, 0.6533758489567854, 35.9/fb
    # nobs,bg,bgerr,lumi = 0, 1., 0.2, 35.9/fb
    # nobs,bg,bgerr,lumi = 0, 0.001, 0.01, 35.9/fb
    nobs,bg,bgerr,lumi = 3905,3658.3,238.767, 35.9/fb
    dataset = SimpleStatsDataSet ( nobs, bg, bgerr, lumi )
    computer = StatsComputer ( dataset, 1. )
    ul = computer.poi_upper_limit ( expected = False, limit_on_xsec = True )
    print ( "ul", ul )
    ule = computer.poi_upper_limit ( expected = True, limit_on_xsec = True )
    print ( "ule", ule )
