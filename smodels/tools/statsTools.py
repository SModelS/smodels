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
    def __init__ ( self, dataset, nsig : Union[float,list],
                   deltas_rel : Union[None,float] = None,
                   marginalize : bool = False ):
        """ initialise with dataset.
        :param dataset: a smodels (combined)dataset
        :param nsig: signal yield, either as float or as list
        :param deltas_rel: relative error on signal. currently unused
        """
        #if deltas_rel != None:
        #    logger.warning("Relative uncertainty on signal not supported for a single region.")
        if dataset.getType() not in [ "efficiencyMap", "combined" ]:
            logger.error ( f"I do not recognize the dataset type {dataset.getType()}" )

        self.dataset = dataset
        self.nsig = nsig
        self.getComputer ( nsig, deltas_rel, marginalize )

    def getComputer( self, nsig : Union [ float, np.ndarray ],
                     deltas_rel : Union [ float, None ],
                     marginalize : bool = False ):
        """ retrieve the statistical model """
        if self.dataset.getType() == "efficiencyMap":
            return self.getComputerSingleBin ( nsig, deltas_rel, marginalize )
        # dataset.getType() is "combined"
        assert self.dataset.type in ("simplified", "pyhf" ), \
            f"I do not recognize the datatype {self.dataset.type}. it is not one of simplified, pyhf"
        if self.dataset.type == "simplified":
            return self.getComputerMultiBinSL ( nsig, deltas_rel, marginalize )
        # self.getComputerPyhf ( nsig, deltas_rel, marginalize )

    def getComputerSingleBin(self, nsig: Union[float, np.ndarray],
            delta_sys : Union[None,float] = None,
            marginalize : bool = False
#            allow_negative_signal : bool = False
    ):
        """
        Create computer from a single bin
        :param nsig: signal yields.
        """
        if delta_sys == None:
            delta_sys = 0.
        dataset = self.dataset
        from smodels.tools.simplifiedLikelihoods import LikelihoodComputer, UpperLimitComputer, Data
        data = Data( dataset.dataInfo.observedN, dataset.dataInfo.expectedBG,
                     dataset.dataInfo.bgError**2, deltas_rel = delta_sys )
        self.data = data
        self.marginalize = marginalize
        self.likelihoodComputer = LikelihoodComputer ( data )
        self.upperLimitComputer = UpperLimitComputer ( )

    def getComputerMultiBinSL(self, nsig: Union[float, np.ndarray],
            delta_sys : Union[None,float] = None,
            marginalize : bool = False
    ):
        """
        Create computer from a multi bin SL result
        :param nsig: signal yields.
        """
        if delta_sys == None:
            delta_sys = 0.
        dataset = self.dataset
        cov = dataset.globalInfo.covariance
        nobs = [ x.dataInfo.observedN for x in dataset._datasets ]
        bg = [ x.dataInfo.expectedBG for x in dataset._datasets ]
        from smodels.tools.simplifiedLikelihoods import LikelihoodComputer, UpperLimitComputer, Data
        data = Data( nobs, bg, cov, third_moment=None, nsignal = nsig,
                     deltas_rel = delta_sys, lumi=dataset.getLumi() )
        self.data = data
        self.marginalize = marginalize
        self.likelihoodComputer = LikelihoodComputer ( data )
        self.upperLimitComputer = UpperLimitComputer ( ntoys = 10000 )

    def likelihood ( self, poi_test : float, expected : Union[bool,Text],
                            return_nll : bool ) -> float:
        """ simple frontend to individual computers """
        self.likelihoodComputer.transform ( expected )
        return self.likelihoodComputer.likelihood ( poi_test * self.nsig, nll = return_nll, marginalize = self.marginalize )

    def maximize_likelihood ( self, expected : Union[bool,Text],
           allowNegativeSignals : bool = True,
           return_nll : bool = False  ) -> dict:
        """ simple frontend to the individual computers, later spey
        :param return_nll: if True, return negative log likelihood
        :param allowNegativeSignals: allow also negative muhats
        :returns: Dictionary of llhd (llhd at mu_hat),
                  muhat, sigma_mu (sigma of mu_hat),
                  optionally also theta_hat
        """
        self.likelihoodComputer.transform ( expected )
        lmax = self.likelihoodComputer.lmax ( nll = return_nll,
               allowNegativeSignals = allowNegativeSignals,
               marginalize = self.marginalize )
        self.sigma_mu = self.likelihoodComputer.sigma_mu
        self.muhat = self.likelihoodComputer.muhat
        ret = { "llhd": lmax, "muhat": self.likelihoodComputer.muhat,
                "sigma_mu": self.likelihoodComputer.sigma_mu }
        if hasattr ( self.likelihoodComputer, "theta_hat" ):
            ret["theta_hat"] = self.likelihoodComputer.theta_hat
        return ret

    def poi_upper_limit ( self, expected : Union [ bool, Text ],
           limit_on_xsec : bool = False ) -> float:
        """ simple frontend to the upperlimit computers, later
            to spey::poi_upper_limit
        :param limit_on_xsec: if True, then return the limit on the
                              cross section
        """
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
