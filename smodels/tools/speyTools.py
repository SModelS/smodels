#!/usr/bin/env python3

"""
.. module:: statsTools
   :synopsis: a module that contains the class responsible for
              all statistical computations. Designed to
              eventually become simply a frontend for spey

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "SpeyComputer" ]

from typing import Union, Text, Dict, List
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
from smodels.tools.simplifiedLikelihoods import LikelihoodComputer, UpperLimitComputer, Data
from smodels.tools.pyhfInterface import PyhfData, PyhfUpperLimitComputer
from smodels.tools.truncatedGaussians import TruncatedGaussians
from smodels.tools.analysesCombinations import AnaCombLikelihoodComputer
from smodels.experiment.datasetObj import DataSet,CombinedDataSet
from typing import Union, Text

class SpeyComputer:
    def __init__ ( self, dataObject : Union['DataSet','CombinedDataSet', list], 
                   backendType : str,
                   nsig : Union[None,float,List] = None, 
                   deltas_rel : Union[None,float] = None, ## possibly not implemented by spey
                   allowNegativeSignals : bool = False):
        """
         Initialise.
        :param dataObject: a smodels (combined)dataset or a list of theory predictions (for combination of analyses)
        :param backendType: which backend to be used (SL, pyhf, ML, ... )
        :param nsig: signal yield, either as float or as list
        :param deltas_rel: relative error on signal. currently unused
        :allowNegativeSignals: if True, negative values for the signal (mu) are allowed.
        """

        self.dataObject = dataObject
        self.nsig = nsig
        self.deltas_sys = deltas_rel
        if self.deltas_sys is None:
            self.deltas_sys = 0.
        self.allowNegativeSignals = allowNegativeSignals

    def get_five_values ( self, expected : Union [ bool, Text ],
                      return_nll : bool = False,
                      check_for_maxima : bool = False )-> Dict:
        """ return the Five Values: l(bsm), l(sm), muhat, l(muhat), sigma(mu_hat) 
        :param check_for_maxima: if true, then check lmax against l(sm) and l(bsm)
             correct, if necessary
        """
        ## call spey to compute this

    def likelihood ( self, poi_test : float, expected : Union[bool,Text],
                            return_nll : bool ) -> float:
        """ simple frontend to individual computers """
        ## call spey to compute this

    def CLs ( self, poi_test : float = 1., expected : Union[bool,Text] = False ) -> Union[float,None]:
        """ compute CLs value for a given value of the poi """
        # now call spey to compute this

    def maximize_likelihood ( self, expected : Union[bool,Text],
           return_nll : bool = False ) -> dict:
        """ simple frontend to the individual computers, later spey
        :param return_nll: if True, return negative log likelihood 
        :returns: Dictionary of llhd (llhd at mu_hat), \
                  muhat, sigma_mu (sigma of mu_hat), \
                  optionally also theta_hat
        """
        # now call spey to compute this

    def poi_upper_limit ( self, expected : Union [ bool, Text ],
           limit_on_xsec : bool = False ) -> float:
        """ simple frontend to the upperlimit computers, later
            to spey::poi_upper_limit
        :param limit_on_xsec: if True, then return the limit on the
                              cross section
        """
        # now call spey to compute this

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
            self.id = "SpeyStatsDataSet"
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
    computer = SpeyComputer ( dataset, 1. )
    ul = computer.poi_upper_limit ( expected = False, limit_on_xsec = True )
    print ( "ul", ul )
    ule = computer.poi_upper_limit ( expected = True, limit_on_xsec = True )
    print ( "ule", ule )
