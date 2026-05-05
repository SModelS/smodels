#!/usr/bin/env python3

"""
.. module:: nnInterface
   :synopsis: Code that delegates the computation of limits and
   likelihoods to machine learned models

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from typing import Union, Text, Tuple, Callable, Dict, Optional
import copy, os
import numpy as np
import sys
import onnxruntime
from smodels.base.smodelsLogging import logger
from smodels.base.physicsUnits import UnitXSec
from smodels.statistics.basicStats import determineBrentBracket, CLsfromNLL, \
         exponentiateNLL, observed, apriori, aposteriori, NllEvalType, \
         CLsWithErrorsfromNLL
from smodels.statistics.exceptions import SModelSStatisticsError as SModelSError
from scipy import optimize, differentiate
from smodels_utils.helper.terminalcolors import *
from smodels.statistics.nnAdapter import NNAdapter
from smodels.tools.caching import roundCache, lru_cache
from smodels.matching.theoryPrediction import mu_digits

nninfo = {
    "hasgreeted": False,
    "repeat": 0
}

nnSettings = {
    "errs_on_min": True
}

def writeOutYields ( theoryPred,
        filename : Union[os.PathLike,None] = "yields.json" ):
    """ a function for debugging only: writes the actual NN input
    into a file called filename

    :param filename: output file name, if None, then it is
    yields_<massparams>,json

    """

    from smodels.base.physicsUnits import GeV
    masses = []
    for node in theoryPred.smsList[0].nodes:
        if node.particle.isSM:
            continue
        masses.append ( float(node.particle.mass.asNumber(GeV)) )
    if filename == None:
        filename = f"yields_{'_'.join(map(str(masses)))}.json"
    nsig = theoryPred.statsComputer.nsig
    computer = theoryPred.statsComputer.upperLimitComputer
    models = computer.adaptors.keys()
    modelToUse = computer.mostSensitiveModel
    gI = theoryPred.dataset.globalInfo
    Dict = { "anaId": gI.id, "masses": masses, "nsignals": nsig,
             "model": modelToUse,
             "txnames":list( set(map(str,theoryPred.txnames))) }
    dicts = []
    if False: # modelToUse == None:
        for m in models:
            yields = computer.totalYieldsFromSignals( m, 1. )
            scaled_yields = computer.scaleYields ( yields, m )
            nn_input = scaled_yields.tolist()
            Dict["model"]=m
            Dict["nn_input"]=nn_input
            dicts.append ( Dict )

    with open ( filename, "wt" ) as f:
        import json
        d = json.dumps ( dicts, indent=4 )
        f.write ( d )
        f.close()

def clsRootFunc( mu : float, return_type: Text,
             modelToUse : Union[None,str], obj : Callable,
             evaluationType : NllEvalType,
             nll_min : float, nll_minA : float, mu_hat : float,
             nSigma : int, pmSigma : int,
             s_nll_min : Optional[float], s_nll_minA : Optional[float] ) \
        -> Union[float,Tuple]:
    """ the cls root function for the ml models.
    i had to put it as a separate function because i want to be
    able to monkey patch it """
    # at - infinity this should be .95,
    # at + infinity it should -.05
    # Make sure to always compute the correct llhd value (from
    # theoryPrediction)
    # and not used the cached value (which is constant for mu~=1 an mu~=0)
    nllA = obj.nll(mu, modelToUse = modelToUse, asimov = 1,
           pmSigma = 0 )
    s_nllA, s_nll = None, None
    if pmSigma != 0:
        # s_nllA = 0. # abs ( obj.nll( mu, modelToUse = modelToUse, asimov = 1,
                   #       pmSigma = 1 ) - nllA )
        s_nllA = abs ( obj.nll( mu, modelToUse = modelToUse, asimov = 1,
                       pmSigma = 1 ) - nllA )
    if False: # asimov == 1 and evaluationType == observed:
        nll = nllA
        if pmSigma != 0:
            s_nll = s_nllA
    else:
        nll = obj.nll (mu, evaluationType=evaluationType,
            modelToUse = modelToUse, asimov = None,
            pmSigma = 0 )
        if pmSigma != 0:
            s_nll = abs ( obj.nll ( mu, evaluationType=evaluationType,
                              modelToUse = modelToUse, asimov = None,
                              pmSigma = 1 ) - nll )
    if pmSigma == 0:
        ret = None
        if nll is not None and nllA is not None:
            ret = CLsfromNLL( nllA, nll_minA, nll, nll_min, (mu_hat > mu), \
                              return_type=return_type, nSigma = nSigma )
        return ret
    if nll is None or nllA is None:
        ret = None, None
    else:
        ret = CLsWithErrorsfromNLL(nllA, nll_minA, nll, nll_min, \
                   s_nllA, s_nll_minA, s_nll, s_nll_min, (mu_hat > mu), \
                   return_type=return_type, nSigma = nSigma )
    if False:
        print ( )
        print ( f"@@NNI so for mu={mu:.3g} rt={return_type} pmSigma {pmSigma} nll={nll:.2f}+-{s_nll:.2f} nll_minA={nll_minA:.2f}+-{s_nll_minA:.3g} nllA={nllA:.2f}+-{s_nllA:.3g} nll_min={nll_min:.2f}+-{s_nll_min:.3g} big_muhat {mu_hat>mu}" )
        print ( f"@@NNI we got {ret[0]:.2f}+-{ret[1]:.4f}, we will return {ret[0]+pmSigma*ret[1]:.3g}" )
        ret2 = CLsfromNLL( nllA, nll_minA, nll, nll_min, (mu_hat > mu), \
                          return_type=return_type, nSigma = nSigma )
        print ( f"@@NNI old was {ret2}" )
    return ret[0]+pmSigma*ret[1]

class NNData:
    """
    Holds data for use in the machine learned models
    :ivar nsignals: signal predictions list divided into sublists,
    one for each json file
    """

    def __init__(self, nsignals, dataObject ):
        self.nsignals = nsignals  # fb
        self.getTotalYield()
        self.dataObject = dataObject
        self.globalInfo = dataObject.globalInfo
        datasets = [ds.getID() for ds in self.dataObject.origdatasets]
        self.origDataSetOrder = datasets

    def getTotalYield ( self ):
        """ the total yield in all signal regions """
        S = 0
        for signal in self.nsignals.values():
            if isinstance(signal, list):
                for sig in signal:
                    S += sig
            else:
                S += signal
        self.totalYield = S

class NNUpperLimitComputer:
    """
    Class that computes the upper limit using the jsons files and signal
    informations in the 'data' instance of 'NNData'
    """

    def __init__( self, data, cl=0.95, lumi=None,
                  pyhfComputer = None ):
        """

        :param data: instance of 'NNData' holding the signals information
        :param cl: confdence level at which the upper limit is desired
        to be computed
        :param lumi: optional luminosity
        :param pyhfComputer: optional pyhfComputer, if we need to fall back
        for some super regions
        :ivar data: created from data
        :ivar nsignals: signal predictions list divided into sublists,
        one for each json file
        :ivar zeroSignalsFlag: list boolean flags in case all signals are zero
        for a specific json
        :ivar cl: created from cl
        :ivar alreadyBeenThere: boolean flag that identifies when nsignals
        accidentally passes twice at two identical values
        """
        self.pyhfComputer = pyhfComputer
        self.data = data
        # first thing we do, we determine whats the most sensitive model
        self.adaptors = {}
        session_options = {}
        import onnxruntime
        version = onnxruntime.__version__.split(".")
        if int(version[0])<=1 and int(version[1])<=20:
            session_options={ "inter_op_num_threads": 1,
                              "intra_op_num_threads": 1 }
        for onnxfilename,onnxb in self.data.globalInfo.onnxes.items():
            self.adaptors[onnxfilename]=NNAdapter ( onnxb,
                    onnxfilename, session_options = session_options )
        # del self.data.globalInfo.onnxes # we wont need that, thank you
        self.lumi = lumi
        self.nsignals = copy.deepcopy ( self.data.nsignals )
        self.determineMostSensitiveModel()
        logger.debug( f"Signals : {self.nsignals}" )
        # self.zeroSignalsFlag = self.data.zeroSignalsFlag
        self.cl = cl

        self.alreadyBeenThere = (
            False  # boolean to detect wether self.signals has returned to an older value
        )
        self.welcome()

    @lru_cache
    def determineMostSensitiveModel ( self ):
        """ determines the most sensitive model, stores all the ULs
        that were needed to compute that.
        """
        onnxfiles = list(self.adaptors.keys())
        if len(onnxfiles)==1:
            self.mostSensitiveModel = onnxfiles[0]
            ulmu = self.getUpperLimitOnMu ( evaluationType=apriori,
                    modelToUse = self.mostSensitiveModel )
            return
        mumin,mostSensitiveModel=float("inf"),None
        for model in onnxfiles:
            ulmu = float("inf")
            try:
                ulmu = self.getUpperLimitOnMu ( evaluationType=apriori,
                        modelToUse = model )
            except SModelSError as e:
                continue
            if ulmu < mumin:
                mumin = ulmu
                mostSensitiveModel = model
        indices = []
        for idx,model in enumerate(self.data.globalInfo.jsonFiles.keys()):
            if model in self.data.globalInfo.mlModels.keys():
                indices.append ( idx )

        ## the most sensitive model and its upper limit we store separately
        if self.pyhfComputer is not None and indices is not []:
            i_idx, ul, jsonf_ = self.pyhfComputer.upperLimitComputer.getBestCombinationIndex( tuple ( indices ) )
            if ul < mumin:
                mumin = ul
                mostSensitiveModel = i_idx # list(self.data.globalInfo.jsonFiles.keys())[i_idx]
        self.mostSensitiveModel = mostSensitiveModel
        self.mumin = mumin # the smallest expected UL


    def isControlRegion ( self, srname : str, modelToUse : str ) -> bool:
        """ check if srname is control region
        :returns: true if srname is control region
        """
        def isCR ( dct : Dict )-> bool:
            """ tiny helper """
            if not "type" in dct:
                return False
            return dct["type"]=="CR"

        for dct in self.data.globalInfo.mlModels[modelToUse]:
            name = dct["onnx"]
            if name == srname:
                return isCR ( dct )
            pname = name+"-0"
            if pname == srname:
                return isCR ( dct )
        return False

    def totalYieldsFromSignals ( self, modelToUse : str, poi_test : float ) -> list :
        """ given the signal yields self.nsignals, return the total
        yields, signal + background

        :param poi_test: signal strength multiplier to test
        :returns: list of total yields
        """

        yields = []
        for srname,smyield in self.adaptors[modelToUse].onnxMeta["smYields"].items():
            p1 = srname.rfind("-")
            realname = srname[:p1]
            if not realname in self.nsignals:
                realname = f"{realname}[{srname[p1+1:]}]"
                if not realname in self.nsignals:
                    continue
                assert realname in self.nsignals, \
                  f"nnInterface: cannot find sr name {realname} in '{' '.join(self.nsignals.keys())}'"
            # smodelsname = self.data.globalInfo
            signal = float ( self.nsignals[realname]*poi_test )
            if self.isControlRegion ( srname, modelToUse ):
                if hasattr ( self.data.globalInfo, "includeCRs" ) and self.data.globalInfo.includeCRs == False:
                    continue
                obsyield = self.adaptors[modelToUse].onnxMeta["obsYields"][srname]
                ## seems like a CR! replaced bkgexpected with observed (postfit)
                smyield = obsyield
            tot = smyield + signal
            yields.append ( tot )
        return yields

    @roundCache(argname='mu',argpos=1,digits=mu_digits)
    def _actual_nll(self, poi_test : float,
        modelToUse : Union[None,str] = None,
        outputType : str = "extended" ):
        """ the method that really wraps around the llhd computation.
        :param modelToUse: if given, compute the nll for that model.
        If None compute for most sensitive analysis.
        :param outputType: if 'extended' return dictionary with all
        values, if 'observed' return nll_obs_1, if 'expected' return
        nll_exp_1, if 'asimov' return nllA_obs_1, if 'asimov_exp'
        return nllA_exp_1

        :returns: dictionary with nlls, obs and exp, mu=0 and 1
        """
        try:
            poi_test = poi_test[0]
        except (TypeError,IndexError) as e:
            pass

        if modelToUse == None:
            modelToUse = self.mostSensitiveModel
        if modelToUse == None:
            return None
        if modelToUse in self.adaptors:
            # from signal yields compute total yields
            yields = self.totalYieldsFromSignals( modelToUse, poi_test )
            ret = self.adaptors[modelToUse].predict(yields)
        else:
            nll = self.pyhfComputer.upperLimitComputer.nll ( poi_test,
                    modelToUse, evaluationType = observed )
            nllA = self.pyhfComputer.upperLimitComputer.nll ( poi_test,
                    modelToUse, evaluationType = observed, asimov = 1. )
            nllE = self.pyhfComputer.upperLimitComputer.nll ( poi_test,
                    modelToUse, evaluationType = apriori )
            nllEA = self.pyhfComputer.upperLimitComputer.nll ( poi_test,
                    modelToUse, evaluationType = apriori, asimov = 1. )
            nll0 = self.pyhfComputer.upperLimitComputer.nll ( 0.,
                    modelToUse, evaluationType = observed )
            nllA0 = self.pyhfComputer.upperLimitComputer.nll ( 0.,
                    modelToUse, evaluationType = observed, asimov = 1. )
            nllE0 = self.pyhfComputer.upperLimitComputer.nll ( 0.,
                    modelToUse, evaluationType = apriori )
            nllEA0 = self.pyhfComputer.upperLimitComputer.nll ( 0.,
                    modelToUse, evaluationType = apriori, asimov = 1. )
            ret = { "nll_obs_1": nll, "nll_exp_1": nllE,
                    "nllA_obs_1": nllA, "nllA_exp_1": nllEA,
                    "nll_obs_0": nll0, "nll_exp_0": nllE0,
                    "nllA_obs_0": nllA0, "nllA_exp_0": nllEA0 }

        # we return what has been asked
        if outputType == "observed":
            return ret["nll_obs_1"]
        if outputType == "expected":
            return ret["nll_exp_1"]
        if outputType == "asimov":
            return ret["nllA_obs_1"]
        if outputType == "asimov_exp":
            return ret["nllA_exp_1"]
        if outputType != "extended":
            logger.error ( f"outputType {outputType} unknown. should be one of 'observed', 'expected', 'extended'." )
            sys.exit(-1)
        return ret

    def welcome(self):
        """
        greet the world
        """

        if nninfo["hasgreeted"]:
            return
        ver = onnxruntime.__version__
        logger.info( f"NN interface, we are using onnxruntime v{ver}" )
        nninfo["hasgreeted"] = True

    @roundCache(argname='mu',argpos=1,digits=mu_digits)
    def CLs( self, mu : float = 1.0, evaluationType : NllEvalType = observed,
              modelToUse : Union[None,str] = None, asimov : Optional[int] = None,
              allowNegativeSignals : bool = True,
              nSigma : int = 0, pmSigma : int = 0 ) -> Optional[float]:
        mu_hat, sigma_mu, clsRoot, nll_min, nll_minA, \
            s_nll_min, s_nll_minA = self.getCLsRootFunc(
                evaluationType=evaluationType,
                allowNegativeSignals=allowNegativeSignals,
                asimov = asimov,
                modelToUse = modelToUse,
                nSigma = nSigma, pmSigma = pmSigma )
        # we set these errors to zero, they should be strongly correlated
        # with s_nll and s_nllA
        # s_nll_min, s_nll_minA = 0., 0. ## FIXME
        if mu_hat is None:
            return float("inf")
        clsRootArgs = { "mu": mu, "return_type": "CLs",
            "modelToUse": modelToUse,
            "obj": self, "evaluationType" : evaluationType,
            "nll_min": nll_min, "nll_minA": nll_minA, "mu_hat": mu_hat,
            "asimov": asimov,
            "nSigma": nSigma, "pmSigma": pmSigma,
            "s_nll_min": s_nll_min, "s_nll_minA": s_nll_minA }
        return float ( clsRoot ( **clsRootArgs ) )

    @roundCache(argname='mu',argpos=1,digits=mu_digits)
    def nll( self, mu : float = 1.0, evaluationType : NllEvalType = observed,
              modelToUse : Union[None,str] = None, asimov : Optional[int] = None,
              pmSigma : int = 0 ) -> Optional[float]:
        """
        Returns the value of the likelihood. \
        Inspired by the 'pyhf.infer.mle' module but for non-log likelihood

        :param return_nll: if true, return nll, not llhd
        :param evaluationType: one of: observed, apriori, aposteriori
        :param modelToUse: if given, compute likelihood for that model.
        :param asimov: if true, compute for asimov data
        :param pmSigma: usually 0, +1 or -1. get the likelihood
        plus that number of sigmas
        If None compute for most sensitive analysis.
        """
        assert asimov in [ 1, 0, False, None ], \
            f"asimov should be one of: 1, 0, False, None"
        ret = self._actual_nll(mu,modelToUse=modelToUse)
        if ret == None:
            return None
        delta = 0
        poi_test = 1
        if abs(mu)<1e-10:
            poi_test = 0
        # evaluationType == observed and \
        if pmSigma != 0 and not "sigma_obs" in ret:
            ## probably we fell back to pyhf likelihoods,
            # so we return Nones
            return None
        if evaluationType == observed:
            if asimov not in [ False, None ]:
                nll = ret[ f'nllA_obs_{int(asimov)}']
                if pmSigma != 0:
                    delta = ret["sigma_obsA"]
            else:
                nll = ret[ f'nll_obs_{poi_test}']
                if pmSigma != 0:
                    delta = ret["sigma_obs"]
        else:
            if asimov not in [ False, None ]:
                nll = ret[ f'nllA_exp_{int(asimov)}']
                if pmSigma != 0:
                    delta = ret["sigma_expA"]
            else:
                nll = ret[ f'nll_exp_{poi_test}']
                if pmSigma != 0 and poi_test != 0:
                    delta = ret["sigma_exp"]
        nll += pmSigma * delta

        logger.debug( f"Calling likelihood")
        return nll

    def likelihood( self, mu=1.0, return_nll=False, evaluationType=observed,
              modelToUse : Union[None,str] = None, asimov : Optional[int] = None,
              pmSigma : int = 0 ) -> Optional[float]:
        """ over the long run we will want to phase out .likelihood
        interfaces entirely """
        nll = self.nll ( mu=mu, evaluationType=evaluationType,
                modelToUse = modelToUse, asimov = asimov, pmSigma = pmSigma )
        if return_nll:
            return nll
        return self.exponentiateNLL ( nll, True )


    @lru_cache
    def nll_min( self, evaluationType=observed,
              allowNegativeSignals=True,
              modelToUse : Union[None,str] = None,
              asimov : Optional[int] = None ):
        """
        Returns the (negative log) max likelihood

        :param evaluationType: one of: observed, apriori, aposteriori
        :param allowNegativeSignals: if False, then negative nsigs are
        replaced with 0.
        :param modelToUse: if given, compute lmax for that model.
        :param asimov: if true, compute for asimov data
        If None compute for most sensitive analysis.
        """
        if modelToUse == None:
            modelToUse = self.mostSensitiveModel
        # FIXME maximize this one function
        if modelToUse == None:
            modelToUse = self.determineMostSensitiveModel()
        if modelToUse is None:
            print ( f"[nnInterface] no most sensitive model found" )
            # return None
            for model in self.adaptors.keys():
            # for model in self.data.globalInfo.onnxMeta.keys():
                ulmu = self.getUpperLimitOnMu ( evaluationType=apriori, modelToUse = model )
                print ( f"[nnInterface] ulmu({model})={ulmu}" )
            return None
        if not modelToUse in self.adaptors.keys():
            if asimov not in [ False, None ]:
                print ( f"[nnInterface] FIXME fix asimov" )
            return self.pyhfComputer.upperLimitComputer.nll_min ( modelToUse,
                    evaluationType, allowNegativeSignals )
            #print ( f"[nnInterface] no {modelToUse} in {', '.join(self.adaptors.keys())}" )
            #return None
        muhat,nllmin = self.adaptors[modelToUse].onnxMeta["nLL_obs_max"]
        if asimov not in [ False, None ]:
            muhat,nllmin = self.adaptors[modelToUse].onnxMeta["nLLA_obs_max"]
            if evaluationType != observed:
                muhat,nllmin = self.adaptors[modelToUse].onnxMeta["nLLA_exp_max"]
        elif evaluationType != observed:
            muhat,nllmin = self.adaptors[modelToUse].onnxMeta["nLL_exp_max"]

        outputType = "observed"
        if evaluationType == apriori:
            outputType = "expected"
        if evaluationType == aposteriori:
            outputType = "asimov"
        if asimov not in [ False, None ]:
            outputType = "asimov"
        options = { "disp": False, "maxiter": 200 }

        ## FIXME compute sigma_mu, compute via nllA

        def myNLL ( x ):
            if type(x) in [ list, np.array, np.ndarray ]:
                ret = []
                for xi in x:
                    ret.append ( myNLL ( xi ) )
                return np.array ( ret )
            ret = self._actual_nll ( x, modelToUse=modelToUse,
                                     outputType=outputType )
            return ret

        method = "Nelder-Mead"
        initx0s = [ 0., .1, -.1, .3, -.3, 1., -1., 3., -3., 10., -10., 100,-100 ]
        for x0 in initx0s:
            bounds=[(-100,100)] if allowNegativeSignals else [(0,100)]
            if bounds[0][0] > x0:
                bounds = [(x0,100)]
            if bounds[0][1] < x0:
                bounds = [(bounds[0][0],x0)]

            o = optimize.minimize ( self._actual_nll, x0=x0,
                    args=(modelToUse,outputType), tol=1e-8, options = options,
                    method = method, bounds=bounds )
            if o.success == True and o.fun>0:
                muhat, nllmin = o.x[0], o.fun
                o = differentiate.hessian ( myNLL, np.array ( [ muhat ] ) )
                hessian = o.ddf[0][0][0]
                sigma_mu = 0.
                if hessian > 0.:
                    sigma_mu = np.sqrt ( 1. / hessian )

                nll0 = self.nll ( mu=0., evaluationType = evaluationType,
                        modelToUse = modelToUse, asimov = asimov )
                if nll0 < nllmin: # if SM is lower, use it
                    nllmin = nll0
                ret = { "nll_min": float ( nllmin ), "muhat": float ( muhat ),
                        "sigma_mu": float ( sigma_mu ) }
                return ret
            if x0 == initx0s:
                method = "L-BFGS-B"
        logger.warning ( f"could not find nll_min!" )
        return None

    def lmax( self, return_nll=False, evaluationType=observed,
              allowNegativeSignals=True,
              modelToUse : Union[None,str] = None,
              asimov : Optional[int] = None ):
        """
        obsolete, use nll_min
        :param return_nll: if true, return nll, not llhd
        """
        nll_min =  self.nll_min( evaluationType = evaluationType,
                           allowNegativeSignals = allowNegativeSignals,
                           modelToUse = modelToUse, asimov = asimov )
        if return_nll:
            return nll_min
        return exponentiateNLL ( nll_min, doIt = True )

    def getUpperLimitOnSigmaTimesEff(self,
            evaluationType : NllEvalType = observed,
            modelToUse : Union[None,str] = None,
            nSigma : int = 0, **kwargs ) -> UnitXSec:
        """
        Compute the upper limit on the fiducial
        cross section sigma times efficiency:

        :param evaluationType: one of: observed, apriori, aposteriori
        :param modelToUse: if given, compute the nll for that model.
        If None compute for most sensitive analysis.
        :param nSigma: the upper limit for central value (0),
        + 1 sigma, - 1 sigma, etc.  For error bands.
        :return: the upper limit on sigma times eff at 'self.cl' level
        (0.95 by default)
        """
        if modelToUse == None:
            modelToUse = self.mostSensitiveModel
        if not modelToUse in self.adaptors:
            ret = self.pyhfComputer.upperLimitComputer.getUpperLimitOnSigmaTimesEff(
                    evaluationType,modelToUse, nSigma )
            return ret
        if self.data.totalYield == 0.:
            return None
        else:
            ul = self.getUpperLimitOnMu( evaluationType=evaluationType,
                                         modelToUse=modelToUse,
                                         nSigma = nSigma, **kwargs )
            if ul == None:
                return ul
            if self.lumi is None:
                logger.error(f"asked for upper limit on fiducial xsec, but no lumi given with the data")
                return ul
            xsec = self.data.totalYield / self.lumi
            return ul * xsec

    @lru_cache
    def getUpperLimitOnMu(self, evaluationType : NllEvalType = observed,
            allowNegativeSignals : bool = False,
            modelToUse : Union[None,str,int] = None,
            nSigma : int = 0, pmSigma : int = 0 ) -> float:
        """
        Compute the upper limit on the signal strength modifier with:
        - by default, the combination of the workspaces in self.workspaces

        :param evaluationType: one of: observed, apriori, aposteriori
        :param modelToUse: if given, compute the nll for that model.
        If None compute for most sensitive analysis.
        :param nSigma: the upper limit for central value (0),
        + 1 sigma, - 1 sigma, etc.  For error bands.
        :return: the upper limit at 'self.cl' level (0.95 by default)
        """
        if modelToUse == None:
            modelToUse = self.mostSensitiveModel
        if modelToUse not in self.adaptors:
            # so its a pyhf one
            ret = self.pyhfComputer.upperLimitComputer.getUpperLimitOnMu(
                    evaluationType,modelToUse, nSigma )
            return ret
        mu_hat, sigma_mu, clsRoot, nll_min, nll_minA, \
            s_nll_min, s_nll_minA = self.getCLsRootFunc(
                evaluationType=evaluationType,
                allowNegativeSignals=allowNegativeSignals,
                modelToUse = modelToUse,
                nSigma = nSigma, pmSigma = pmSigma )
        # we set these errors to zero, they should be strongly correlated
        # with s_nll and s_nllA
        # s_nll_min, s_nll_minA = 0., 0. ## FIXME
        if mu_hat is None:
            return float("inf")
        clsRootArgs = {"return_type": "CLs-alpha", "modelToUse": modelToUse,
            "obj": self, "evaluationType" : evaluationType,
            "nll_min": nll_min, "nll_minA": nll_minA, "mu_hat": mu_hat,
            "nSigma": nSigma, "pmSigma": pmSigma,
            "s_nll_min": s_nll_min, "s_nll_minA": s_nll_minA }
        try:
            a, b = determineBrentBracket(mu_hat, sigma_mu, clsRoot,
                    allowNegative = allowNegativeSignals, args=clsRootArgs,
                        verbose = False )
        except Exception as e:
            logger.debug ( f"exception {e}" )
            return float("inf")
        mu_lim = optimize.brentq(clsRoot, a, b,
                args = tuple(clsRootArgs.values()), rtol=1e-03, xtol=1e-06 )
        return float ( mu_lim )

    def getCLsRootFunc(self, evaluationType: NllEvalType = observed,
            allowNegativeSignals : bool = True,
            modelToUse : Union[None,str] = None,
            asimov : Optional[int] = None,
            nSigma : int = 0, pmSigma : int = 0 ) -> Tuple[float, float, Callable]:
        """
        Obtain the function "CLs-alpha[0.05]" whose root defines the upper limit,
        plus mu_hat and sigma_mu

        :param evaluationType: one of: observed, apriori, aposteriori
        :param modelToUse: if given, compute the nll for that model.
        :param nSigma: the upper limit for central value (0),
        + 1 sigma, - 1 sigma, etc.  For error bands.
        If None compute for most sensitive analysis.
        """
        # a posteriori expected is needed here
        # mu_hat is mu_hat for signal_rel
        fA = self.nll_min( evaluationType = aposteriori,
                           allowNegativeSignals=allowNegativeSignals,
                           modelToUse = modelToUse, asimov = 1 )
        if fA == None:
            return None, None, None, None, None, None, None
        mu_hatA, sigma_muA, nll_minA = fA["muhat"], fA["sigma_mu"], fA["nll_min"]

        fmin = self.nll_min ( evaluationType=evaluationType,
                              allowNegativeSignals=allowNegativeSignals,
                              modelToUse = modelToUse )

        if fmin == None:
            return None, None, None, None, None, None, None

        mu_hat, sigma_mu, nll_min = \
                ( fmin[x] for x in ["muhat", "sigma_mu", "nll_min"] )
        mu_hat = mu_hat if mu_hat is not None else 0.0
        s_nll_minA, s_nll_min = 0., 0.
        if nnSettings["errs_on_min"] and pmSigma != 0:
            # actually we get better coverage if we dont do this!
            # if we want to compute ULs for nlls +- 1 sigma,
            # we compute the usual mu_hat, but add pmSigma times sigma
            # to nll_min(A)
            s_nll_minA = abs ( self.nll ( mu_hatA,
                    evaluationType = aposteriori , modelToUse = modelToUse,
                    pmSigma = 1, asimov = 1 ) - nll_minA )
            s_nll_min = abs ( self.nll ( mu_hat, evaluationType = evaluationType,
                     modelToUse = modelToUse, pmSigma = 1 ) - nll_min )

        #from smodels.base import runtime
        #useTevatron = runtime.experimentalFeature ( "tevatroncls" )
        #if useTevatron:
        #    return mu_hat, sigma_mu, clsRootTevatron
        return mu_hat, sigma_mu, clsRootFunc, nll_min, nll_minA,\
                s_nll_min, s_nll_minA
