#!/usr/bin/env python3

"""
.. module:: nnInterface
   :synopsis: Code that delegates the computation of limits and likelihoods to machine learned models

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from typing import Union, Text, Tuple, Callable, Dict
import copy, os
import numpy as np
import sys
import onnxruntime
from smodels.base.smodelsLogging import logger
from smodels.statistics.basicStats import determineBrentBracket, CLsfromNLL, CLsfromLsb, getNumericalVariance
from smodels.statistics.exceptions import SModelSStatisticsError as SModelSError
from scipy import optimize, differentiate
from smodels_utils.helper.terminalcolors import *
from smodels.statistics.nnAdapter import NNAdapter

nninfo = {
    "hasgreeted": False,
    "repeat": 0
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
    models = computer.data.globalInfo.onnxMeta.keys()
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
    # import sys, IPython; IPython.embed( colors = "neutral" ) # ; sys.exit()


class NNData:
    """
    Holds data for use in the machine learned models
    :ivar nsignals: signal predictions list divided into sublists, one for each json file
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

    def __init__(self, data, cl=0.95, lumi=None ):
        """

        :param data: instance of 'NNData' holding the signals information
        :param cl: confdence level at which the upper limit is desired to be computed
        :ivar data: created from data
        :ivar nsignals: signal predictions list divided into sublists, one for each json file
        :ivar zeroSignalsFlag: list boolean flags in case all signals are zero for a specific json
        :ivar cl: created from cl
        :ivar alreadyBeenThere: boolean flag that identifies when nsignals accidentally passes twice at two identical values
        """

        self.data = data
        # first thing we do, we determine whats the most sensitive model
        self.adaptors = {}
        for jsonfilename,onnxb in self.data.globalInfo.onnxes.items():
            self.adaptors[jsonfilename]=NNAdapter ( onnxb )
        self.lumi = lumi
        self.nsignals = copy.deepcopy ( self.data.nsignals )
        self.determineMostSensitiveModel()
        logger.debug("Signals : {}".format(self.nsignals))
        # self.zeroSignalsFlag = self.data.zeroSignalsFlag
        self.cl = cl

        self.alreadyBeenThere = (
            False  # boolean to detect wether self.signals has returned to an older value
        )
        self.welcome()

    def determineMostSensitiveModel ( self ):
        """ determines the most sensitive model, stores all the ULs
        that were needed to compute that.
        """
        jsonfiles = list(self.data.globalInfo.onnxMeta.keys())
        # print ( f"@@NN66 determineMostSensitiveModel jsonfiles {jsonfiles}" )
        self.cachedULs= { None: {} }
        if len(jsonfiles)==1:
            self.mostSensitiveModel = jsonfiles[0]
            ulmu = self.getUpperLimitOnMu ( expected=True, modelToUse = self.mostSensitiveModel )
            self.cachedULs[self.mostSensitiveModel] = {}
            self.cachedULs[self.mostSensitiveModel][True]=ulmu
            self.cachedULs[None][True]=ulmu
            return
        mumin,mostSensitiveModel=float("inf"),None
        # print ( f"@@NN66 more than one!" )
        for model in jsonfiles:
            ulmu = float("inf")
            try:
                # print ( f"@@NN99 get ulmu for {model}" )
                ulmu = self.getUpperLimitOnMu ( expected=True, modelToUse = model )
                self.cachedULs[model] = { True: ulmu }
                # print ( f"@@NN54 determineMostSensitiveModel for {model} we have ulmu={mu}" )
            except SModelSError as e:
                self.cachedULs[model] = { True: ulmu }
                # print ( f"@@NN54 determineMostSensitiveModel for {model} we have ulmu={mu}" )
                continue
            if ulmu < mumin:
                mumin = ulmu
                mostSensitiveModel = model
        self.cachedULs[None][True]=mumin # the smallest expected UL
        ## the most sensitive model and its upper limit we store separately
        self.mostSensitiveModel = mostSensitiveModel
        self.mumin = mumin # the smallest expected UL
        # print ( f"@@NN35 determineMostSensitiveModel after {model} {self.mostSensitiveModel}" )

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
            # print ( f"@@NN12 {dct}" )
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
        for srname,smyield in self.data.globalInfo.onnxMeta[modelToUse]["smYields"].items():
            p1 = srname.rfind("-")
            realname = srname[:p1]
            if not realname in self.nsignals:
                realname = f"{realname}[{srname[p1+1:]}]"
                assert realname in self.nsignals, \
                  f"nnInterface: cannot find sr name {realname} in {' '.join(self.nsignals.keys())}"
            # smodelsname = self.data.globalInfo
            signal = float ( self.nsignals[realname]*poi_test )
            if self.isControlRegion ( srname, modelToUse ):
                if hasattr ( self.data.globalInfo, "includeCRs" ) and self.data.globalInfo.includeCRs == False:
                    continue
                obsyield = self.data.globalInfo.onnxMeta[modelToUse]["obsYields"][srname]
                ## seems like a CR! replaced bkgexpected with observed (postfit)
                smyield = obsyield
            tot = smyield + signal
            yields.append ( tot )
        return yields

    def negative_log_likelihood(self, poi_test : float,
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

        # from signal yields compute total yields
        yields = self.totalYieldsFromSignals( modelToUse, poi_test )
        """
        print ( f"@@NNX we evaluate at {yields[:2]}" )
        # we scale these yields
        scaled_yields = self.scaleYields ( yields, modelToUse )
        # we send this through the network
        arr = self.predict ( scaled_yields, modelToUse )
        # from the output we compute the NLLs
        ret = self.nllsFromPrediction ( arr, modelToUse, poi_test )
        """
        ret = self.adaptors[modelToUse].predict(yields)

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
        # print ( f"@@NN22 ret {ret} oldret {self.negative_log_likelihood_old(poi_test,modelToUse,outputType)}" )
        return ret

    def welcome(self):
        """
        greet the world
        """

        if nninfo["hasgreeted"]:
            return
        logger.info( f"NN interface, we are using onnxruntime v{onnxruntime.__version__}" )
        nninfo["hasgreeted"] = True

    def likelihood( self, mu=1.0, return_nll=False, expected=False,
              modelToUse : Union[None,str] = None, asimov : bool = False ):
        """
        Returns the value of the likelihood. \
        Inspired by the 'pyhf.infer.mle' module but for non-log likelihood

        :param return_nll: if true, return nll, not llhd
        :param expected: if False, compute expected values, if True, \
            compute a priori expected, if "posteriori" compute posteriori \
            expected
        :param modelToUse: if given, compute likelihood for that model.
        :param asimov: if true, compute for asimov data
        If None compute for most sensitive analysis.
        """
        ret = self.negative_log_likelihood(mu,modelToUse=modelToUse)
        if ret == None:
            return None
        if expected:
            if asimov:
                nll = ret['nllA_exp_1']
            else:
                nll = ret['nll_exp_1']
        else:
            nll = ret['nll_obs_1']
            if asimov:
                nll = ret['nllA_obs_1']

        logger.debug( f"Calling likelihood")
        return self.exponentiateNLL ( nll, not return_nll )

    def exponentiateNLL(self, nll, doIt = True ):
        """if doIt, then compute likelihood from nll,
        else return nll"""
        if nll == None:
            return None
            #if doIt:
            #    return 0.0
            #return 9000.0
        if doIt:
            return np.exp(-nll )
        return nll

    def lmax( self, return_nll=False, expected=False,
              allowNegativeSignals=True,
              modelToUse : Union[None,str] = None,
              asimov : bool = False ):
        """
        Returns the (negative log) max likelihood

        :param return_nll: if true, return nll, not llhd
        :param expected: if False, compute expected values, if True,
        compute a priori expected, if "posteriori" compute posteriori
        expected
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
            for model in self.data.globalInfo.onnxMeta.keys():
                ulmu = self.getUpperLimitOnMu ( expected=True, modelToUse = model )
                print ( f"[nnInterface] ulmu({model})={ulmu}" )
            return None
        if not modelToUse in self.data.globalInfo.onnxMeta:
            print ( f"[nnInterface] no {modelToUse} in {', '.join(self.data.globalInfo.onnxMeta.keys())}" )
            return None
        muhat,nllmin = self.data.globalInfo.onnxMeta[modelToUse]["nLL_obs_max"]
        if asimov:
            muhat,nllmin = self.data.globalInfo.onnxMeta[modelToUse]["nLLA_obs_max"]
            if expected:
                muhat,nllmin = self.data.globalInfo.onnxMeta[modelToUse]["nLLA_exp_max"]
        elif expected:
            muhat,nllmin = self.data.globalInfo.onnxMeta[modelToUse]["nLL_exp_max"]

        outputType = "observed"
        if expected == True:
            outputType = "expected"
        if expected == "posteriori":
            outputType = "asimov"
        #print ( f"@@NNO outputType {outputType} " )
        #bounds=[(-1,10)]
        #def callme ( args ):
        #    print ( f"@@NN7 callme {args}" )
        options = { "disp": False, "maxiter": 200 }

        ## FIXME compute sigma_mu, compute via nllA

        def myNLL ( x ):
            # print ( f"@@NN8 myNLL x={x}" )
            if type(x) in [ list, np.array, np.ndarray ]:
                ret = []
                for xi in x:
                    ret.append ( myNLL ( xi ) )
                return np.array ( ret )
            ret = self.negative_log_likelihood ( x, modelToUse=modelToUse,
                                                 outputType=outputType )
            # print ( f"@@NN90 myNLL ret={x}" )
            return ret

        method = "Nelder-Mead"
        initx0s = [ 0., .1, -.1, .3, -.3, 1., -1., 3., -3., 10., -10., 100,-100 ]
        bounds=[(-100,100)]
        if not allowNegativeSignals:
            bounds=[(0,100)]
        for x0 in initx0s:
            o = optimize.minimize ( self.negative_log_likelihood, x0=x0,
                    args=(modelToUse,outputType), tol=1e-8, options = options,
                    method = method, bounds=bounds )
            #if o.fun < 0:
            #    print ( f"@@NN44 o {o}" )
            if o.success == True and o.fun>0:
                muhat, nllmin = o.x[0], o.fun
                # import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
                o = differentiate.hessian ( myNLL, np.array ( [ muhat ] ) )
                hessian = o.ddf[0][0][0]
                sigma_mu = 0.
                if hessian > 0.:
                    sigma_mu = np.sqrt ( 1. / hessian )

                ret = { "nll_min": nllmin, "muhat": muhat, "sigma_mu": sigma_mu }
                # print ( f"@@NN45 muhat {muhat} nllmin {nllmin} hessian {hessian} allowNegativeSignals {allowNegativeSignals}" )
                # print ( f"@@NN47 negative_log_likelihood {self.negative_log_likelihood(muhat)}" )
                return ret
            if x0 == initx0s:
                method = "L-BFGS-B"
        #print ( f"could not find nll_min" )
        #print ( f"nll(0) {self.negative_log_likelihood(0.)}" )
        #print ( f"nll(.1) {self.negative_log_likelihood(.1)}" )
        #print ( f"nll(-.1) {self.negative_log_likelihood(-.1)}" )
        #sys.exit()
        logger.warning ( f"could not find nll_min!" )
        return None

    def getUpperLimitOnSigmaTimesEff(self, expected=False,
            modelToUse : Union[None,str] = None ):
        """
        Compute the upper limit on the fiducial cross section sigma times efficiency:

        :param expected:  - if set to 'True': uses expected SM backgrounds as signals
                          - else: uses 'self.nsignals'
        :param modelToUse: if given, compute the nll for that model.
        If None compute for most sensitive analysis.
        :return: the upper limit on sigma times eff at 'self.cl' level (0.95 by default)
        """
        if self.data.totalYield == 0.:
            return None
        else:
            ul = self.getUpperLimitOnMu( expected=expected, modelToUse=modelToUse )
            if ul == None:
                return ul
            if self.lumi is None:
                logger.error(f"asked for upper limit on fiducial xsec, but no lumi given with the data")
                return ul
            xsec = self.data.totalYield / self.lumi
            return ul * xsec

    def getUpperLimitOnMu(self, expected=False, allowNegativeSignals = False,
            modelToUse : Union[None,str] = None ) -> float:
        """
        Compute the upper limit on the signal strength modifier with:
            - by default, the combination of the workspaces contained into self.workspaces

        :param expected:  - if set to 'True': uses expected SM backgrounds as signals
                          - else: uses 'self.nsignals'
        :param modelToUse: if given, compute the nll for that model.
        If None compute for most sensitive analysis.
        :return: the upper limit at 'self.cl' level (0.95 by default)
        """
        #nninfo["repeat"]=nninfo["repeat"]+1
        #if nninfo["repeat"]>10:
        #    sys.exit()
        # print ( f"@@NN13 getUpperLimitOnMu modelToUse {modelToUse}" )
        if modelToUse in self.cachedULs:
            if expected in self.cachedULs[modelToUse]:
                return self.cachedULs[modelToUse][expected]
        else:
            self.cachedULs[modelToUse]={}
        if modelToUse == None:
            modelToUse = self.mostSensitiveModel
        mu_hat, sigma_mu, clsRoot = self.getCLsRootFunc(expected=expected,
                              allowNegativeSignals=allowNegativeSignals,
                              modelToUse = modelToUse)
        if mu_hat is None:
            return float("inf")
        # print ( f"@@NN76 clsRoot expected {expected} modelToUse {modelToUse}" )
        clsRootArgs = {"return_type": "CLs-alpha", "modelToUse": modelToUse }
        try:
            a, b = determineBrentBracket(mu_hat, sigma_mu, clsRoot,
                    allowNegative = allowNegativeSignals, args=clsRootArgs,
                        verbose = True )
        except Exception as e:
            return float("inf")
        mu_lim = optimize.brentq(clsRoot, a, b, args = tuple(clsRootArgs.values()), rtol=1e-03, xtol=1e-06)
        if False: # False and expected == "posteriori":
            print ( f"@@NN473 getUpperLimitOnMu r={1./mu_lim:.3f} expected {expected}" )
        if not modelToUse in self.cachedULs:
            self.cachedULs[modelToUse]={}
        self.cachedULs[modelToUse][expected]=mu_lim # store
        return mu_lim

    def getCLsRootFunc(self, expected: bool = False,
            allowNegativeSignals : bool = True,
            modelToUse : Union[None,str] = None ) -> Tuple[float, float, Callable]:
        """
        Obtain the function "CLs-alpha[0.05]" whose root defines the upper limit,
        plus mu_hat and sigma_mu

        :param expected: if True, compute expected likelihood, else observed
        :param modelToUse: if given, compute the nll for that model.
        If None compute for most sensitive analysis.
        """
        # print ( f"@@NN88 getCLsRootFunc modelToUse {modelToUse}" )
        # a posteriori expected is needed here
        # mu_hat is mu_hat for signal_rel
        fmh = self.lmax(expected="posteriori", allowNegativeSignals=allowNegativeSignals,
                             return_nll=True, modelToUse = modelToUse )
        if fmh == None:
            return None, None, None
        mu_hat, sigma_mu, nll0A = fmh["muhat"], fmh["sigma_mu"], fmh["nll_min"]

        nll0 = nll0A

        if True: # expected != "posteriori":
            fmh = self.lmax(expected=expected, allowNegativeSignals=allowNegativeSignals, modelToUse = modelToUse )
            mu_hat, sigma_mu, nll0 = fmh["muhat"], fmh["sigma_mu"], fmh["nll_min"]
            mu_hat = mu_hat if mu_hat is not None else 0.0
        if False: # expected == "posteriori":
            fmh = self.lmax(expected=expected, allowNegativeSignals=allowNegativeSignals, modelToUse = modelToUse )
            mu_hat, sigma_mu, nll0 = fmh["muhat"], fmh["sigma_mu"], fmh["nll_min"]
            mu_hat = mu_hat if mu_hat is not None else 0.0
            #nll02 = self.likelihood(.298, expected=False, return_nll=True,
            #        modelToUse = modelToUse )
            #print ( f"@@nnInterface.D nll0 {nll0} nll(.298) {nll02}" )

        # logger.error ( f"COMB nll0A {nll0A:.3f} mu_hatA {mu_hatA:.3f}" )
        # return 1.
        # print ( f"@@NNX getCLsRootFunc expected {expected} nll0 {nll0:.1f} nll0A {nll0A:.1f}" )

        def clsRootTevatron( mu: float, return_type: Text = "CLs-alpha",
                     modelToUse : Union[None,str] = None ) -> float:
            # at - infinity this should be .95,
            # at + infinity it should -.05
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            nll_b = self.likelihood(0., return_nll=True, expected=expected, modelToUse = modelToUse )
            if nll_b is None:
                return None
            sigma2_b = getNumericalVariance ( self, 0. )
            nll_sb = self.likelihood(mu, return_nll=True, expected=expected, modelToUse = modelToUse )
            if nll_sb is None:
                return None
            sigma2_sb = getNumericalVariance ( self, 1. )
            CLs = CLsfromLsb(nll_sb, nll_b, sigma2_sb, sigma2_b, return_type=return_type)
            return CLs

        def clsRootAsimov( mu: float, return_type: Text = "CLs-alpha",
                     modelToUse : Union[None,str] = None ) -> float:
            # at - infinity this should be .95,
            # at + infinity it should -.05
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            # print ( f"@@NN732 clsRootAsimov modelToUse {modelToUse}" )
            nllA = self.likelihood(mu, return_nll=True, modelToUse = modelToUse, asimov = True )
            nll = nllA
            if expected != "posteriori":
                nll = self.likelihood(mu, return_nll=True, expected=expected, modelToUse = modelToUse, asimov = False )
            ret =  CLsfromNLL(nllA, nll0A, nll, nll0, (mu_hat > mu), return_type=return_type) if (nll is not None and nllA is not None) else None
            if False: #  and expected == "posteriori" and abs(mu-.765)<.1:
                print ( f"@@NN653 {RED}expected {expected} mu {mu:.3f} nllA {nllA:.3f} nll0A {nll0A:.3f} nll {nll:.3f} nll0 {nll0:.3f} muhat {mu_hat:.3f} CLs {ret} model {modelToUse} {RESET}" )
            return ret

        #from smodels.base import runtime
        #useTevatron = runtime.experimentalFeature ( "tevatroncls" )
        #if useTevatron:
        #    return mu_hat, sigma_mu, clsRootTevatron
        return mu_hat, sigma_mu, clsRootAsimov
