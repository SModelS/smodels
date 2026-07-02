#!/usr/bin/env python3

"""
.. module:: nnInterface
   :synopsis: Code that delegates the computation of limits and
   likelihoods to machine learned models

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from typing import Union, Text, Tuple, Callable, Optional
import copy
import numpy as np
import sys
import onnxruntime
from smodels.base.smodelsLogging import logger
from smodels.base.physicsUnits import UnitLumi
from smodels.statistics.basicStats import determineBrentBracket, CLsfromNLL, \
         observed, apriori, aposteriori, NllEvalType, \
         CLsWithErrorsfromNLL
from smodels.statistics.exceptions import SModelSStatisticsError as SModelSError
from scipy import optimize, differentiate
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
        filename : Optional[os.PathLike] = None ):
    """ a function for debugging only: writes the actual NN input
    into a file called filename

    :param filename: output file name, if None, then it is
    yields_<massparams>.json
    """

    from smodels.base.physicsUnits import GeV
    masses = []
    for node in theoryPred.smsList[0].nodes:
        if node.particle.isSM:
            continue
        masses.append ( float(node.particle.mass.asNumber(GeV)) )
    if filename == None:
        filename = f"yields_{'_'.join(map(str,map(int,masses)))}.json"
    gI = theoryPred.dataset.globalInfo
    if "-orig" in gI.id:
        return
    print ( f"[nnInterface] writing yields for {gI.id} to {filename}" )
    dicts = []
    Dict = { "anaId": gI.id, "masses": masses,
             "txnames":list( set(map(str,theoryPred.txnames))) }
    ms = theoryPred.statsComputer.getMostSensitiveModel()
    # import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
    Dict["most_sensitive"]=ms.name
    Dict["ul_min"]=ms.getUpperLimitOnMu()
    Dict["nll0"]=theoryPred.nll ( mu=0., writeYields = False )
    Dict["nll0A"]=theoryPred.nll ( mu=0., evaluationType = observed, asimov = 0 )
    Dict["nll"]=theoryPred.nll ( mu=1., writeYields = False )
    Dict["nllA"]=theoryPred.nll ( mu=1., evaluationType = observed, asimov = 0 )
    Dict["nll_mu5"]=theoryPred.nll ( mu=5., writeYields = False )
    Dict["nllA_mu5"]=theoryPred.nll ( mu=5., evaluationType = observed, asimov = 0, writeYields = False )
    dicts.append ( Dict )

    def removeZeros ( nsig : dict ) -> dict:
        newD = {}
        for k,v in nsig.items():
            if v > 0.:
                newD[k]=v
        return newD

    for computer in theoryPred.statsComputer.subComputers:
        if not hasattr ( computer, "totalYieldsFromSignals" ):
            continue
        Dict = {}
        # m = computer.data
        yields_0 = computer.totalYieldsFromSignals( 0. )
        yields_1 = computer.totalYieldsFromSignals( 1. )
        yields_5 = computer.totalYieldsFromSignals( 5. )
        # scaled_yields = computer.scaleYields ( yields, m )
        # nn_input = scaled_yields.tolist()
        Dict["model"]=computer.name
        Dict["nsignals"]=removeZeros ( computer.nsignals )
        Dict["yields_mu0"]= yields_0
        Dict["yields_mu1"]= yields_1
        Dict["yields_mu5"]= yields_5
        # Dict["nn_input"]=nn_input
        dicts.append ( Dict )

    with open ( filename, "wt" ) as f:
        import json
        d = json.dumps ( dicts, indent=4 )
        f.write ( d )
        f.close()


def clsRootFunc( mu : float, return_type: Text, obj : Callable,
        evaluationType : NllEvalType, nll_min : float, nll_minA : float,
        mu_hat : float, nSigma : int, pmSigma : int,
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
    nllA = obj.nll(mu, asimov = 0, pmSigma = 0, evaluationType = observed )
    # nllA = obj.nll(mu, asimov = 0, pmSigma = 0, evaluationType = evaluationType )
    s_nllA, s_nll = None, None
    if pmSigma != 0:
        # s_nllA = 0. # abs ( obj.nll( mu, asimov = 0,
        #       pmSigma = 1 ) - nllA )
        s_nllA = abs ( obj.nll( mu, asimov = 0, pmSigma = 1,
                       evaluationType = observed ) - nllA )
    if evaluationType == aposteriori:
        nll = nllA
        if pmSigma != 0:
            s_nll = s_nllA
    else:
        nll = obj.nll ( mu, evaluationType=evaluationType,
                        asimov = None, pmSigma = 0 )
        if pmSigma != 0:
            s_nll = abs ( obj.nll ( mu, evaluationType=evaluationType,
                                    asimov = None, pmSigma = 1 ) - nll )
    if pmSigma == 0:
        ret = None
        if nll is not None and nllA is not None:
            ret = CLsfromNLL( nllA, nll_minA, nll, nll_min, (mu_hat > mu), \
                              return_type=return_type, nSigma = nSigma )
        if True and evaluationType == aposteriori: # and abs(mu-1)<.01:
            CLs = ret
            if return_type == "CLs-alpha":
                CLs += 0.05
            print ( f"@@NNIa CLs({mu:.3g}) for nllA {nllA} nll {nll} nll_minA {nll_minA} nll_min {nll_min} mu {mu} mu_hat {mu_hat} nSigma {nSigma} eType {evaluationType} return_type {return_type} pmSigma {pmSigma} CLs {CLs}" )
        return ret
    if nll is None or nllA is None:
        ret = None, None
    else:
        ret = CLsWithErrorsfromNLL(nllA, nll_minA, nll, nll_min, \
                   s_nllA, s_nll_minA, s_nll, s_nll_min, (mu_hat > mu), \
                   return_type=return_type, nSigma = nSigma )
    if True and evaluationType == aposteriori:
        print ( f"@@NNIa CLs for nllA {nllA} nll {nll} nll_minA {nll_minA} nll_min {nll_min} mu {mu} mu_hat {mu_hat} nSigma {nSigma} return_type {return_type} pmSigma {pmSigma} {ret}" )
    return ret[0]+pmSigma*ret[1]

class NNData:
    """
    Holds data for use in the machine learned models
    :ivar nsignals: signal predictions list divided into sublists,
    one for each json file
    """

    def __init__(self, nsignals, dataObject ):
        self.nsignals = nsignals  # fb
        self.dataObject = dataObject
        self.globalInfo = dataObject.globalInfo
        datasets = [ds.getID() for ds in self.dataObject.origdatasets]
        self.origDataSetOrder = datasets
        self.getTotalYield()

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

    def getTotalXSec ( self ):
        return self.totalYield / self.globalInfo.lumi

class NNUpperLimitComputer:
    """
    Class that computes the upper limit using the jsons files and signal
    informations in the 'data' instance of 'NNData'
    """

    def __init__( self, data : NNData, cl : float = 0.95,
            lumi : Optional[UnitLumi] = None, onnxfilename : str = "?" ):
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
        self.data = data
        self.name = onnxfilename
        self.allowNegativeSignals = False
        self.dataType = "nn"
        # first thing we do, we determine whats the most sensitive model
        session_options = {}
        import onnxruntime
        version = onnxruntime.__version__.split(".")
        if int(version[0])<=1 and int(version[1])<=20:
            session_options={ "inter_op_num_threads": 1,
                              "intra_op_num_threads": 1 }
        if onnxfilename not in self.data.globalInfo.cachedModels:
            logger.error ( f"could not find {onnxfilename} among cached models")
            sys.exit(-1)
        onnxb = data.globalInfo.cachedModels[onnxfilename]
        #for onnxfilename,onnxb in self.data.globalInfo.cachedModels.items():
        #    if not onnxfilename.endswith ( ".onnx" ):
        #        continue
        self.adaptor = NNAdapter ( onnxb, onnxfilename,
             session_options = session_options )

        # del self.data.globalInfo.onnxes # we wont need that, thank you
        self.lumi = lumi
        self.nsignals = copy.deepcopy ( self.data.nsignals )
        logger.debug( f"Signals : {self.nsignals}" )
        # self.zeroSignalsFlag = self.data.zeroSignalsFlag
        self.cl = cl

        self.alreadyBeenThere = (
            False  # boolean to detect wether self.signals has returned to an older value
        )
        self.welcome()
        self.checkConsistencyMu0()

    def myNLL ( self, x : list ) -> np.array:
        if type(x) in [ list, np.array, np.ndarray ]:
            ret = []
            for xi in x:
                ret.append ( self.myNLL ( xi ) )
            return np.array ( ret )
        d = self._actual_nll ( x )
        ret = d[ self.str_nll ]
        return ret

    def getSigmaMu ( self, mu : float, str_nll : str ) -> float:
        """ get sigma at mu
        :param str_nll: what nll to ask, e.g. nll_obs_1

        :returns: sigma
        """
        self.str_nll = str_nll
        o = differentiate.hessian ( self.myNLL, np.array ( [ mu ]  ) )
        hessian = o.ddf[0][0][0]
        sigma_mu = 0.
        if hessian > 0.:
            sigma_mu = np.sqrt ( 1. / hessian )
        return sigma_mu

    def getTotalXSec ( self ):
        return self.data.getTotalXSec()

    def isControlRegion ( self, srname : str ) -> bool:
        """ check if srname is control region
        :returns: true if srname is control region
        """
        srname0 = srname
        if srname.endswith ( "-0" ):
            srname0 = srname[:-2]
        for region in self.data.globalInfo.srMappings.values():
            if srname == region["onnx"]:
                ret = region["type"]=="CR"
                return ret
            if srname0 == region["onnx"]:
                ret = region["type"]=="CR"
                return ret
        return False

    def checkConsistencyMu0 ( self ):
        """ when getting predictions for bkg_yields (SRs) and obs_yields (CRs),
        nll_*_mu0 == nll_*_mu1. check for this.
        """
        tolerance = 1e-2 # FIXME eventually we need to lower this number
        nlls = self._actual_nll ( poi_test = 0. )
        errors = {}
        for label in [ "nll_exp", "nll_obs", "nllA_exp", "nllA_obs" ]:
            nll0, nll1 = nlls[f"{label}_0"], nlls[f"{label}_1"]
            err = 2. * abs ( nll1-nll0 ) / ( nll1+nll0 )
            errors[label]=err
            if err > tolerance:
                raise SModelSError ( f"error for {self.name} {label} for mu=0 is too large: {err:.2g}>{tolerance:.1g}" )
        if False:
            print ( f"[nnInterface] consistency check for {self.name}: max relative error={max(errors.values()):.3g}" )
            print ( f"[nnInterface] consistency check for {self.name}: errors: {errors}" )
            print ( f"[nnInterface] consistency check for {self.name}: nlls: {nlls}" )

    def totalYieldsFromSignals ( self, poi_test : float ) -> list :
        """ given the signal yields self.nsignals, return the total
        yields, signal + background

        :param poi_test: signal strength multiplier to test
        :returns: list of total yields
        """

        yields = []
        for srname,smyield in self.adaptor.onnxMeta["bkg_yields"].items():
            p1 = srname.rfind("-")
            realname = srname[:p1]
            if realname not in self.nsignals:
                realname = f"{realname}[{srname[p1+1:]}]"
                if realname not in self.nsignals:
                    continue
                assert realname in self.nsignals, \
                  f"nnInterface: cannot find sr name {realname} in '{' '.join(self.nsignals.keys())}'"
            # smodelsname = self.data.globalInfo
            signal = float ( self.nsignals[realname]*poi_test )
            if self.isControlRegion ( srname ):
                if hasattr ( self.data.globalInfo, "includeCRs" ) and self.data.globalInfo.includeCRs == False:
                    continue
                obsyield = self.adaptor.onnxMeta["obs_yields"][srname]
                ## seems like a CR! replaced bkgexpected with observed (postfit)
                smyield = obsyield
            tot = smyield + signal
            yields.append ( tot )
        return yields

    @roundCache(argname='mu',argpos=1,digits=mu_digits)
    def _actual_nll(self, poi_test : float,
        outputType : Optional[str] = None ):
        """ the method that really wraps around the llhd computation.
        :param outputType: if None return dictionary with all
        values, else supply string, e.g. nll_obs_1, for observed at mu=1.

        :returns: dictionary with nlls, obs and exp, mu=0 and 1
        """
        try:
            poi_test = poi_test[0]
        except (TypeError,IndexError):
            pass

        # from signal yields compute total yields
        yields = self.totalYieldsFromSignals( poi_test )
        ret = self.adaptor.predict(yields)
        if False and poi_test == 0.0:
            print ( f"@@X yields {yields} ret {ret}" )
            print ( f"@@name {self.name}" )

        if outputType == None:
            return ret
        # we return what has been asked
        if outputType in ret:
            return ret[ outputType ]
        logger.error ( f"outputType {outputType} unknown. should be one of: None, nll_obs_1, etc" )
        sys.exit(-1)

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
             allowNegativeSignals : bool = False,
             nSigma : int = 0, pmSigma : int = 0 ) -> float:
        mu_hat, sigma_mu, clsRoot, nll_min, nll_minA, \
            s_nll_min, s_nll_minA = self.getCLsRootFunc(
                evaluationType=evaluationType,
                allowNegativeSignals=allowNegativeSignals,
                nSigma = nSigma, pmSigma = pmSigma )
        # we set these errors to zero, they should be strongly correlated
        # with s_nll and s_nllA
        # s_nll_min, s_nll_minA = 0., 0. ## FIXME
        if mu_hat is None:
            return float("inf")
        clsRootArgs = { "mu": mu, "return_type": "CLs",
            "obj": self, "evaluationType" : evaluationType,
            "nll_min": nll_min, "nll_minA": nll_minA, "mu_hat": mu_hat,
            "nSigma": nSigma, "pmSigma": pmSigma,
            "s_nll_min": s_nll_min, "s_nll_minA": s_nll_minA }
        ret = float ( clsRoot ( **clsRootArgs ) )
        if False and evaluationType == aposteriori:
            print ( f"@@NNI1 CLs({mu:.2g})={ret:.3f} nll_min {nll_min} nll_minA {nll_minA}" )
        return ret

    @roundCache(argname='mu',argpos=1,digits=mu_digits)
    def nll( self, mu : float = 1.0, evaluationType : NllEvalType = observed,
             asimov : Optional[int] = None,
             pmSigma : int = 0 ) -> Optional[float]:
        """
        Returns the value of the likelihood. \
        Inspired by the 'pyhf.infer.mle' module but for non-log likelihood

        :param evaluationType: one of: observed, apriori, aposteriori
        :param asimov: if true, compute for asimov data
        :param pmSigma: usually 0, +1 or -1. get the likelihood
        plus that number of sigmas
        If None compute for most sensitive analysis.
        """
        assert asimov in [ 0, None ], \
            f"asimov {asimov} not acccepted. should be one of: 0, None"
        ret = self._actual_nll(mu)
        if ret == None:
            return None
        poi_test = 1
        if abs(mu)<1e-6:
            poi_test = 0
        # evaluationType == observed and \
        if pmSigma != 0 and "sigma_obs" not in ret:
            ## probably we fell back to pyhf likelihoods,
            # so we return Nones
            return None
        # print ( f"@@0 getting nll asimov={asimov} mu={poi_test} etype={evaluationType}" )
        if evaluationType == observed:
            if asimov not in [ None ]:
                c_label = f'nllA_obs_{poi_test}'
                s_label = "sigma_obsA"
            else:
                c_label = f'nll_obs_{poi_test}'
                s_label = "sigma_obs"
        elif evaluationType == aposteriori:
            ## when asked for a posteriori, we return nllA_obs_1
            c_label = f'nllA_obs_{poi_test}'
            s_label = "sigma_obsA"
        else: # evaluationType = apriori
            if asimov not in [ None ]:
                c_label = f'nllA_exp_{poi_test}'
                s_label = "sigma_expA"
            else:
                c_label = f'nll_exp_{poi_test}'
                s_label = "sigma_exp"
        nll = ret[ c_label ]
        if pmSigma != 0:
            nll += pmSigma * ret[ s_label ]

        logger.debug( "Calling likelihood")
        return nll

    @lru_cache
    def nll_min( self, evaluationType : NllEvalType = observed,
              allowNegativeSignals : bool = True,
              asimov : Optional[int] = None ) -> dict:
        """
        Returns the (negative log) max likelihood

        :param evaluationType: one of: observed, apriori, aposteriori
        :param allowNegativeSignals: if False, then negative nsigs are
        replaced with 0.
        :param asimov: if true, compute for asimov data
        If None compute for most sensitive analysis.
        """
        if evaluationType == aposteriori:
            nll_min = self.nll ( 0., evaluationType = observed, asimov = 0 )
            ret = { "muhat": 0., "nll_min":  nll_min, "sigma_mu": 1. }
            return ret
        obs_v_exp = "obs"
        if evaluationType == apriori:
        # if evaluationType != observed:
            obs_v_exp = "exp"
        A = ""
        if asimov not in [ None ]: #  or evaluationType == aposteriori:
            A = "A"
        str_nll = f"nll{A}_{obs_v_exp}_1"
        self.str_nll = str_nll
        # print ( f"@@NNI nll_min for eType {evaluationType} asimov {asimov}: str_nll {str_nll}" )
        ## these values are global nll_min
        # str_nll = f"nLL{A}_{obs_v_exp}"
        # str_nll_min = f"{str_nll}_max"
        #muhat_g,nllmin_g = self.adaptor.onnxMeta[ str_nll_min ]
        #print ( f"@@  --- muhat_g {muhat_g} nllmin_g {nllmin_g}" )
        options = { "disp": False, "maxiter": 200 }

        ## FIXME compute sigma_mu, compute via nllA


        method = "Nelder-Mead"
        initx0s = [ 0., .1, -.1, .3, -.3, 1., -1., 3., -3., 10., -10., 100,-100 ]
        for x0 in initx0s:
            bounds=[(-100,100)] if allowNegativeSignals else [(0,100)]
            if bounds[0][0] > x0:
                bounds = [(x0,100)]
            if bounds[0][1] < x0:
                bounds = [(bounds[0][0],x0)]

            o = optimize.minimize ( self._actual_nll, x0=x0,
                    args=(str_nll,), tol=1e-8, options = options,
                    method = method, bounds=bounds )
            if o.success == True and o.fun>0:
                muhat, nllmin = o.x[0], o.fun
                sigma_mu = self.getSigmaMu ( muhat, str_nll )

                #nll0 = self.nll ( mu=0., evaluationType = evaluationType,
                #                  asimov = asimov )
                #if nll0 < nllmin: # if SM is lower, use it
                #    nllmin = nll0
                #    muhat = 0.
                ret = { "nll_min": float ( nllmin ), "muhat": float ( muhat ),
                        "sigma_mu": float ( sigma_mu ) }
                # print ( f"@@NNr nll_min ret {ret} allowNegativeSignals {allowNegativeSignals} evaluationType {evaluationType} asimov {asimov} name {self.name}" )
                return ret
            if x0 == initx0s:
                method = "L-BFGS-B"
        logger.warning ( "could not find nll_min!" )
        return None


    @lru_cache
    def getUpperLimitOnMu(self, evaluationType : NllEvalType = observed,
            allowNegativeSignals : bool = False,
            nSigma : int = 0, pmSigma : int = 0 ) -> float:
        """
        Compute the upper limit on the signal strength modifier with:
        - by default, the combination of the workspaces in self.workspaces

        :param evaluationType: one of: observed, apriori, aposteriori
        :param nSigma: the upper limit for central value (0),
        + 1 sigma, - 1 sigma, etc.  For error bands.
        :return: the upper limit at 'self.cl' level (0.95 by default)
        """
        mu_hat, sigma_mu, clsRoot, nll_min, nll_minA, \
            s_nll_min, s_nll_minA = self.getCLsRootFunc(
                evaluationType=evaluationType,
                allowNegativeSignals=allowNegativeSignals,
                nSigma = nSigma, pmSigma = pmSigma )
        # we set these errors to zero, they should be strongly correlated
        # with s_nll and s_nllA
        # s_nll_min, s_nll_minA = 0., 0. ## FIXME
        if mu_hat is None:
            return float("inf")
        clsRootArgs = {"return_type": "CLs-alpha", "obj": self,
            "evaluationType" : evaluationType,
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
        if False and evaluationType == aposteriori:
            print ()
            print ( f"@@NNI0 UL: {self.name} allowN {allowNegativeSignals}" )
            print ( f"@@NNI1 UL: mu_hat {mu_hat} nll_min {nll_min} nll_minA {nll_minA}" )
            print ( f"@@NNI3 get UL {mu_lim}" )
        return float ( mu_lim )

    def getCLsRootFunc(self, evaluationType: NllEvalType = observed,
            allowNegativeSignals : bool = False,
            nSigma : int = 0, pmSigma : int = 0 ) -> Tuple[float, float, Callable]:
        """
        Obtain the function "CLs-alpha[0.05]" whose root defines the upper limit,
        plus mu_hat and sigma_mu

        :param evaluationType: one of: observed, apriori, aposteriori
        :param nSigma: the upper limit for central value (0),
        + 1 sigma, - 1 sigma, etc.  For error bands.
        If None compute for most sensitive analysis.
        """
        # a posteriori expected is needed here
        # mu_hat is mu_hat for signal_rel
        """
        fA = self.nll_min( evaluationType = observed,
            allowNegativeSignals=allowNegativeSignals, asimov = 0 )
        if fA == None:
            return None, None, None, None, None, None, None
        mu_hatA, sigma_muA, nll_minA = fA["muhat"], fA["sigma_mu"], fA["nll_min"]
        print ( f"@@02 mu_hatA {mu_hatA} nll_minA {nll_minA}" )
        """
        s_nll_minA, s_nll_min = 0., 0.
        mu_hatA = 0
        nll_minA = self.nll ( mu = 0, evaluationType = observed, asimov = 0 )
        if evaluationType == aposteriori:
            n2 = self.nll_min ( evaluationType = observed, asimov = 0 )
            print ( f"@@min_A nll0A {nll_minA} n2 {n2}" )
        sigma_muA = 1
        #fmin = self.nll_min ( evaluationType=evaluationType,
        #                      allowNegativeSignals=allowNegativeSignals )
        #if fmin == None:
        #    return None, None, None, None, None, None, None
        mu_hat = 0. # fmin["muhat"]

        if evaluationType == aposteriori:
            # this is not true, right?
            nll_min = nll_minA
            #mu hat = 0.
            # mu_hat = fmin["muhat"]
            sigma_mu = sigma_muA
        else:
            fmin = self.nll_min ( evaluationType=evaluationType,
                                  allowNegativeSignals=allowNegativeSignals )

            #if fmin == None:
            #    return None, None, None, None, None, None, None

            mu_hat, sigma_mu, nll_min = \
                    ( fmin[x] for x in ["muhat", "sigma_mu", "nll_min"] )
            mu_hat = mu_hat if mu_hat is not None else 0.0
        if nnSettings["errs_on_min"] and pmSigma != 0:
            # actually we get better coverage if we dont do this!
            # if we want to compute ULs for nlls +- 1 sigma,
            # we compute the usual mu_hat, but add pmSigma times sigma
            # to nll_min(A)
            s_nll_minA = abs ( self.nll ( mu_hatA,
                    evaluationType = observed,
                    pmSigma = 1, asimov = 0 ) - nll_minA )
            if evaluationType == aposteriori:
                s_nll_min = s_nll_minA
            else:
                s_nll_min = abs ( self.nll ( mu_hat,
                     evaluationType = evaluationType,
                     pmSigma = 1 ) - nll_min )

        #from smodels.base import runtime
        #useTevatron = runtime.experimentalFeature ( "tevatroncls" )
        #if useTevatron:
        #    return mu_hat, sigma_mu, clsRootTevatron
        return mu_hat, sigma_mu, clsRootFunc, nll_min, nll_minA,\
                s_nll_min, s_nll_minA
