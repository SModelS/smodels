#!/usr/bin/env python3

"""
.. module:: testStatistics
   :synopsis: Tests the statistics functionality.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Jory Sonneveld <jory@opmijnfiets.nl>

"""
import sys

sys.path.insert(0, "../")
import unittest

# from smodels.tools import statistics
from smodels.tools.simplifiedLikelihoods import UpperLimitComputer, LikelihoodComputer, Data
from smodels.tools.srCombinations import getCombinedUpperLimitFor, getCombinedPyhfStatistics
from smodels.tools.statistics import TruncatedGaussians
from smodels.theory.theoryPrediction import theoryPredictionsFor, _getCombinedResultFor, _getDataSetPredictions
from smodels.tools.statistics import determineBrentBracket
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from databaseLoader import database
from smodels.experiment.databaseObj import Database
from smodels.theory import decomposer
from math import floor, log10
import numpy as np
import math
from spey import get_uncorrelated_region_statistical_model, get_multi_region_statistical_model, ExpectationType, AvailableBackends
from spey.hypothesis_testing.utils import find_root_limits, compute_confidence_level # to be removed in the end
from spey.hypothesis_testing.test_statistics import compute_teststatistics # to be removed in the end
from scipy import optimize
from pyhf import Workspace
from pyhf.infer import hypotest
from pyhf.infer.mle import fit, fixed_poi_fit # remove fixed_poi_fit in the end
from jsonpatch import apply_patch
from smodels.tools.physicsUnits import fb,pb


class StatisticsTest(unittest.TestCase):
    def lLHDFromLimits(self):
        """to do some statistics on the chi2"""
        allow_negative_signal = False
        for backendNumber in [1,2]:
            if backendNumber == 1:
                print("PYHF BACKEND.")
            elif backendNumber == 2:
                print("SIMPLIFIED LIKELIHOOD BACKEND.")
            nsig = 1.0
            nobs, nbg = 100.0, 100.0
            statModel = get_uncorrelated_region_statistical_model(observations=nobs,
                                                                    backgrounds=nbg,
                                                                    background_uncertainty=np.sqrt(0.001),
                                                                    signal_yields=nsig,
                                                                    xsection=None,
                                                                    analysis="UnitTest",
                                                                    backend=AvailableBackends(backendNumber)
                                                                )
            config = statModel.backend.model.config()
            bounds = config.suggested_bounds
            if allow_negative_signal:
                bounds[config.poi_index] = (config.minimum_poi, 100)
            else:
                bounds[config.poi_index] = (0, 100)
            ulobs = statModel.poi_upper_limit(expected=ExpectationType.observed,par_bounds=bounds)
            ulexp = statModel.poi_upper_limit(expected=ExpectationType.aposteriori,par_bounds=bounds)
            print("ulobs", ulobs)
            print("ulexp", ulexp)
            f = open("llhds.csv", "wt")
            dx = 0.5
            totdir, totlim, totmarg = 0.0, 0.0, 0.0
            for nsig in np.arange(0.1, 100.0, dx):
                print()
                print("nsig=", nsig)
                statModel = get_uncorrelated_region_statistical_model(observations=nobs,
                                                                        backgrounds=nbg,
                                                                        background_uncertainty=np.sqrt(0.001),
                                                                        signal_yields=nsig,
                                                                        xsection=None,
                                                                        analysis="UnitTest",
                                                                        backend=AvailableBackends(backendNumber)
                                                                    )
                llhddir = statModel.likelihood(nsig,return_nll=False)
                chi2dir = statModel.chi2()
                llhdmarg = None
                ### Marginalization does not work but we do not care for the moment.
                # if backendNumber == 2: # pyhf cannot marginalize?
                #     llhdmarg = np.exp(-statModel.likelihood(nsig, marginalize=True))
                #     args = {"marginalize": True}
                #     # Spey doesn't allow for marginalized chi2 yet
                #     chi2marg = statModel.chi2(**args)
                #     print("llhd marg", llhdmarg, chi2marg)
                computer = TruncatedGaussians ( ulobs, ulexp, nsig )
                ret = computer.likelihood ( mu=1.)
                llhdlim, muhat, sigma_mu = ret["llhd"], ret["muhat"], ret["sigma_mu"]
                chi2lim = computer.chi2 ( llhdlim )
                print("llhd from limits (NOT A SPEY FEATURE)", llhdlim, chi2lim)
                totdir += llhddir * dx
                totlim += llhdlim * dx
                # if backendNumber == 2:
                #     totmarg += llhdmarg * dx
                f.write("%s,%s,%s,%s\n" % (nsig, llhddir, llhdlim, llhdmarg))
            print("total direct", totdir)
            print("total limit", totlim)
            # if backendNumber == 1:
            #     print("No marginalization with pyhf.")
            # elif backendNumber == 2:
            #     print("total marg", totmarg)
            f.close()

    def testUpperLimitOnSigmaTimesEff_SL_1(self):
        """
        Test the cross-section upper limit computation with simplified likelhood backend from SModelS and Spey packages.
        """
        database='./testFilesTim/test_db/'
        # Set the path to the database
        database = Database(database)
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2018-41"])[0]

        filename = './testFilesTim/wino_Spectrum_640_95'

        model = Model(BSMList, SMList)
        model.updateParticles(filename)
        smstoplist = decomposer.decompose(model, sigmacut=0.)
        expected = False
        allow_negative_signal = False

        #SL from SModelS
        dataSetResults = []
        # Compute predictions for each data set (for UL analyses there is one single set)
        for dataset in expRes.datasets:
            predList = _getDataSetPredictions(dataset, smstoplist, maxMassDist=0.2)
            if predList:
                dataSetResults.append(predList)
        combinedRes = _getCombinedResultFor(dataSetResults, expRes, marginalize=False)

        srNsigDict = dict( [ [pred.dataset.getID(),(pred.xsection.value * pred.dataset.getLumi()).asNumber()] for pred in combinedRes.datasetPredictions] )
        nsig = [srNsigDict[ds.getID()] if ds.getID() in srNsigDict else 0.0 for ds in combinedRes.dataset.origdatasets]
        computer = UpperLimitComputer(ntoys=10000)
        dataset = combinedRes.dataset
        cov = dataset.globalInfo.covariance
        nobs = [x.dataInfo.observedN for x in dataset._datasets]
        bg = [x.dataInfo.expectedBG for x in dataset._datasets]
        d = Data(
            observed=nobs,
            backgrounds=bg,
            covariance=cov,
            third_moment=None,
            nsignal=nsig,
            deltas_rel=0.2,
            lumi=dataset.getLumi(),
        )
        xsec = sum(d.nsignal) / d.lumi

        mu_hat_SL, sigma_mu_SL, clsRoot_SL = computer.getCLsRootFunc(d, marginalize=False, expected=expected)
        a, b = determineBrentBracket(mu_hat_SL, sigma_mu_SL, clsRoot_SL, allowNegative=allow_negative_signal )
        mu_ul_SL = optimize.brentq(clsRoot_SL, a, b, rtol=1e-03, xtol=1e-06)
        xsec_ul_SL = mu_ul_SL*xsec
        computer = LikelihoodComputer(d)
        llhdMax_SL = computer.likelihood(mu=mu_hat_SL,nll=False)

        # Spey computation
        expectedDict = {False:ExpectationType.observed,
                        True:ExpectationType.apriori,
                        "posteriori":ExpectationType.aposteriori}
        statModel = get_multi_region_statistical_model(analysis=dataset.globalInfo.id,
                                                        signal=nsig,
                                                        observed=nobs,
                                                        covariance=cov,
                                                        nb=bg,
                                                        third_moment=None,
                                                        delta_sys=0.2,
                                                        xsection=xsec
                                                        )
        config = statModel.backend.model.config()
        bounds = config.suggested_bounds
        if allow_negative_signal:
            bounds[config.poi_index] = (config.minimum_poi, 100)
        else:
            bounds[config.poi_index] = (0, 100)
        mu_ul_spey = statModel.poi_upper_limit(expected=expectedDict[expected],allow_negative_signal=allow_negative_signal,par_bounds=bounds)
        xsec_ul_spey = mu_ul_spey*statModel.xsection

        mu_hat_spey, llhdMax_spey = statModel.maximize_likelihood(allow_negative_signal=allow_negative_signal, expected=ExpectationType.observed, return_nll=False)

        self.assertAlmostEqual(xsec_ul_SL._value,xsec_ul_spey._value,3)
        self.assertAlmostEqual(mu_hat_SL,mu_hat_spey,4)
        self.assertAlmostEqual(llhdMax_SL,llhdMax_spey,4)

    def testUpperLimitOnSigmaTimesEff_SL_2(self):
        """
        Test the cross-section upper limit computation with simplified likelhood backend from SModelS and Spey packages.
        """
        database='./testFilesTim/test_db/'
        # Set the path to the database
        database = Database(database)
        expRes = database.getExpResults(analysisIDs=["CMS-SUS-21-002-agg"])[0]

        filename = './testFilesTim/wino_Spectrum_640_95'

        model = Model(BSMList, SMList)
        model.updateParticles(filename)
        smstoplist = decomposer.decompose(model, sigmacut=0.)
        expected = False
        allow_negative_signal = True

        #SL from SModelS
        dataSetResults = []
        # Compute predictions for each data set (for UL analyses there is one single set)
        for dataset in expRes.datasets:
            predList = _getDataSetPredictions(dataset, smstoplist, maxMassDist=0.2)
            if predList:
                dataSetResults.append(predList)
        combinedRes = _getCombinedResultFor(dataSetResults, expRes, marginalize=False)

        srNsigDict = dict( [ [pred.dataset.getID(),(pred.xsection.value * pred.dataset.getLumi()).asNumber()] for pred in combinedRes.datasetPredictions] )
        nsig = [srNsigDict[ds.getID()] if ds.getID() in srNsigDict else 0.0 for ds in combinedRes.dataset.origdatasets]

        computer = UpperLimitComputer(ntoys=10000)
        dataset = combinedRes.dataset
        cov = dataset.globalInfo.covariance
        nobs = [x.dataInfo.observedN for x in dataset._datasets]
        bg = [x.dataInfo.expectedBG for x in dataset._datasets]

        d = Data(
            observed=nobs,
            backgrounds=bg,
            covariance=cov,
            third_moment=None,
            nsignal=nsig,
            deltas_rel=0.2,
            lumi=dataset.getLumi(),
        )
        xsec = sum(d.nsignal) / d.lumi

        mu_hat_SL, sigma_mu_SL, clsRoot_SL = computer.getCLsRootFunc(d, marginalize=False, expected=expected)
        a, b = determineBrentBracket(mu_hat_SL, sigma_mu_SL, clsRoot_SL, allowNegative=allow_negative_signal )
        mu_ul_SL = optimize.brentq(clsRoot_SL, a, b, rtol=1e-03, xtol=1e-06)
        xsec_ul_SL = mu_ul_SL*xsec

        computer = LikelihoodComputer(d)
        llhdMax_SL = computer.likelihood(mu=mu_hat_SL,nll=False)


        # Spey computation
        expectedDict = {False:ExpectationType.observed,
                        True:ExpectationType.apriori,
                        "posteriori":ExpectationType.aposteriori}
        statModel = get_multi_region_statistical_model(analysis=dataset.globalInfo.id,
                                                        signal=nsig,
                                                        observed=nobs,
                                                        covariance=cov,
                                                        nb=bg,
                                                        third_moment=None,
                                                        delta_sys=0.2,
                                                        xsection=xsec
                                                        )
        config = statModel.backend.model.config()
        bounds = config.suggested_bounds
        if allow_negative_signal:
            bounds[config.poi_index] = (config.minimum_poi, 100)
        else:
            bounds[config.poi_index] = (0, 100)

        mu_ul_spey = statModel.poi_upper_limit(expected=ExpectationType.observed,allow_negative_signal=allow_negative_signal,par_bounds=bounds)
        xsec_ul_spey = mu_ul_spey*statModel.xsection

        mu_hat_spey, llhdMax_spey = statModel.maximize_likelihood(allow_negative_signal=allow_negative_signal, expected=ExpectationType.observed, return_nll=False)

        self.assertAlmostEqual(xsec_ul_SL._value,xsec_ul_spey._value,3)
        self.assertAlmostEqual(mu_hat_SL,mu_hat_spey,4)
        self.assertAlmostEqual(llhdMax_SL,llhdMax_spey,4)

        # Check consistency on CLs
        theta_hat0, _ = computer.findThetaHat( 0. )
        sigma_mu = computer.getSigmaMu(mu_hat_SL, theta_hat0)

        nll0 = computer.likelihood( mu_hat_SL, marginalize=False, nll=True)
        import copy
        aModel = copy.deepcopy(d)
        aModel.observed = np.array([x + y for x, y in zip(d.backgrounds, theta_hat0)])
        aModel.name = aModel.name + "A"
        compA = LikelihoodComputer(aModel, 10000)
        mu_hatA = compA.findMuHat()
        nll0A = compA.likelihood( mu=mu_hatA, marginalize=False, nll=True)

        nll = computer.likelihood(mu=1., marginalize=False, nll=True)
        nllA = compA.likelihood(mu=1., marginalize=False, nll=True)

        from smodels.tools.statistics import CLsfromNLL
        CLs_SL = CLsfromNLL(nllA, nll0A, nll, nll0, return_type="1-CLs")

        CLs_spey = statModel.exclusion_confidence_level(allow_negative_signal=allow_negative_signal)[0]
        self.assertAlmostEqual(CLs_SL,CLs_spey,4)

        # Check consitency on qmu and qmuA
        qmu = 0.0 if 2 * (nll - nll0) < 0.0 else 2 * (nll - nll0)
        sqmu_SL = np.sqrt(qmu)
        qA = 2 * (nllA - nll0A)
        if qA < 0.0:
            qA = 0.0
        sqA_SL = np.sqrt(qA)

        def _tmu(mu, max_logpdf, logpdf) -> float:
            return -2 * (logpdf(mu) - max_logpdf)
        def _tmu_tilde(mu, muhat, min_logpdf, logpdf):
            return -2 * (logpdf(mu) - (min_logpdf if muhat >= 0.0 else logpdf(0.0)))
        def qmu_tilde(mu, muhat, min_logpdf, logpdf):
            return 0.0 if muhat > mu else _tmu_tilde(mu, muhat, min_logpdf, logpdf)
        def qmu(mu, muhat, max_logpdf, logpdf) -> float:
            return 0.0 if muhat > mu else _tmu(mu, max_logpdf, logpdf)

        def get_test_statistic(test_stat):
            """
            Retrieve test statistic function

            :raises UnknownTestStatistics: if input doesn't match any available function.
            """
            # syntax with mu is not necessary so accomodate both
            if test_stat in ["qmu", "q"]:
                test_stat = "q"
            elif test_stat in ["qtilde", "qmutilde"]:
                test_stat = "qmutilde"
            options = {"qmutilde": qmu_tilde, "q": qmu}

            if options.get(test_stat, False) is False:
                raise UnknownTestStatistics(f"Requested test statistics {test_stat} does not exist.")

            return options[test_stat]

        def compute_teststatistics(mu, maximum_likelihood, logpdf, maximum_asimov_likelihood, asimov_logpdf, teststat):
            """
            Compute the test statistic for the observed data under the studied model.

            :param mu (`float`): Signal strength
            :param maximum_likelihood (`Tuple[float, float]`): muhat and minimum negative log-likelihood
            :param logpdf (`Callable[[float], float]`): log of the full density
            :param maximum_asimov_likelihood (`Tuple[float, float]`): muhat and minimum negative
                                                                      log-likelihood for asimov data
            :param asimov_logpdf (`Callable[[float], float]`): log of the full density for asimov data
            :param teststat (`Text`): `"qmutilde"`, `"q"` or `"q0"`
            :return `Tuple[float, float, float]`: sqrt(qmu), sqrt(qmuA) and distance between them
            :raises `UnknownTestStatistics`: if input doesn't match any available function.
            """
            teststat_func = get_test_statistic(teststat)

            muhat, min_nll = maximum_likelihood
            muhatA, min_nllA = maximum_asimov_likelihood

            # min_logpdf = -min_nll
            qmu = teststat_func(mu, muhat, -min_nll, logpdf)
            qmuA = teststat_func(mu, muhatA, -min_nllA, asimov_logpdf)
            sqrt_qmu = np.sqrt(qmu)
            sqrt_qmuA = np.sqrt(qmuA)

            if teststat in ["q", "q0", "qmu"]:
                delta_teststat = sqrt_qmu - sqrt_qmuA
            else:
                delta_teststat = (
                    sqrt_qmu - sqrt_qmuA
                    if sqrt_qmu <= sqrt_qmuA
                    else np.true_divide(qmu - qmuA, 2.0 * sqrt_qmuA)
                )

            return sqrt_qmu, sqrt_qmuA, delta_teststat

        test_stat = "q" if allow_negative_signal else "qmutilde"
        (maximum_likelihood, logpdf, maximum_asimov_likelihood, asimov_logpdf) = statModel._prepare_for_hypotest(expected=ExpectationType.observed, test_statistics=test_stat)
        sqrt_qmu_spey, sqrt_qmuA_spey, delta_teststat = compute_teststatistics(mu=1.,maximum_likelihood=maximum_likelihood,logpdf=logpdf ,maximum_asimov_likelihood=maximum_asimov_likelihood,asimov_logpdf=asimov_logpdf,teststat=test_stat)

        self.assertAlmostEqual(sqmu_SL,sqrt_qmu_spey,4)
        self.assertAlmostEqual(sqA_SL,sqrt_qmuA_spey,2)

    def testUpperLimitOnSigmaTimesEff_pyhf_1(self):
        """
        Test the cross-section upper limit computation with simplified likelhood backend from SModelS and Spey packages.
        """
        database='./testFilesTim/test_db/'
        # Set the path to the database
        database = Database(database)
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2019-08"])[0]

        filename = './testFilesTim/TChiWH_400_100_400_100.slha'
        model = Model(BSMList, SMList)
        model.updateParticles(filename)
        smstoplist = decomposer.decompose(model, sigmacut=0.)
        expected = False
        allow_negative_signal = False

        # pyhf from SModelS
        dataSetResults = []
        # Compute predictions for each data set (for UL analyses there is one single set)
        for dataset in expRes.datasets:
            predList = _getDataSetPredictions(dataset, smstoplist, maxMassDist=0.2)
            if predList:
                dataSetResults.append(predList)
        combinedRes = _getCombinedResultFor(dataSetResults, expRes, marginalize=False)
        dataset = combinedRes.dataset

        srNsigDict = dict([[pred.dataset.getID(),(pred.xsection.value * pred.dataset.getLumi()).asNumber(),] for pred in combinedRes.datasetPredictions])
        nsig = [srNsigDict[ds.getID()] if ds.getID() in srNsigDict else 0.0 for ds in dataset.origdatasets]

        nsignals = list()
        datasets = [ds.getID() for ds in dataset.origdatasets]
        for jsName in dataset.globalInfo.jsonFiles:
            subSig = list()
            for srName in dataset.globalInfo.jsonFiles[jsName]:
                sig = nsig[datasets.index(srName)]
                subSig.append(sig)
            nsignals.append(subSig)

        def getWSInfo(jsons):
            """
            Getting informations from the json files

            :ivar channelsInfo: list of dictionaries (one dictionary for each json file) containing useful information about the json files
                - :key signalRegions: list of dictonaries with 'json path' and 'size' (number of bins) of the 'signal regions' channels in the json files
                - :key otherRegions: list of strings indicating the path to the control and validation region channels
            """
            # Identifying the path to the SR and VR channels in the main workspace files
            channelsInfo = []  # workspace specifications
            for ws in dataset.globalInfo.jsons:
                wsChannelsInfo = {}
                wsChannelsInfo["signalRegions"] = []
                wsChannelsInfo["otherRegions"] = []
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
                channelsInfo.append(wsChannelsInfo)
            return channelsInfo

        def patchMaker(jsons, nsignals, includeCRs):
            """
            Method that creates the list of patches to be applied to the `self.inputJsons` workspaces, one for each region given the `self.nsignals` and the informations available in `self.channelsInfo` and the content of the `self.inputJsons`
            NB: It seems we need to include the change of the "modifiers" in the patches as well

            :return: the list of patches, one for each workspace
            """
            # Constructing the patches to be applied on the main workspace files
            patches = []
            channelsInfo = getWSInfo(jsons)
            for ws, info, subSig in zip(jsons, channelsInfo, nsignals):
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

                if not includeCRs:
                    for path in info["otherRegions"]:
                        patch.append({"op": "remove", "path": path})
                patches.append(patch)
            return patches

        includeCRs = False
        if hasattr(dataset.globalInfo, "includeCRs"):
            includeCRs = dataset.globalInfo.includeCRs

        patches = patchMaker(dataset.globalInfo.jsons, nsignals, includeCRs)

        def CLs_root_pyhf(mu, data, model, args):
            return 0.05 - float(hypotest(mu, data, model, **args))

        res_pyhf = {}
        for jsName, json, patch in zip(dataset.globalInfo.jsonFiles, dataset.globalInfo.jsons, patches):
            wsDict = apply_patch(json, patch)
            ws = Workspace(wsDict)
            msettings = {
                "normsys": {"interpcode": "code4"},
                "histosys": {"interpcode": "code4p"},
            }
            model = ws.model(modifier_settings=msettings)
            data = ws.data(model)

            low = 0.1
            high = 1.
            bounds = model.config.suggested_bounds()
            args = {}
            args["par_bounds"] = bounds
            while True:
                CLs_low = CLs_root_pyhf(low, data, model, args)
                CLs_high = CLs_root_pyhf(high, data, model, args)
                if CLs_low < 0 and CLs_high > 0:
                    break
                elif CLs_low < 0 and CLs_high < 0:
                    high *= 2.
                elif CLs_low > 0 and CLs_high > 0:
                    low *= 1/2.
                else:
                    print(f'CLs_low > 0 and CLs_high < 0: CLs_root_pyhf({low}) = {CLs_root_pyhf(low, data, model, args)}; CLs_root_pyhf({high})={CLs_root_pyhf(high, data, model, args)} for {jsName}')
                bounds[model.config.poi_index] = (0., high+0.1)
                args["par_bounds"] = bounds
            mu_ul_pyhf = optimize.brentq(lambda mu: CLs_root_pyhf(mu, data, model, args), low, high, rtol=1e-3, xtol=1e-3)
            mu_hat_pyhf, minNllh_pyhf = fit(data, model, return_fitted_val=True, par_bounds = bounds)
            xsec = combinedRes.xsection.value
            res_pyhf[jsName] = {"mu_ul":mu_ul_pyhf,"xsec_ul":mu_ul_pyhf*xsec, "mu_hat":mu_hat_pyhf[model.config.poi_index], "llhdMax":np.exp(-minNllh_pyhf/2.)}


        # Computation of statistical quantities with Spey
        xsec_ul_spey = getCombinedUpperLimitFor(dataset, nsig, expected=False, deltas_rel=combinedRes.deltas_rel, allowNegativeSignals=allow_negative_signal)
        res = getCombinedPyhfStatistics(dataset=dataset, nsig=nsig, marginalize=False, deltas_rel=combinedRes.deltas_rel, nll=False, expected=False, allowNegativeSignals=allow_negative_signal)
        mu_hat_spey = res["muhat"]
        llhdMax_spey = res["lmax"]

        self.assertAlmostEqual(res_pyhf['BkgOnly.json']['xsec_ul'].asNumber(fb),xsec_ul_spey.asNumber(fb),3)
        self.assertAlmostEqual(res_pyhf['BkgOnly.json']['mu_hat'],mu_hat_spey,4)
        self.assertAlmostEqual(res_pyhf['BkgOnly.json']['llhdMax'],llhdMax_spey,4)

    def testUpperLimitOnSigmaTimesEff_pyhf_2(self):
        """
        Test the cross-section upper limit computation with simplified likelhood backend from SModelS and Spey packages.
        """
        database='official'
        # Set the path to the database
        database = Database(database)
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2018-31"])[1]

        filename = './testFilesTim/T6bbHH_1000_600_60_1000_600_60.slha'
        filename = '/home/pascal/SModelS/smodels/test/testFilesTim/T6bbHH_1000_600_60_1000_600_60.slha'
        model = Model(BSMList, SMList)
        model.updateParticles(filename)
        smstoplist = decomposer.decompose(model, sigmacut=0.)
        expected = False
        allow_negative_signal = False

        # pyhf from SModelS
        dataSetResults = []
        # Compute predictions for each data set (for UL analyses there is one single set)
        for dataset in expRes.datasets:
            predList = _getDataSetPredictions(dataset, smstoplist, maxMassDist=0.2)
            if predList:
                dataSetResults.append(predList)
        combinedRes = _getCombinedResultFor(dataSetResults, expRes, marginalize=False)
        dataset = combinedRes.dataset

        srNsigDict = dict([[pred.dataset.getID(),(pred.xsection.value * pred.dataset.getLumi()).asNumber(),] for pred in combinedRes.datasetPredictions])
        nsig = [srNsigDict[ds.getID()] if ds.getID() in srNsigDict else 0.0 for ds in dataset.origdatasets]

        nsignals = list()
        datasets = [ds.getID() for ds in dataset.origdatasets]
        for jsName in dataset.globalInfo.jsonFiles:
            subSig = list()
            for srName in dataset.globalInfo.jsonFiles[jsName]:
                sig = nsig[datasets.index(srName)]
                subSig.append(sig)
            nsignals.append(subSig)

        def getWSInfo(jsons):
            """
            Getting informations from the json files

            :ivar channelsInfo: list of dictionaries (one dictionary for each json file) containing useful information about the json files
                - :key signalRegions: list of dictonaries with 'json path' and 'size' (number of bins) of the 'signal regions' channels in the json files
                - :key otherRegions: list of strings indicating the path to the control and validation region channels
            """
            # Identifying the path to the SR and VR channels in the main workspace files
            channelsInfo = []  # workspace specifications
            for ws in dataset.globalInfo.jsons:
                wsChannelsInfo = {}
                wsChannelsInfo["signalRegions"] = []
                wsChannelsInfo["otherRegions"] = []
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
                channelsInfo.append(wsChannelsInfo)
            return channelsInfo

        def patchMaker(jsons, nsignals, includeCRs):
            """
            Method that creates the list of patches to be applied to the `self.inputJsons` workspaces, one for each region given the `self.nsignals` and the informations available in `self.channelsInfo` and the content of the `self.inputJsons`
            NB: It seems we need to include the change of the "modifiers" in the patches as well

            :return: the list of patches, one for each workspace
            """
            # Constructing the patches to be applied on the main workspace files
            patches = []
            channelsInfo = getWSInfo(jsons)
            for ws, info, subSig in zip(jsons, channelsInfo, nsignals):
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

                if not includeCRs:
                    for path in info["otherRegions"]:
                        patch.append({"op": "remove", "path": path})
                patches.append(patch)
            return patches

        includeCRs = False
        if hasattr(dataset.globalInfo, "includeCRs"):
            includeCRs = dataset.globalInfo.includeCRs

        patches = patchMaker(dataset.globalInfo.jsons, nsignals, includeCRs)

        def CLs_root_pyhf(mu, data, model, args):
            return 0.05 - float(hypotest(mu, data, model, **args))

        def findBracketBounds(model, data, jsName, minimum_poi):
            if allow_negative_signal:
                low = minimum_poi/5
            else:
                low = 0.1
            high = 1.
            bounds = model.config.suggested_bounds()
            if allow_negative_signal:
                bounds[config.poi_index] = (config.minimum_poi, 100)
            else:
                bounds[config.poi_index] = (0, 100)
            args = {}
            args["par_bounds"] = bounds
            while True:
                CLs_low = CLs_root_pyhf(low, data, model, args)
                CLs_high = CLs_root_pyhf(high, data, model, args)
                if CLs_low < 0 and CLs_high > 0:
                    break
                elif CLs_low < 0 and CLs_high < 0:
                    high *= 2.
                elif CLs_low > 0 and CLs_high > 0:
                    if allow_negative_signal:
                        if low-0.5 < minimum_poi:
                            low = minimum_poi
                        else:
                            low = low-0.5
                    low *= 1/2.
                else:
                    print(f'CLs_low > 0 and CLs_high < 0: CLs_root_pyhf({low}) = {CLs_root_pyhf(low, data, model, args)}; CLs_root_pyhf({high})={CLs_root_pyhf(high, data, model, args)} for {jsName}')
                bounds[model.config.poi_index] = (low, high+0.1)
                args["par_bounds"] = bounds
            return low, high, args

        res_pyhf = {"observed":{},"apriori":{}}
        for expectationType in [ExpectationType.observed,ExpectationType.apriori]:
            expected = "observed" if expectationType==ExpectationType.observed else "apriori"
            for jsName, json, patch in zip(dataset.globalInfo.jsonFiles, dataset.globalInfo.jsons, patches):
                # if a patch has no signal, skip
                if all([sig==0. for sig in patch[0]['value']['data']]):
                    continue
                statModel = get_multi_region_statistical_model(analysis=dataset.globalInfo.id,
                                                                signal=patch,
                                                                observed=json,
                                                                xsection=sum(nsig)/dataset.getLumi()
                                                                )
                config = statModel.backend.model.config()
                bounds = config.suggested_bounds
                if allow_negative_signal:
                    bounds[config.poi_index] = (config.minimum_poi, 100)
                else:
                    bounds[config.poi_index] = (0, 100)
                wsDict = apply_patch(json, patch)
                if expected == "apriori":
                    for id,sample in enumerate(wsDict['channels'][0]['samples']):
                        if sample['name']=="Bkg":
                            index_bkg = id
                    for obs in wsDict['observations']:
                        if obs['name']==wsDict['channels'][0]['name']:
                            obs['data']=wsDict['channels'][0]['samples'][index_bkg]['data']

                msettings = {
                    "normsys": {"interpcode": "code4"},
                    "histosys": {"interpcode": "code4p"},
                }
                ws = Workspace(wsDict)
                model = ws.model(modifier_settings=msettings)
                data = ws.data(model)
                low, high, args = findBracketBounds(model,data,jsName,config.minimum_poi)
                # print(f'CLs_root_pyhf({low}) = {CLs_root_pyhf(low, data, model, args)}; CLs_root_pyhf({high})={CLs_root_pyhf(high, data, model, args)} for {jsName}')

                mu_ul_pyhf = optimize.brentq(lambda mu: CLs_root_pyhf(mu, data, model, args), low, high, rtol=1e-3, xtol=1e-3)
                xsec = combinedRes.xsection.value.asNumber(fb)*fb
                xsec_ul_pyhf = mu_ul_pyhf*xsec

                mu_hat_pyhf, minNllh_pyhf = fit(data, model, return_fitted_val=True, par_bounds = bounds)
                llhdMax_pyhf = np.exp(-minNllh_pyhf/2.)
                mu_hat_pyhf = mu_hat_pyhf[model.config.poi_index]

                res_pyhf[expected][jsName] = {"mu_ul":mu_ul_pyhf,"xsec_ul":xsec_ul_pyhf, "mu_hat":mu_hat_pyhf, "llhdMax":llhdMax_pyhf}


                mu_ul_spey = statModel.poi_upper_limit(expected=expectationType,allow_negative_signal=allow_negative_signal,par_bounds=bounds)
                xsec_ul_spey = mu_ul_spey*statModel.xsection

                mu_hat_spey, llhdMax_spey = statModel.maximize_likelihood(allow_negative_signal=allow_negative_signal, expected=expectationType, return_nll=False)

                # For prefit dataset, mu_hat must be 0
                if expected == 'apriori':
                    self.assertAlmostEqual(mu_hat_pyhf,0.,3)
                    self.assertAlmostEqual(mu_hat_spey,0.,3)
                self.assertAlmostEqual(mu_ul_pyhf*xsec._value,xsec_ul_spey._value,3)
                self.assertAlmostEqual(mu_hat_pyhf,mu_hat_spey,4)
                self.assertAlmostEqual(llhdMax_pyhf,llhdMax_spey,4)

        # Computation of statistical quantities with Spey
        xsec_ul_spey = getCombinedUpperLimitFor(dataset, nsig, expected=False, deltas_rel=combinedRes.deltas_rel, allowNegativeSignals=allow_negative_signal)
        res = getCombinedPyhfStatistics(dataset=dataset, nsig=nsig, marginalize=False, deltas_rel=combinedRes.deltas_rel, nll=False, expected=False, allowNegativeSignals=allow_negative_signal)
        mu_hat_spey = res["muhat"]
        llhdMax_spey = res["lmax"]

        # Find the "best" set of SRs, i.e. the one that gives the lowest expected upper limit on mu
        bestJSName = None
        mu_hat_exp_min = np.inf
        for JSName in res_pyhf['apriori']:
            if res_pyhf['apriori'][JSName]['mu_ul'] < mu_hat_exp_min:
                mu_hat_exp_min = res_pyhf['apriori'][JSName]['mu_ul']
                bestJSName = JSName

        self.assertAlmostEqual(res_pyhf['observed'][bestJSName]['xsec_ul']._value,xsec_ul_spey._value,3)
        self.assertAlmostEqual(res_pyhf['observed'][bestJSName]['mu_hat'],mu_hat_spey,4)
        self.assertAlmostEqual(res_pyhf['observed'][bestJSName]['llhdMax'],llhdMax_spey,4)


    def testUnderfluctuatingLlhdsFromLimits(self):
        """test the likelihoods from limits allowing for underfluctuations"""
        comparisons = {
            False: {0: 0.2869456402153424, 5: 0.05697393370613237},
            True: {0: 0.43038397431221154, 5: 0.030362667100323704},
        }
        print("Truncated Gaussian is not a spey feature yet.")
        for nsig in [0, 5]:
            for allowNegatives in [False, True]:
                computer = TruncatedGaussians ( 4.5, 5.45, nsig )
                ret = computer.likelihood ( mu=1.,
                       nll = False, allowNegativeMuhat = allowNegatives )
                llhdlim, muhat, sigma_mu = ret["llhd"], ret["muhat"], ret["sigma_mu"]
                c = comparisons[allowNegatives][nsig]
                self.assertAlmostEqual(llhdlim, c, 2)

    def testCorrectedLlhdsFromLimits(self):
        """test the likelihoods from limits allowing for underfluctuations"""
        comparisons = {
            0.0: {0: 0.094661, 3: 0.151640, 5: 0.12555219},
            0.6: {0: 0.1233638, 3: 0.1504459, 5: 0.11380},
        }
        print("Truncated Gaussian is not a spey feature yet.")
        for nsig in [0, 3, 5]:
            for x in [0.0, 0.6]:
                computer = TruncatedGaussians ( 8.52, 6.18, nsig, corr = x )
                ret = computer.likelihood ( 1., False, False )
                llhdlim, muhat, sigma_mu = ret["llhd"], ret["muhat"], ret["sigma_mu"]
                c = comparisons[x][nsig]
                self.assertAlmostEqual(llhdlim, c, 2)

    def testChi2FromLimits(self):
        """test the chi2 value that we obtain from limits"""
        allow_negative_signal = False
        print("Truncated Gaussian is not a spey feature yet.")
        for backendNumber in [1,2]:
            #backendNumber == 1: pyhf backend
            #backendNumber == 2: simplified likelihood backend
            nsig = 35.0
            nobs, nbg = 110., 100.0
            lumi = 1.
            statModel = get_uncorrelated_region_statistical_model(observations=nobs,
                                                                    backgrounds=nbg,
                                                                    background_uncertainty=np.sqrt(0.001),
                                                                    signal_yields=nsig,
                                                                    xsection=None,
                                                                    analysis="UnitTest",
                                                                    backend=AvailableBackends(backendNumber)
                                                                )
            config = statModel.backend.model.config()
            bounds = config.suggested_bounds
            if allow_negative_signal:
                bounds[config.poi_index] = (config.minimum_poi, 100)
            else:
                bounds[config.poi_index] = (0, 100)
            mu_ul_obs = statModel.poi_upper_limit(expected=ExpectationType.observed, par_bounds=bounds)
            mu_ul_exp = statModel.poi_upper_limit(expected=ExpectationType.aposteriori, par_bounds=bounds)
            xsec = nsig/lumi
            xsec_ul_obs = xsec*mu_ul_obs
            xsec_ul_exp = xsec*mu_ul_exp
            llhddir = statModel.likelihood(poi_test=1.,return_nll=False)
            chi2dir = statModel.chi2()
            ### Marginalization does not work but we do not care for the moment.
            # if backendNumber == 2: # pyhf cannot marginalize
            #     llhdmarg = statModel.likelihood(mu=1., marginalize=True)
            #     chi2marg = statModel.chi2(marginalize=True)
            computer = TruncatedGaussians ( xsec_ul_obs, xsec_ul_exp, nsig )
            ret = computer.likelihood ( mu=1. )
            llhdlim, muhat, sigma_mu = ret["llhd"], ret["muhat"], ret["sigma_mu"]
            self.assertAlmostEqual(llhdlim,0.0034205732477661462,4) # Allowing differences below 5 places was a bit too restrictive, ok with only 4 (0.2% difference is ok)
            self.assertAlmostEqual(muhat,0.23328649242374602,3)
            self.assertAlmostEqual(sigma_mu,0.3383372145700653,3) # Allowing differences below 4 places was a bit too restrictive, ok with only 3 (0.08% difference is ok)
            ### Marginalization does not work but we do not care for the moment.
            # if backendNumber == 2: # pyhf cannot marginalize?
            #     chi2lim = computer.chi2 ( ) # llhdlim )
            #     ## relative error on chi2, for this example is about 4%
            #     rel = abs(chi2lim - chi2marg) / chi2marg
            #     self.assertAlmostEqual(rel, 0.04, 1)

    def testUpperLimit(self):
        allow_negative_signal = True
        for backendNumber in [1,2]:
            #backendNumber == 1: pyhf backend
            #backendNumber == 2: simplified likelihood backend
            nsig = 1.0
            nobs, nbg = 100.0, 100.0
            bgError = round(np.sqrt(0.001),5) # pyhf convergence fails for small values of the background error, so fix it at a small value that converges correclty.
            statModel = get_uncorrelated_region_statistical_model(observations=nobs,
                                                                    backgrounds=nbg,
                                                                    background_uncertainty=bgError,
                                                                    signal_yields=nsig,
                                                                    xsection=None,
                                                                    analysis="UnitTest",
                                                                    backend=AvailableBackends(backendNumber)
                                                                )
            config = statModel.backend.model.config()
            bounds = config.suggested_bounds
            if allow_negative_signal:
                bounds[config.poi_index] = (config.minimum_poi, 100)
            else:
                bounds[config.poi_index] = (0, 100)
            re = statModel.poi_upper_limit(par_bounds=bounds)
            self.assertAlmostEqual(re/(1.06*20.), 1., 1)

    def testApproxGaussian(self):
        print("Truncated Gaussian is not a spey feature yet.")
        ## turn experimental features on
        from smodels.tools import runtime

        runtime._experimental = True
        expRes = database.getExpResults(analysisIDs=["CMS-PAS-SUS-12-026"])
        self.assertTrue(len(expRes), 1)
        filename = "./testFiles/slha/T1tttt.slha"
        model = Model(BSMList, SMList)
        model.updateParticles(filename)
        smstoplist = decomposer.decompose(model, sigmacut=0)
        prediction = theoryPredictionsFor(expRes[0], smstoplist)[0]
        prediction.computeStatistics()

        c = 0.0
        for muval in np.arange(0.0, 0.2, 0.02):
            llhd = prediction.likelihood(mu=muval)
            c += llhd
        self.assertAlmostEqual(prediction.likelihood(), 1.563288e-35, 3)
        self.assertAlmostEqual(c, 0.011523436957977766, 3)

    def testPredictionInterface(self):
        """A simple test to see that the interface in datasetObj
        and TheoryPrediction to the statistics tools is working correctly
        """
        expRes = database.getExpResults(analysisIDs=["CMS-SUS-13-012"])[0]

        filename = "./testFiles/slha/simplyGluino.slha"
        model = Model(BSMList, SMList)
        model.updateParticles(filename)
        smstoplist = decomposer.decompose(model, sigmacut=0)
        prediction = theoryPredictionsFor(expRes, smstoplist)[0]
        pred_signal_strength = prediction.xsection.value
        nsig = (pred_signal_strength * expRes.globalInfo.lumi).asNumber()

        for backendNumber,backend in zip([1,2],["pyhf","SL"]):
            #backendNumber == 1: pyhf backend
            #backendNumber == 2: simplified likelihood backend
            prediction.computeStatistics(backend=backend)
            ill = math.log(prediction.likelihood())
            ichi2 = prediction.chi2()
            statModel = get_uncorrelated_region_statistical_model(observations=4.,
                                                                    backgrounds=2.2,
                                                                    background_uncertainty=1.1,
                                                                    signal_yields=nsig,
                                                                    xsection=None,
                                                                    analysis="UnitTest",
                                                                    backend=AvailableBackends(backendNumber)
                                                                )
            dll = math.log(statModel.likelihood(poi_test=1.,return_nll=False))
            self.assertTrue(abs((ill-dll)*100./dll) < 3) # Previous test "self.assertAlmostEqual(ill, dll, places=2)" was too restrictive.
                                                         # 3% difference is ok when comparing results from simplified_likelihood backend only (using SModelS) to results obtained with pyhf backend from spey.
                                                         # 0.01% difference is ok when comparing results from simplified_likelihood backend only (using SModelS) to results obtained with simplified_backend backend from spey.
            dchi2 = statModel.chi2()
            # print ( "dchi2,ichi2",dchi2,ichi2)
            self.assertTrue(abs((ichi2-dchi2)*100./ichi2) < 1) # Previous test "self.assertAlmostEqual(ichi2, dchi2, places=2)" was too restrictive.
                                                               # 1% difference is ok when comparing results from simplified_likelihood backend only (using SModelS) to results obtained with pyhf backend from spey.
                                                               # 0.2% difference is ok when comparing results from simplified_likelihood backend only (using SModelS) to results obtained with simplified_backend backend from spey.

    def testZeroLikelihood(self):
        """A test to check if a llhd of 0 is being tolerated"""
        for backendNumber in [1,2]:
            #backendNumber == 1: pyhf backend
            #backendNumber == 2: simplified likelihood backend
            statModel = get_uncorrelated_region_statistical_model(observations=1e20,
                                                                    backgrounds=2.2,
                                                                    background_uncertainty=1.1,
                                                                    signal_yields=2.,
                                                                    xsection=None,
                                                                    analysis="UnitTest",
                                                                    backend=AvailableBackends(backendNumber)
                                                                )
            # Spey returns nan for pyhf backend but 0 for SL backend
            llhd = statModel.likelihood(poi_test=1.,return_nll=False)
            self.assertAlmostEqual(0.0, llhd, places=2)
            try:
                dchi2 = statModel.chi2()
                # Spey gives a negative chi2 for the moment!
                # ichi2 = 4.486108149972863e21
                # self.assertAlmostEqual(dchi2 / ichi2, 1.0, places=4)
            except IndexError as e:
                self.assertTrue("too many indices for array: array is 0-dimensional, but 1 were indexed" in str(e))

    def round_to_sign(self, x, sig=3):
        """
        Round the given number to the significant number of digits.
        """
        if x == 0.0:
            return 0.0
        rounding = sig - int(floor(log10(x))) - 1
        if rounding == 0:
            return int(x)
        else:
            return round(x, rounding)

    def testLikelihood(self):
        """
        Compare the computed chi2 from a given observed
        and expected upper limit and a theory prediction
        with the previously known result for the value of
        the chi2.


        All values of nobs, nsig, nb, deltab come from the
        SModelS database and are for the T1 simplified model
        from efficiency results of one of
        ATLAS-SUSY-2013-02
        ATLAS-CONF-2013-047
        CMS-SUS-13-012
        ATLAS-CONF-2013-054

        """

        expected_values = [
            # mgluino          mlsp          nsig               nobs              nb             deltab           llhd                 chi2
            # ----------       ----------    ---------------    ----------------  ----------     ----------       -------------------- ----------------
            {
                "mgluino": 500,
                "mlsp": 200,
                "nsig": 384.898,
                "nobs": 298.413,
                "nb": 111.0,
                "deltab": 11.0,
                "llhd": 0.00024197,
                "chi2": 7.32614,
            }
        ]
        {
            "mgluino": 500,
            "mlsp": 300,
            "nsig": 185.166,
            "nobs": 223.619,
            "nb": 111.0,
            "deltab": 11.0,
            "llhd": 0.00215989,
            "chi2": 3.67088900,
        },
        {
            "mgluino": 500,
            "mlsp": 400,
            "nsig": 450.820,
            "nobs": 2331.38,
            "nb": 2120.0,
            "deltab": 110.0,
            "llhd": 0.00075499,
            "chi2": 2.85026,
        },
        {
            "mgluino": 600,
            "mlsp": 100,
            "nsig": 476.150,
            "nobs": 437.874,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00100575,
            "chi2": 3.52406595,
        },
        {
            "mgluino": 600,
            "mlsp": 200,
            "nsig": 264.387,
            "nobs": 232.912,
            "nb": 111.0,
            "deltab": 11.0,
            "llhd": 0.00028921,
            "chi2": 7.57051432,
        },
        {
            "mgluino": 600,
            "mlsp": 300,
            "nsig": 171.766,
            "nobs": 211.038,
            "nb": 111.0,
            "deltab": 11.0,
            "llhd": 0.00198061,
            "chi2": 4.02903,
        },
        {
            "mgluino": 600,
            "mlsp": 400,
            "nsig": 66.9991,
            "nobs": 150.393,
            "nb": 111.0,
            "deltab": 11.0,
            "llhd": 0.00845030,
            "chi2": 1.89194013,
        },
        {
            "mgluino": 600,
            "mlsp": 500,
            "nsig": 157.571,
            "nobs": 2167.25,
            "nb": 2120.0,
            "deltab": 110.0,
            "llhd": 0.00217371,
            "chi2": 0.845615,
        },
        {
            "mgluino": 700,
            "mlsp": 100,
            "nsig": 307.492,
            "nobs": 325.060,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00159632,
            "chi2": 3.41955,
        },
        {
            "mgluino": 700,
            "mlsp": 200,
            "nsig": 211.534,
            "nobs": 228.763,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00061670,
            "chi2": 6.26852955,
        },
        {
            "mgluino": 700,
            "mlsp": 300,
            "nsig": 147.084,
            "nobs": 167.631,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00012707,
            "chi2": 10.1408114,
        },
        {
            "mgluino": 700,
            "mlsp": 400,
            "nsig": 420.524,
            "nobs": 2332.28,
            "nb": 2120.0,
            "deltab": 110.0,
            "llhd": 0.00100854,
            "chi2": 2.28195399,
        },
        {
            "mgluino": 700,
            "mlsp": 500,
            "nsig": 186.726,
            "nobs": 2162.70,
            "nb": 2120.0,
            "deltab": 110.0,
            "llhd": 0.00165411,
            "chi2": 1.39601,
        },
        {
            "mgluino": 700,
            "mlsp": 600,
            "nsig": 5.18888,
            "nobs": 24.3271,
            "nb": 37.0,
            "deltab": 6.0,
            "llhd": 0.00545866,
            "chi2": 4.3836,
        },
        {
            "mgluino": 800,
            "mlsp": 100,
            "nsig": 169.670,
            "nobs": 213.312,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00116298,
            "chi2": 5.09803029,
        },
        {
            "mgluino": 800,
            "mlsp": 200,
            "nsig": 152.221,
            "nobs": 212.732,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00223053,
            "chi2": 3.81519907,
        },
        {
            "mgluino": 800,
            "mlsp": 300,
            "nsig": 98.6749,
            "nobs": 175.141,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00298021,
            "chi2": 3.72566221,
        },
        {
            "mgluino": 800,
            "mlsp": 400,
            "nsig": 59.3935,
            "nobs": 141.966,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00262008,
            "chi2": 4.30322736,
        },
        {
            "mgluino": 800,
            "mlsp": 500,
            "nsig": 27.7738,
            "nobs": 123.172,
            "nb": 111.0,
            "deltab": 11.0,
            "llhd": 0.01605069,
            "chi2": 0.903864,
        },
        {
            "mgluino": 800,
            "mlsp": 600,
            "nsig": 6.40339,
            "nobs": 39.2979,
            "nb": 33.0,
            "deltab": 6.0,
            "llhd": 0.04536718,
            "chi2": -0.000501411,
        },
        {
            "mgluino": 800,
            "mlsp": 700,
            "nsig": 4.38635,
            "nobs": 132.824,
            "nb": 125.0,
            "deltab": 10.0,
            "llhd": 0.02525385,
            "chi2": 0.0565348,
        },
        {
            "mgluino": 900,
            "mlsp": 100,
            "nsig": 18.8255,
            "nobs": 14.1228,
            "nb": 4.9,
            "deltab": 1.6,
            "llhd": 0.02122262,
            "chi2": 2.85426343,
        },
        {
            "mgluino": 900,
            "mlsp": 200,
            "nsig": 16.0543,
            "nobs": 6.77062,
            "nb": 4.9,
            "deltab": 1.6,
            "llhd": 0.00187567,
            "chi2": 8.43890579,
        },
        {
            "mgluino": 900,
            "mlsp": 300,
            "nsig": 64.4188,
            "nobs": 142.220,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00181163,
            "chi2": 5.00642964,
        },
        {
            "mgluino": 900,
            "mlsp": 400,
            "nsig": 44.8312,
            "nobs": 140.979,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00692173,
            "chi2": 2.34800741,
        },
        {
            "mgluino": 900,
            "mlsp": 500,
            "nsig": 24.4723,
            "nobs": 120.688,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00601021,
            "chi2": 2.72454478,
        },
        {
            "mgluino": 900,
            "mlsp": 600,
            "nsig": 67.0446,
            "nobs": 2165.25,
            "nb": 2120.0,
            "deltab": 110.0,
            "llhd": 0.00328372,
            "chi2": 0.0371125,
        },
        {
            "mgluino": 900,
            "mlsp": 700,
            "nsig": 1.00167,
            "nobs": 5.0,
            "nb": 5.2,
            "deltab": 1.4,
            "llhd": 0.13964962,
            "chi2": 0.107139,
        },
        {
            "mgluino": 900,
            "mlsp": 800,
            "nsig": 0.86634,
            "nobs": 24.0,
            "nb": 37.0,
            "deltab": 6.0,
            "llhd": 0.01303119,
            "chi2": 2.638323,
        },
        {
            "mgluino": 1000,
            "mlsp": 100,
            "nsig": 11.7426,
            "nobs": 11.8786,
            "nb": 4.9,
            "deltab": 1.6,
            "llhd": 0.05712388,
            "chi2": 1.06934,
        },
        {
            "mgluino": 1000,
            "mlsp": 200,
            "nsig": 9.85815,
            "nobs": 7.98535,
            "nb": 4.9,
            "deltab": 1.6,
            "llhd": 0.03180710,
            "chi2": 2.63593288,
        },
        {
            "mgluino": 1000,
            "mlsp": 300,
            "nsig": 6.80275,
            "nobs": 6.14772,
            "nb": 4.9,
            "deltab": 1.6,
            "llhd": 0.04255251,
            "chi2": 2.25703866,
        },
        {
            "mgluino": 1000,
            "mlsp": 400,
            "nsig": 25.8451,
            "nobs": 120.523,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.00525390,
            "chi2": 2.99137140,
        },
        {
            "mgluino": 1000,
            "mlsp": 500,
            "nsig": 18.6299,
            "nobs": 122.095,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.01056408,
            "chi2": 1.59516468,
        },
        {
            "mgluino": 1000,
            "mlsp": 600,
            "nsig": 10.2636,
            "nobs": 119.968,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.01536934,
            "chi2": 0.84011292,
        },
        {
            "mgluino": 1000,
            "mlsp": 700,
            "nsig": 4.59470,
            "nobs": 121.728,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.02078487,
            "chi2": 0.23333618,
        },
        {
            "mgluino": 1000,
            "mlsp": 800,
            "nsig": 1.91162,
            "nobs": 121.196,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.02193946,
            "chi2": 0.12552973,
        },
        {
            "mgluino": 1000,
            "mlsp": 900,
            "nsig": 0.62255,
            "nobs": 133.0,
            "nb": 125.0,
            "deltab": 10.0,
            "llhd": 0.02290147,
            "chi2": 0.253293,
        },
        {
            "mgluino": 1100,
            "mlsp": 100,
            "nsig": 5.95636,
            "nobs": 5.94515,
            "nb": 4.9,
            "deltab": 1.6,
            "llhd": 0.05236744,
            "chi2": 1.87080349,
        },
        {
            "mgluino": 1100,
            "mlsp": 200,
            "nsig": 5.51938,
            "nobs": 5.92472,
            "nb": 4.9,
            "deltab": 1.6,
            "llhd": 0.06002489,
            "chi2": 1.59429,
        },
        {
            "mgluino": 1100,
            "mlsp": 300,
            "nsig": 3.93082,
            "nobs": 6.07873,
            "nb": 4.9,
            "deltab": 1.6,
            "llhd": 0.09732480,
            "chi2": 0.61881103,
        },
        {
            "mgluino": 1100,
            "mlsp": 400,
            "nsig": 2.80428,
            "nobs": 6.54033,
            "nb": 4.9,
            "deltab": 1.6,
            "llhd": 0.12350249,
            "chi2": 0.0882501,
        },
        {
            "mgluino": 1100,
            "mlsp": 500,
            "nsig": 12.6778,
            "nobs": 125.271,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.01749246,
            "chi2": 0.560614,
        },
        {
            "mgluino": 1100,
            "mlsp": 600,
            "nsig": 8.02475,
            "nobs": 119.742,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.01700322,
            "chi2": 0.64005829,
        },
        {
            "mgluino": 1100,
            "mlsp": 700,
            "nsig": 4.74108,
            "nobs": 120.211,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.01979279,
            "chi2": 0.334838,
        },
        {
            "mgluino": 1100,
            "mlsp": 800,
            "nsig": 1.79622,
            "nobs": 120.858,
            "nb": 126.0,
            "deltab": 13.0,
            "llhd": 0.02187799,
            "chi2": 0.134012,
        },
        {
            "mgluino": 1100,
            "mlsp": 900,
            "nsig": 4.82397,
            "nobs": 2166.20,
            "nb": 2120.0,
            "deltab": 110.0,
            "llhd": 0.00313424,
            "chi2": 0.119045,
        },
        {
            "mgluino": 1100,
            "mlsp": 1000,
            "nsig": 0.1606,
            "nobs": 25.0,
            "nb": 37.0,
            "deltab": 6.0,
            "llhd": 0.01796058,
            "chi2": 1.9806,
        },
        {
            "mgluino": 2100,
            "mlsp": 1000,
            "nsig": 0.1606,
            "nobs": 2.0,
            "nb": 0.7,
            "deltab": 6.0,
            "llhd": 0.108669,
            "chi2": -0.161304,
        }

        for d in expected_values:
            nobs = d["nobs"]
            nsig = d["nsig"]
            nb = d["nb"]
            deltab = d["deltab"]
            # print ("ns="+str(nsig)+"; nobs = "+str(nobs)+"; nb="+str(nb)+"; db="+str(deltab))
            # Chi2 as computed by statistics module:
            for backendNumber in [1,2]:
                #backendNumber == 1: pyhf backend
                #backendNumber == 2: simplified likelihood backend
                statModel = get_uncorrelated_region_statistical_model(observations=nobs,
                                                                        backgrounds=nb,
                                                                        background_uncertainty=deltab,
                                                                        signal_yields=2.,
                                                                        xsection=None,
                                                                        analysis="UnitTest",
                                                                        backend=AvailableBackends(backendNumber)
                                                                    )
                # spey doesn't allow for marginalized chi2 yet
                # chi2_actual = statModel.chi2( marginalize=True)  # , .2*nsig )
                # chi2_expected = d["chi2"]
                # print('chi2exp', chi2_expected)
                # if not chi2_expected == None and not np.isnan(chi2_expected):
                #     #                 chi2_expected = self.round_to_sign(chi2_expected, 2)
                #     # Check that chi2 values agree:
                #     self.assertAlmostEqual(
                #         abs(chi2_actual - chi2_expected) / chi2_expected, 0.0, places=2
                #     )
                # else:
                #     self.assertTrue(chi2_actual == None or np.isnan(chi2_actual))

                # likelihood as computed by statistics module:
                # computer = LikelihoodComputer( nobs, nb, deltab**2 )
                # likelihood_actual = statistics.likelihood( nsig,
                #    nobs, nb, deltab, deltas)
                likelihood_actual = statModel.likelihood(poi_test=1.,return_nll=False)
                # likelihood_actual = statistics.likelihood()
                #             logger.error("llk= "+str(likelihood_actual)+" nsig="+str(nsig)+" nobs = "+str(nobs)+" nb="+str(nb)+"+-"+str(deltab))
                # print('llhdactual', likelihood_actual)
                if not likelihood_actual == None and not np.isnan(likelihood_actual):
                    likelihood_actual = self.round_to_sign(likelihood_actual, 4)

                # The previously computed likelihood:
                # (using: ntoys=100000)
                likelihood_expected = d["llhd"]
                # print('llhdexp', likelihood_expected)
                if not likelihood_expected == None and not np.isnan(likelihood_expected):
                    likelihood_expected = self.round_to_sign(likelihood_expected, 4)

                    # Check that likelihood values agree:
                    self.assertAlmostEqual(likelihood_actual, likelihood_expected, delta=2 * 1e-1)
                else:
                    self.assertTrue(likelihood_actual == None or np.isnan(likelihood_actual))


if __name__ == "__main__":
    unittest.main()
