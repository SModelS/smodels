#!/usr/bin/env python3

"""
.. module:: testStatistics
   :synopsis: Tests the statistics functionality.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Jory Sonneveld <jory@opmijnfiets.nl>

"""
import sys
sys.path.insert(0,"../")
import unittest
#from smodels.tools import statistics
from smodels.tools.simplifiedLikelihoods import UpperLimitComputer, LikelihoodComputer, Data
from smodels.tools.statistics import likelihoodFromLimits, chi2FromLimits
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from databaseLoader import database
from smodels.theory import decomposer
from math import floor, log10
import numpy as np
import math

class StatisticsTest(unittest.TestCase):

    def lLHDFromLimits(self):
        """ to do some statistics on the chi2 """
        nsig = 1.0
        nobs, nbg = 100, 100.
        m = Data( nobs, nbg, 0.001, None, nsig, deltas_rel=0.)
        ulcomp = UpperLimitComputer()
        ulobs = ulcomp.ulSigma(m)
        ulexp = ulcomp.ulSigma(m,expected=True)
        print ( "ulobs", ulobs )
        print ( "ulexp", ulexp )
        f=open("llhds.csv","wt")
        dx=.5
        totdir,totlim,totmarg=0.,0.,0.
        for nsig in np.arange ( 0.1, 100., dx ):
            print ( )
            print ( "nsig=", nsig )
            m = Data( nobs, nbg, 0.001, None, nsig, deltas_rel=0.)
            llhdcomp = LikelihoodComputer ( m )
            llhddir = llhdcomp.likelihood ( nsig )
            chi2dir = llhdcomp.chi2 ( nsig )
            llhdmarg = llhdcomp.likelihood ( nsig, marginalize=True )
            chi2marg = llhdcomp.chi2 ( nsig, marginalize=True )
            print ( "llhd direct", llhddir, chi2dir )
            print ( "llhd marg", llhdmarg, chi2marg )
            llhdlim = likelihoodFromLimits ( ulobs, ulexp, nsig )
            chi2lim = chi2FromLimits ( llhdlim, ulexp )
            print ( "llhd from limits", llhdlim, chi2lim )
            totdir+=llhddir*dx
            totlim+=llhdlim*dx
            totmarg+=llhdmarg*dx
            f.write ( "%s,%s,%s,%s\n" % ( nsig, llhddir, llhdlim, llhdmarg ) )
        print ( "total direct", totdir )
        print ( "total limit", totlim )
        print ( "total marg", totmarg )
        f.close()

    def testChi2FromLimits(self):
        """ test the chi2 value that we obtain from limits """
        nsig = 35.0
        nobs, nbg = 110, 100.
        m = Data( nobs, nbg, 0.001, None, nsig, deltas_rel=0.)
        ulcomp = UpperLimitComputer()
        ulobs = ulcomp.ulSigma(m)
        ulexp = ulcomp.ulSigma(m,expected=True)
        dx=.5
        m = Data( nobs, nbg, 0.001, None, nsig, deltas_rel=0.)
        llhdcomp = LikelihoodComputer ( m )
        llhddir = llhdcomp.likelihood ( nsig )
        chi2dir = llhdcomp.chi2 ( nsig )
        llhdmarg = llhdcomp.likelihood ( nsig, marginalize=True )
        chi2marg = llhdcomp.chi2 ( nsig, marginalize=True )
        llhdlim = likelihoodFromLimits ( ulobs, ulexp, nsig )
        chi2lim = chi2FromLimits ( llhdlim, ulexp )
        ## relative error on chi2, for this example is about 4% 
        rel = abs (chi2lim - chi2marg ) / chi2marg
        self.assertAlmostEqual ( rel, 0.04, 1 )

    def testUpperLimit(self):
        m = Data( 100., 100., 0.001, None, 1.0,deltas_rel=0.)
        comp = UpperLimitComputer()
        re = comp.ulSigma(m)
        self.assertAlmostEqual(re/(1.06*20.),1., 1)

    def testApproxGaussian(self):
        ## turn experimental features on
        from smodels.tools import runtime
        runtime._experimental = True
        expRes = database.getExpResults(analysisIDs=["CMS-PAS-SUS-12-026"])
        self.assertTrue(len(expRes),1)
        filename = "./testFiles/slha/T1tttt.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename)
        smstoplist = decomposer.decompose(model, sigmacut=0)
        prediction = theoryPredictionsFor(expRes[0], smstoplist)[0]
        prediction.computeStatistics()
        import numpy
        c=0.
        for i in numpy.arange(0.,.2 ,.02 ):
            l=prediction.getLikelihood(i)
            c+=l
        self.assertAlmostEqual ( prediction.likelihood, 1.563288e-35, 3 )
        self.assertAlmostEqual ( c, 0.010213526179269492, 3 )

    def testPredictionInterface(self):
        """ A simple test to see that the interface in datasetObj
        and TheoryPrediction to the statistics tools is working correctly
        """
        expRes = database.getExpResults(analysisIDs=['CMS-SUS-13-012'] )[0]

        filename = "./testFiles/slha/simplyGluino.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename)
        smstoplist = decomposer.decompose(model, sigmacut=0)
        prediction = theoryPredictionsFor(expRes, smstoplist)[0]
        pred_signal_strength = prediction.xsection.value
        prediction.computeStatistics()
        ill = math.log(prediction.likelihood)
        ichi2 = prediction.chi2
        nsig = (pred_signal_strength*expRes.globalInfo.lumi).asNumber()
        m = Data(4, 2.2, 1.1**2, None, nsignal=nsig,deltas_rel=0.2)
        computer = LikelihoodComputer(m)
        dll = math.log(computer.likelihood(nsig, marginalize=False ) )
        self.assertAlmostEqual(ill, dll, places=2)
        dchi2 = computer.chi2( nsig, marginalize=False )
        # print ( "dchi2,ichi2",dchi2,ichi2)
        self.assertAlmostEqual(ichi2, dchi2, places=2)

    def testZeroLikelihood(self):
        """ A test to check if a llhd of 0 is being tolerated
        """
        nsig = 2
        m = Data(1e20, 2.2, 1.1**2, None, nsignal=nsig, deltas_rel=0.2 )
        computer = LikelihoodComputer(m)
        llhd = computer.likelihood(2, marginalize=False )
        nll = computer.likelihood(2, marginalize=False, nll=True )
        self.assertAlmostEqual(0., llhd, places=2)
        dchi2 = computer.chi2( nsig, marginalize=False )
        ichi2= 4.486108149972863e+21
        self.assertAlmostEqual(ichi2, dchi2, places=2)

    def round_to_sign(self, x, sig=3):
        """
        Round the given number to the significant number of digits.
        """
        if x==0.: return 0.
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
        {'mgluino':  500, 'mlsp':  200, 'nsig':   384.898, 'nobs':  298.413, 'nb':  111.0 , 'deltab': 11.0 , 'llhd': 0.00024197 , 'chi2': 7.32614 }]
        {'mgluino':  500, 'mlsp':  300, 'nsig':   185.166, 'nobs':  223.619, 'nb':  111.0 , 'deltab': 11.0 , 'llhd': 0.00215989 , 'chi2': 3.67088900 },
        {'mgluino':  500, 'mlsp':  400, 'nsig':   450.820, 'nobs':  2331.38, 'nb':  2120.0, 'deltab': 110.0, 'llhd': 0.00075499 , 'chi2': 2.85026 },
        {'mgluino':  600, 'mlsp':  100, 'nsig':   476.150, 'nobs':  437.874, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00100575 , 'chi2': 3.52406595 },
        {'mgluino':  600, 'mlsp':  200, 'nsig':   264.387, 'nobs':  232.912, 'nb':  111.0 , 'deltab': 11.0 , 'llhd': 0.00028921 , 'chi2': 7.57051432 },
        {'mgluino':  600, 'mlsp':  300, 'nsig':   171.766, 'nobs':  211.038, 'nb':  111.0 , 'deltab': 11.0 , 'llhd': 0.00198061 , 'chi2': 4.02903 },
        {'mgluino':  600, 'mlsp':  400, 'nsig':   66.9991, 'nobs':  150.393, 'nb':  111.0 , 'deltab': 11.0 , 'llhd': 0.00845030 , 'chi2': 1.89194013 },
        {'mgluino':  600, 'mlsp':  500, 'nsig':   157.571, 'nobs':  2167.25, 'nb':  2120.0, 'deltab': 110.0, 'llhd': 0.00217371 , 'chi2': 0.845615 },
        {'mgluino':  700, 'mlsp':  100, 'nsig':   307.492, 'nobs':  325.060, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00159632 , 'chi2': 3.41955 },
        {'mgluino':  700, 'mlsp':  200, 'nsig':   211.534, 'nobs':  228.763, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00061670 , 'chi2': 6.26852955 },
        {'mgluino':  700, 'mlsp':  300, 'nsig':   147.084, 'nobs':  167.631, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00012707 , 'chi2': 10.1408114 },
        {'mgluino':  700, 'mlsp':  400, 'nsig':   420.524, 'nobs':  2332.28, 'nb':  2120.0, 'deltab': 110.0, 'llhd': 0.00100854 , 'chi2': 2.28195399 },
        {'mgluino':  700, 'mlsp':  500, 'nsig':   186.726, 'nobs':  2162.70, 'nb':  2120.0, 'deltab': 110.0, 'llhd': 0.00165411 , 'chi2': 1.39601 },
        {'mgluino':  700, 'mlsp':  600, 'nsig':   5.18888, 'nobs':  24.3271, 'nb':  37.0  , 'deltab': 6.0  , 'llhd': 0.00545866 , 'chi2': 4.3836     },
        {'mgluino':  800, 'mlsp':  100, 'nsig':   169.670, 'nobs':  213.312, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00116298 , 'chi2': 5.09803029 },
        {'mgluino':  800, 'mlsp':  200, 'nsig':   152.221, 'nobs':  212.732, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00223053 , 'chi2': 3.81519907 },
        {'mgluino':  800, 'mlsp':  300, 'nsig':   98.6749, 'nobs':  175.141, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00298021 , 'chi2': 3.72566221 },
        {'mgluino':  800, 'mlsp':  400, 'nsig':   59.3935, 'nobs':  141.966, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00262008 , 'chi2': 4.30322736 },
        {'mgluino':  800, 'mlsp':  500, 'nsig':   27.7738, 'nobs':  123.172, 'nb':  111.0 , 'deltab': 11.0 , 'llhd': 0.01605069 , 'chi2': 0.903864 },
        {'mgluino':  800, 'mlsp':  600, 'nsig':   6.40339, 'nobs':  39.2979, 'nb':  33.0  , 'deltab': 6.0  , 'llhd': 0.04536718 , 'chi2': -0.000501411 },
        {'mgluino':  800, 'mlsp':  700, 'nsig':   4.38635, 'nobs':  132.824, 'nb':  125.0 , 'deltab': 10.0 , 'llhd': 0.02525385 , 'chi2': 0.0565348 },
        {'mgluino':  900, 'mlsp':  100, 'nsig':   18.8255, 'nobs':  14.1228, 'nb':  4.9   , 'deltab': 1.6  , 'llhd': 0.02122262 , 'chi2': 2.85426343 },
        {'mgluino':  900, 'mlsp':  200, 'nsig':   16.0543, 'nobs':  6.77062, 'nb':  4.9   , 'deltab': 1.6  , 'llhd': 0.00187567 , 'chi2': 8.43890579 },
        {'mgluino':  900, 'mlsp':  300, 'nsig':   64.4188, 'nobs':  142.220, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00181163 , 'chi2': 5.00642964 },
        {'mgluino':  900, 'mlsp':  400, 'nsig':   44.8312, 'nobs':  140.979, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00692173 , 'chi2': 2.34800741 },
        {'mgluino':  900, 'mlsp':  500, 'nsig':   24.4723, 'nobs':  120.688, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00601021 , 'chi2': 2.72454478 },
        {'mgluino':  900, 'mlsp':  600, 'nsig':   67.0446, 'nobs':  2165.25, 'nb':  2120.0, 'deltab': 110.0, 'llhd': 0.00328372 , 'chi2': 0.0371125 },
        {'mgluino':  900, 'mlsp':  700, 'nsig':   1.00167, 'nobs':  5.0    , 'nb':  5.2   , 'deltab': 1.4  , 'llhd': 0.13964962 , 'chi2': 0.107139   },
        {'mgluino':  900, 'mlsp':  800, 'nsig':   0.86634, 'nobs':  24.0   , 'nb':  37.0  , 'deltab': 6.0  , 'llhd': 0.01303119 , 'chi2': 2.638323   },
        {'mgluino': 1000, 'mlsp':  100, 'nsig':   11.7426, 'nobs':  11.8786, 'nb':  4.9   , 'deltab': 1.6  , 'llhd': 0.05712388 , 'chi2': 1.06934 },
        {'mgluino': 1000, 'mlsp':  200, 'nsig':   9.85815, 'nobs':  7.98535, 'nb':  4.9   , 'deltab': 1.6  , 'llhd': 0.03180710 , 'chi2': 2.63593288 },
        {'mgluino': 1000, 'mlsp':  300, 'nsig':   6.80275, 'nobs':  6.14772, 'nb':  4.9   , 'deltab': 1.6  , 'llhd': 0.04255251 , 'chi2': 2.25703866 },
        {'mgluino': 1000, 'mlsp':  400, 'nsig':   25.8451, 'nobs':  120.523, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.00525390 , 'chi2': 2.99137140 },
        {'mgluino': 1000, 'mlsp':  500, 'nsig':   18.6299, 'nobs':  122.095, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.01056408 , 'chi2': 1.59516468 },
        {'mgluino': 1000, 'mlsp':  600, 'nsig':   10.2636, 'nobs':  119.968, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.01536934 , 'chi2': 0.84011292 },
        {'mgluino': 1000, 'mlsp':  700, 'nsig':   4.59470, 'nobs':  121.728, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.02078487 , 'chi2': 0.23333618 },
        {'mgluino': 1000, 'mlsp':  800, 'nsig':   1.91162, 'nobs':  121.196, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.02193946 , 'chi2': 0.12552973 },
        {'mgluino': 1000, 'mlsp':  900, 'nsig':   0.62255, 'nobs':  133.0  , 'nb':  125.0 , 'deltab': 10.0 , 'llhd': 0.02290147 , 'chi2': 0.253293 },
        {'mgluino': 1100, 'mlsp':  100, 'nsig':   5.95636, 'nobs':  5.94515, 'nb':  4.9   , 'deltab': 1.6  , 'llhd': 0.05236744 , 'chi2': 1.87080349 },
        {'mgluino': 1100, 'mlsp':  200, 'nsig':   5.51938, 'nobs':  5.92472, 'nb':  4.9   , 'deltab': 1.6  , 'llhd': 0.06002489 , 'chi2': 1.59429 },
        {'mgluino': 1100, 'mlsp':  300, 'nsig':   3.93082, 'nobs':  6.07873, 'nb':  4.9   , 'deltab': 1.6  , 'llhd': 0.09732480 , 'chi2': 0.61881103 },
        {'mgluino': 1100, 'mlsp':  400, 'nsig':   2.80428, 'nobs':  6.54033, 'nb':  4.9   , 'deltab': 1.6  , 'llhd': 0.12350249 , 'chi2': 0.0882501 },
        {'mgluino': 1100, 'mlsp':  500, 'nsig':   12.6778, 'nobs':  125.271, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.01749246 , 'chi2': 0.560614 },
        {'mgluino': 1100, 'mlsp':  600, 'nsig':   8.02475, 'nobs':  119.742, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.01700322 , 'chi2': 0.64005829 },
        {'mgluino': 1100, 'mlsp':  700, 'nsig':   4.74108, 'nobs':  120.211, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.01979279 , 'chi2': 0.334838 },
        {'mgluino': 1100, 'mlsp':  800, 'nsig':   1.79622, 'nobs':  120.858, 'nb':  126.0 , 'deltab': 13.0 , 'llhd': 0.02187799 , 'chi2': 0.134012 },
        {'mgluino': 1100, 'mlsp':  900, 'nsig':   4.82397, 'nobs':  2166.20, 'nb':  2120.0, 'deltab': 110.0, 'llhd': 0.00313424 , 'chi2': 0.119045 },
        {'mgluino': 1100, 'mlsp': 1000, 'nsig':   0.1606 , 'nobs':  25.0 ,   'nb':  37.0  , 'deltab': 6.0 ,  'llhd': 0.01796058 , 'chi2': 1.9806},
        {'mgluino': 2100, 'mlsp': 1000, 'nsig':   0.1606 , 'nobs':  2.0 ,   'nb':  0.7  , 'deltab': 6.0 ,  'llhd': 0.108669 , 'chi2': -0.161304}




        for d in expected_values:
            nobs = d['nobs']
            nsig = d['nsig']
            nb = d['nb']
            deltab = d['deltab']
            # print ("ns="+str(nsig)+"; nobs = "+str(nobs)+"; nb="+str(nb)+"; db="+str(deltab))
            # Chi2 as computed by statistics module:
            m = Data(nobs, nb, deltab**2,deltas_rel=0.2)
            computer = LikelihoodComputer(m)
            chi2_actual = computer.chi2(nsig, marginalize=True ) ## , .2*nsig )
            chi2_expected = d['chi2']
            #print('chi2exp', chi2_expected)
            if not chi2_expected==None and not np.isnan(chi2_expected):
#                 chi2_expected = self.round_to_sign(chi2_expected, 2)
                # Check that chi2 values agree:
                self.assertAlmostEqual(abs(chi2_actual-chi2_expected)/chi2_expected,0., places=2 )
            else:
                self.assertTrue(chi2_actual == None or np.isnan(chi2_actual))


            # likelihood as computed by statistics module:
            # computer = LikelihoodComputer( nobs, nb, deltab**2 )
            #likelihood_actual = statistics.likelihood( nsig,
            #    nobs, nb, deltab, deltas)
            likelihood_actual = computer.likelihood(nsig, marginalize=False )
            # likelihood_actual = statistics.likelihood()
#             logger.error("llk= "+str(likelihood_actual)+" nsig="+str(nsig)+" nobs = "+str(nobs)+" nb="+str(nb)+"+-"+str(deltab))
            #print('llhdactual', likelihood_actual)
            if not likelihood_actual==None and not np.isnan(likelihood_actual):
                likelihood_actual = self.round_to_sign(likelihood_actual, 4)

            # The previously computed likelihood:
            # (using: ntoys=100000)
            likelihood_expected = d['llhd']
            #print('llhdexp', likelihood_expected)
            if not likelihood_expected==None and not np.isnan(likelihood_expected):
                likelihood_expected = self.round_to_sign(likelihood_expected, 4)

                # Check that likelihood values agree:
                self.assertAlmostEqual(likelihood_actual, likelihood_expected,
                        delta=2*1e-1)
            else:
                self.assertTrue(likelihood_actual == None or np.isnan(likelihood_actual))

if __name__ == "__main__":
    unittest.main()
