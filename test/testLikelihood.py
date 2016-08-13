#!/usr/bin/env python

"""
.. module: testLikelihood
   :synopsis: Test likelihood functions module.

.. moduleauthor:: Jory Sonneveld <jory@opmijnfiets.nl>

"""

import sys
sys.path.append('..')
import unittest
from smodels.tools import likelihood as like
import math
import numpy
import redirector

numpy.seterr ( all="ignore" )



class TestLikelihood(unittest.TestCase):
    """
    Unittest class to test likelihood functions.
    """

    def mest_mean_std(self):
        """
        Test if mean and standard deviation found
        are values that in turn would yield
        0.95 when testing the confidence level
        as given by error functions.
        See 1202.3415.
        This is check is built into the likelihood
        module.
        """
        ## print("\ntest_mean_std")
        exp_cl = 0.95

        expected = (None, None)
        actual = like.mean_std(0, 0)
        self.assertEqual(expected, actual)

        mean, std = like.mean_std(1, None, lumi=20) # fb, fb-1
        if not mean == None:
            actual_cl = round(like._cl_from_erfs(mean, std, 1), 2)
            self.assertTrue(exp_cl - 0.01 <= actual_cl <= exp_cl + 0.01)

        mean, std = like.mean_std(1, None, lumi=20) # fb, fb-1
        if not mean == None:
            actual_cl = round(like._cl_from_erfs(mean, std, 1/1.05), 2)
            self.assertTrue(exp_cl - 0.01 <= actual_cl <= exp_cl + 0.01)

        mean, std = like.mean_std(1, exp=1.96, lumi=20000) # pb, pb-1
        actual_cl = round(like._cl_from_erfs(mean, std, 1), 2)
        self.assertTrue(exp_cl - 0.01 <= actual_cl <= exp_cl + 0.01)


    def mest_cl_from_erf(self):
        """
        Test if _cl_from_erfs gives results as expected from
        error functions.
        """
        ## print("\ntest_cl_from_erf")
        expected = math.erf(1)
        actual = like._cl_from_erfs(0, 1/math.sqrt(2), 1)
        self.assertEqual(expected, actual)

        expected = math.erf(0.5)
        actual = like._cl_from_erfs(0, 1/math.sqrt(2), 0.5)
        self.assertEqual(expected, actual)

        denom = math.sqrt(2)
        expected = (math.erf((0.8-0.5)/denom) + math.erf(0.5/denom))/\
                (1+math.erf(0.5/denom))
        actual = like._cl_from_erfs(0.5, 1, 0.8)
        self.assertEqual(expected, actual)



    def mest_real_vs_estimated_mean(self):
        """
        Test in reverse: from a certain expectation of
        behavior of the likelihood function to its
        result.

        Specifically, assume expected are 200 events,
        observed are 220, mu is around half so theory
        expects 40 events, and the lumi is 20 inverse
        femtobarn. Then the likelihood should peak
        around mu is half.

        Another example:
        If one assumes expected are 200 events,
        observed are 220, mu is around 1 so theory
        expects 20 events, and the lumi is 20 inverse
        femtobarn. Then the likelihood should peak
        around 20 events also, or mu = 1.

        Note that the error tolerance is high, mu +- 0.1*mu.
        """
        ## print("\ntest_real_mean_vs_estimated_mean")
        b = 200 # background
        n = 220 # observed
        lumi = 20.0 # fb-1
        s = 40 # theory prediction for signal
        ## expect mu_max (mean) is half, since events seen is
        ## halfway between background (null hyp) and theory signal

        sigma_expected = like._simple_signal_calculator(b, b)/lumi
        sigma_observed = like._simple_signal_calculator(b, n)/lumi

        mean, std = like.mean_std(sigma_observed, sigma_expected, lumi=lumi)
        if mean == None:
            print('Mean was None for b, n, s, lumi:', b, n, s, lumi)
        expected_mean, expected_std = like._mean_std_direct(lumi, b, n)
        actual_mean = mean
        # check that computed mean is within twenty precent of expected mean:
        self.assertAlmostEqual(expected_mean, actual_mean, delta=0.1*expected_mean)

    def mest_unit_handling(self):
        """
        Test the unit handling in the likelihood module.
        Compute a chi2 using upper limits and luminosity
        carrying equal (except luminosity) units.

        The likelihood computation uses mathematics modules in scipy
        and has to proceed with floats.
        For this reason units are stripped off the given quantities
        before computing likelihoods.
        Important is that everything is in the same unit (pb and pb-1,
        or fb and fb-1, etc.).
        """
        ## print("test_unit_handling\n")

        from smodels.tools.physicsUnits import fb, pb, GeV, TeV
        theo =  5.72E-01 * pb # theory prediction
        obs =  2.42E-01 * pb # observed upper limit
        exp =  1.21E-01 * pb # expected upper limit
        lumi = 20.3 / fb # luminosity
        chi2_actual = like.chi2(theo=theo, obs=obs, exp=exp, lumi=lumi)

        # Strip units
        theo = theo.asNumber(pb)
        obs = obs.asNumber(pb)
        exp = exp.asNumber(pb)
        lumi = lumi.asNumber(1/pb)
        chi2_expected = like.chi2(theo=theo, obs=obs, exp=exp, lumi=lumi)

        self.assertAlmostEqual(chi2_actual, chi2_expected, places=4)


    def test_different_unit_handling(self):
        """
        Test the unit handling in the likelihood module.
        Compute a chi2 using upper limits and luminosity
        of varying units.

        The likelihood computation uses mathematics modules in scipy
        and has to proceed with floats.
        For this reason units are stripped off the given quantities
        before computing likelihoods.
        Important is that everything is in the same unit (pb and pb-1,
        or fb and fb-1, etc.).
        """
        # print("test_different_unit_handling\n")

        from smodels.tools.physicsUnits import fb, pb, GeV, TeV
        theo =  5.72E-01 * pb # theory prediction
        obs =  2.42E+02 * fb # observed upper limit
        exp =  1.21E+02 * fb # expected upper limit
        lumi = 20.3 / fb # luminosity
        with redirector.stderr_redirected():
            chi2_actual = like.chi2(theo=theo, obs=obs, exp=exp, lumi=lumi)

        # Strip units
        theo = theo.asNumber(pb)
        obs = obs.asNumber(pb)
        exp = exp.asNumber(pb)
        lumi = lumi.asNumber(1/pb)
        chi2_expected = like.chi2(theo=theo, obs=obs, exp=exp, lumi=lumi)

        self.assertAlmostEqual(chi2_actual, chi2_expected, places=4)




    def test_computes_as_expected(self):
        """
        Compare the computed chi2 from a given observed
        and expected upper limit and a theory prediction
        with the previously known result for the value of
        the chi2.

        All values come from CMS-SUS-13-012 and are for the
        SUSY T2 simplified model.

        """

        # Some computed chi2 for the cMSSM per m0, m12 mass:
        cmssm = {(360, 320): {'theo': 0.159063876816674, 'ul':
                0.0241183906043336, 'chi2': 396.58185129759, 'chi2same':
                167.091726448567, 'exp': 0.0144994229179552},
            (220, 400): {'theo': 0.000107010958268822,
                'ul': 0.0156837974746417, 'chi2': 2.64747547067593,
                'chi2same': 0.00017738124919616, 'exp': 0.00926164467441977},
            (320, 420): {'theo': 0.0468244069320068, 'ul': 0.0133231969697019,
                'chi2': 101.760215892775, 'chi2same': 47.4496400369738,
                'exp': 0.00779852391849953},
            (680, 340): {'theo': 0.00262701906839619, 'ul': 0.0112503806366831,
                'chi2': 0.91178864361843, 'chi2same': 0.209411432544899,
                'exp': 0.00645267360013492},
            (300, 320): {'theo': 0.192082854327662, 'ul': 0.028065328061141,
                'chi2': 390.51781228915, 'chi2same': 179.94739230636,
                'exp': 0.0177837644860081},
            (740, 300): {'theo': 0.000940339604199469,
                'ul': 0.0114050240268197, 'chi2': 2.17881032984429,
                'chi2same': 0.0260972810558525, 'exp': 0.00653143569784816}}

        # Some computed chi2 for SUSY T2 per squark and LSP mass:
        susyt2 = {(920, 240): {'theo': 0.0074949587,
                'ul': 0.014833663776516, 'chi2': 0.002623592,
                'chi2same': 0.98063145, 'exp': 0.008875729143619},
            (600, 200): {'theo': 0.2022680048, 'ul': 0.094371646642684,
                'chi2': 22.444613026, 'chi2same': 17.647051507,
                'exp': 0.070692755281925},
            (580, 200): {'theo': 0.2448195834, 'ul': 0.111287331581116,
                'chi2': 23.092230949, 'chi2same': 18.590966972,
                'exp': 0.085996641218662},
            (780, 640): {'theo': 0.0287726581, 'ul': 0.481117188930511,
                'chi2': 2.243092366, 'chi2same': 0.013726668,
                'exp': 0.280160522460937},
            (920, 60): {'theo': 0.0073324571, 'ul': 0.01252143420279,
                'chi2': 0.063742229, 'chi2same': 1.317230221,
                'exp': 0.00720074782148},
            (760, 300): {'theo': 0.0351058525, 'ul': 0.043850513547658,
                'chi2': 1.08977845, 'chi2same': 2.462019886,
                'exp': 0.027029524743557},
            (300, 20): {'theo': 14.398782, 'ul': 1.02101625204086,
                'chi2': 786.312657357, 'chi2same': 764.005750521,
                'exp': 1.00471384525299}}

        for known_values in [cmssm, susyt2]:
            for key in known_values:

                di = known_values[key]
                # The observed upper limit:
                observed_ul = di['ul']
                # The expected upper limit:
                expected_ul = di['exp']
                # The given theory prediction:
                theo = di['theo']
                # The previously computed chi2:
                chi2_expected = di['chi2']
                # The previously computed chi2 for expected_ul = observed_ul:
                # (This is used in the absence of expected upper limits)
                chi2same_expected = di['chi2same']

                # Compute the chi2:
                with redirector.stderr_redirected( ):
                    with redirector.stdout_redirected ():
                        chi2_actual = like.chi2(theo, observed_ul, expected_ul)
                        chi2same_actual = like.chi2(theo, observed_ul, observed_ul)

                self.assertAlmostEqual(chi2_actual, chi2_expected, places=4)
                self.assertAlmostEqual(chi2same_actual, chi2same_expected,
                    places=4)

if __name__ == '__main__':
    unittest.main(exit=False)
