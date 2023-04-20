#!/usr/bin/env python3

"""
.. module:: testSpey
   :synopsis: Test spey against the numbers we got from CMS-SUS-21-002

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
        
import unittest
import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.particlesLoader import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from smodels.theory import decomposer
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
import numpy as np

import time

class SpeyTest(unittest.TestCase):
    def test21002 ( self ):
        db = Database ( "./dbSUS21002/" )
        er = db.getExpResults ()[0]
        defaults = {}
        ### this is what we get
        defaults [ "TChiWZ_600_200_600_200" ] = { "obs": 1.141, "exp": 1.163 }
        defaults [ "TChiWH_600_300_600_300" ] = { "obs": 0.5076, "exp": 0.7924 }
        defaults [ "TChiWH_1000_350_1000_350" ] = { "obs": 0.55713, "exp": 0.8857 }
        defaults [ "TChiWZ_900_1_900_1" ] = { "obs": 0.6109, "exp": 0.9915 }
        if True:
        ### these are the ones given by CMS-SUS-21-002
            defaults [ "TChiWZ_600_200_600_200" ] = { "obs": 1.0705, "exp": 1.0711 }
            defaults [ "TChiWH_600_300_600_300" ] = { "obs": 0.4825, "exp": 0.7151 }
            defaults [ "TChiWH_1000_350_1000_350" ] = { "obs": 0.6222, "exp": 0.9412 }
            defaults [ "TChiWZ_900_1_900_1" ] = { "obs": 0.6315, "exp": 0.9771 }
        verbose = True
        t0 = time.time()
        deltas = []
        for slhaname in defaults.keys():
            model = Model(BSMparticles=BSMList, SMparticles=SMList)
            fname = f"./testFiles/21002/{slhaname}.slha"
            model.updateParticles(inputFile=fname)
            toplist = decomposer.decompose(model)
            predictions = theoryPredictionsFor(er, toplist, combinedResults=True)
            for p in predictions:
                robs,rexp = ( p.getRValue( expected=x) for x in [ False, True ] )
                r = { "obs": p.getRValue(), "exp": p.getRValue ( expected=True ) }
                base = defaults[slhaname]
                for exp in [ "obs", "exp" ]:
                    delta = 2. * abs ( r[exp] - base[exp] ) / ( r[exp]+base[exp] )
                    deltas.append ( delta )
                    if verbose:
                        print ( f"{slhaname} {exp}: r {r[exp]:.3f} r_base {base[exp]} delta {delta:.3f}" )
                    if delta > .15:
                        line = f"mismatch for {slhaname}({exp}): base={base[exp]}, computed={r[exp]}"
                        logger.error ( line )
                        self.assertTrue ( delta < .15 )
        t = time.time() - t0
        if verbose:
            print ( f"took {t:.3f}s, average delta {np.mean(deltas):.4f}" )
        
    def testSingleRegion ( self ):
        from smodels.tools.speyTools import SpeyComputer, SimpleSpeyDataSet
        nobs,bg,bgerr,lumi,reful,refule = 3905, 3658.3, 238.767, 35.9/fb, 21.5336356*fb, 13.455494*fb
        cases = [ { "nobs": 3905, "bg": 3658.3, "bgerr": 238.767, "lumi": 35.9/fb,
                    "reful": 21.5336356*fb, "refule": 13.455494*fb } ]
        cases.append ( { "nobs": 0, "bg": .001, "bgerr": .01, "lumi": 35.9/fb,
                         "reful": 0.053678853586*fb, "refule": 0.0537152307*fb } )
        cases.append ( { "nobs": 3, "bg": 4.1, "bgerr": .6533758489, "lumi": 35.9/fb,
                         "reful": 0.123728081503*fb, "refule": 0.1520938597094264*fb } )
        for case in cases:
            dataset = SimpleSpeyDataSet ( case["nobs"], case["bg"], case["bgerr"], 
                                          case["lumi"] )
            computer = SpeyComputer ( dataset, 1. )
            ul = computer.poi_upper_limit ( expected = False, limit_on_xsec = True )
            ule = computer.poi_upper_limit ( expected = True, limit_on_xsec = True )
            if False:
                print ( "ul", ul.asNumber(fb) )
                print ( "ule", ule.asNumber(fb) )
            self.assertAlmostEqual ( abs((ul - case["reful"]).asNumber(fb)), 0., 3 )
            self.assertAlmostEqual ( abs((ule - case["refule"]).asNumber(fb)), 0.,3 )
        


if __name__ == "__main__":
    unittest.main()
