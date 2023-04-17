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

import time

class SpeyTest(unittest.TestCase):
    def test21002 ( self ):
        db = Database ( "./dbSUS21002/" )
        er = db.getExpResults ()[0]
        defaults = {}
        ### these are the ones given by CMS-SUS-21-002
        # defaults [ "TChiWZ_600_200_600_200" ] = { "obs": 1.0705, "exp": 1.0711 }
        # defaults [ "TChiWH_600_300_600_300" ] = { "obs": 0.4825, "exp": 0.7151 }
        # defaults [ "TChiWH_1000_350_1000_350" ] = { "obs": 0.6222, "exp": 0.9412 }
        # defaults [ "TChiWZ_900_1_900_1" ] = { "obs": 0.6315, "exp": 0.9771 }
        ### this is what we get
        defaults [ "TChiWZ_600_200_600_200" ] = { "obs": 1.141, "exp": 1.163 }
        defaults [ "TChiWH_600_300_600_300" ] = { "obs": 0.5076, "exp": 0.7924 }
        defaults [ "TChiWH_1000_350_1000_350" ] = { "obs": 0.55713, "exp": 0.8857 }
        defaults [ "TChiWZ_900_1_900_1" ] = { "obs": 0.6109, "exp": 0.9915 }
        verbose = True
        t0 = time.time()
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
                    if verbose:
                        print ( f"{exp}: r {r[exp]:.3f} r_base {base[exp]} delta {delta:.3f}" )
                    if delta > .1:
                        line = f"mismatch for {slhaname}({exp}): base={base[exp]}, computed={r[exp]}"
                        logger.error ( line )
                        self.assertTrue ( delta < .1 )
        t = time.time() - t0
        if verbose:
            print ( f"took {t:.3f}s" )
        

if __name__ == "__main__":
    unittest.main()
