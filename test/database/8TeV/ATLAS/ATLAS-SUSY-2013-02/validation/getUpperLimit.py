#!/usr/bin/python

import sys,os
home = os.path.expanduser("~")
sys.path.insert(0,os.path.join(home,"smodels-utils/"))
sys.path.insert(0,os.path.join(home,"smodels/"))

from smodels.experiment.databaseObjects import Database
from validation.plotProducer import ValidationPlot, getExpIdFromPath
from smodels.tools.physicsUnits import pb, fb, GeV
from smodels_utils.dataPreparation import origPlotObjects
NAN=float('nan')

axes="2*Eq(mother,x)_Eq(inter0,0.5*x+0.5*y)_Eq(lsp,y)"
origP=origPlotObjects.OrigPlot.fromString(axes)


database = Database(os.path.join(home,"smodels-database"))
expRes = database.getExpResults(analysisIDs=[ getExpIdFromPath() ],datasetIDs=[None])

masses=[]
masses.append ( [[600*GeV, 300*GeV, 70*GeV],[ 600*GeV, 300*GeV, 70*GeV ]] )
masses.append ( [[600*GeV, 400*GeV, 100*GeV],[ 600*GeV, 400*GeV, 100*GeV ]] )
masses.append ( [[580*GeV, 330*GeV, 80*GeV],[ 580*GeV, 330*GeV, 80*GeV ]] )
masses.append ( [[585.*GeV, 332.5*GeV, 80.*GeV],[ 585.*GeV, 332.5*GeV, 80.*GeV ]] )
masses.append ( [[585.*GeV, 332.*GeV, 80.*GeV],[ 585.*GeV, 332.*GeV, 80.*GeV ]] )

unitless = origP.getParticleMasses ( 585.0, 80.0 )
masses.append ( [[m*GeV for m in mm] for mm in unitless ] )

for mass in masses:
    print mass, expRes.getUpperLimitFor(txname="T5WW",mass=mass )
