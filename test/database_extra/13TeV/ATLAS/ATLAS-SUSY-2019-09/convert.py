#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse

argparser = argparse.ArgumentParser(description =
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath',
help = 'path to the package smodels_utils',\
type = str )
argparser.add_argument ('-smodelsPath', '--smodelsPath',
help = 'path to the package smodels_utils',\
type = str )
args = argparser.parse_args()

if args.utilsPath:
    utilsPath = args.utilsPath
else:
    databaseRoot = '../../../'
    sys.path.append(os.path.abspath(databaseRoot))
    from utilsPath import utilsPath
    utilsPath = databaseRoot + utilsPath
if args.smodelsPath:
    sys.path.append(os.path.abspath(args.smodelsPath))

sys.path.append(os.path.abspath(utilsPath))
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z



#+++++++ global info block ++++++++++++++
info 			= MetaInfoInput('ATLAS-SUSY-2019-09')
info.url 		= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/'
info.sqrts 		= 13
info.lumi 		= 139
info.publication = "https://link.springer.com/content/pdf/10.1140/epjc/s10052-021-09749-7.pdf"
info.prettyName = '3 leptons EW-ino'
info.private 	= False
info.arxiv 		= 'https://arxiv.org/abs/2106.01676'
info.contact    = 'atlas-phys-susy-conveners@cern.ch'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)


#+++++++ next txName block ++++++++++++++
TChiWZ                      = dataset.addTxName('TChiWZ')
TChiWZ.checked 				= 'no'
TChiWZ.constraint           = "[[['W']],[['Z']]]"
TChiWZ.condition            = None
TChiWZ.massConstraint       = None
TChiWZ.conditionDescription = None
TChiWZ.source               = "ATLAS"

# off-shell part

TChiWZoff             = dataset.addTxName('TChiWZoff')
TChiWZoff.checked     = 'no'
TChiWZoff.constraint  = "71.*([[['mu+','mu-']],[['l','nu']]] + [[['e+','e-']],[['l','nu']]])"
TChiWZoff.condition   = "cGtr([[['mu+','mu-']],[['l','nu']]],[[['e+','e-']],[['l','nu']]])"
TChiWZoff.massConstraint = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.source = "ATLAS"

#+++++++ next txName block ++++++++++++++
TChiWH                      = dataset.addTxName('TChiWH')
TChiWH.checked 				= 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.condition            = None
TChiWH.massConstraint       = None
TChiWH.conditionDescription = None
TChiWH.source               = "ATLAS"

#+++++++ next mass plane block ++++++++++++++ - For mass plane x,y
TChiWZ_plus_on                  = TChiWZ.addMassPlane(2*[[x, y]])
TChiWZ_plus_on.figure           = 'Aux. Fig. 6a'
TChiWZ_plus_on.figureUrl        = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/figaux_06a.png'
TChiWZ_plus_on.dataUrl          = 'https://www.hepdata.net/record/ins1866951?table=AuxFig%208a%20WZ%20Excl.%20Upper%20Limit%20Obs.%20Wino-bino(%2b)%20($\Delta%20m$)'
TChiWZ_plus_on.exclusionDataUrl = 'https://www.hepdata.net/record/ins1866951?version=1&table=Fig%2016b%20WZ%20Exclusion%3A%20Wino-bino(%2B)%20(%24%5CDelta%20m%24)%2C%20Obs'
TChiWZ_plus_on.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'upperLimits', 'expectedUpperLimits' ],
					units 		= [None, None, 'pb','pb'],
					dataFiles 	= ['orig/Fig16aWZExclusion:Wino-bino(+),onshell_Exp.csv','orig/Fig16aWZExclusion:Wino-bino(+),onshell_Obs.csv','orig/newAuxFig8aWZExcl.UpperLimitObs.Wino-bino(+)(Deltam).csv','orig/newAuxFig8bWZExcl.UpperLimitExp.Wino-bino(+)(Deltam).csv'],
					dataFormats	= ['csv','csv','csv','csv'])

# TChiWZ_plus_on.setSources(dataLabels 	= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits', 'expectedUpperLimits' ],
# 					units 		= [None, None, None, None, None, None, 'pb','pb'],
# 					dataFiles 	= ['orig/Fig16aWZExclusion:Wino-bino(+),Exp.csv','orig/Fig16aWZExclusion:Wino-bino(+),Exp_Down.csv','orig/Fig16aWZExclusion:Wino-bino(+),Exp_Up.csv','orig/Fig16aWZExclusion:Wino-bino(+),Obs.csv','orig/Fig16aWZExclusion:Wino-bino(+),Obs_Down.csv','orig/Fig16aWZExclusion:Wino-bino(+),Obs_Up.csv','orig/newAuxFig8aWZExcl.UpperLimitObs.Wino-bino(+)(Deltam).csv', 'orig/newAuxFig8bWZExcl.UpperLimitExp.Wino-bino(+)(Deltam).csv'],
# 					dataFormats	= ['csv','csv','csv','csv','csv','csv','csv','csv'])


#+++++++ next mass plane block ++++++++++++++ - For mass plane x,x-y
TChiWZ_plus_on                  = TChiWZ.addMassPlane(2*[[x, x-y]])
TChiWZ_plus_on.figure           = 'Aux. Fig. 6a'
TChiWZ_plus_on.figureUrl        = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/figaux_06a.png'
TChiWZ_plus_on.dataUrl          = 'https://www.hepdata.net/record/ins1866951?table=AuxFig%208a%20WZ%20Excl.%20Upper%20Limit%20Obs.%20Wino-bino(%2b)%20($\Delta%20m$)'
TChiWZ_plus_on.exclusionDataUrl = 'https://www.hepdata.net/record/ins1866951?version=1&table=Fig%2016b%20WZ%20Exclusion%3A%20Wino-bino(%2B)%20(%24%5CDelta%20m%24)%2C%20Obs'
TChiWZ_plus_on.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'upperLimits', 'expectedUpperLimits' ],
					units 		= [None, None, 'pb','pb'],
					dataFiles 	= ['orig/newFig16aWZExclusion:Wino-bino(+),onshell_Exp.csv','orig/newFig16aWZExclusion:Wino-bino(+),onshell_Obs.csv','orig/AuxFig8aWZExcl.UpperLimitObs.Wino-bino(+)(Deltam).csv','orig/AuxFig8bWZExcl.UpperLimitExp.Wino-bino(+)(Deltam).csv'],
					dataFormats	= ['csv','csv','csv','csv'])

# TChiWZ_plus_on.setSources(dataLabels 	= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits', 'expectedUpperLimits' ],
# 					units 		= [None, None, None, None, None, None, 'pb','pb'],
# 					dataFiles 	= ['orig/newFig16aWZExclusion:Wino-bino(+),Exp.csv','orig/newFig16aWZExclusion:Wino-bino(+),Exp_Down.csv','orig/newFig16aWZExclusion:Wino-bino(+),Exp_Up.csv','orig/newFig16aWZExclusion:Wino-bino(+),Obs.csv','orig/newFig16aWZExclusion:Wino-bino(+),Obs_Down.csv','orig/newFig16aWZExclusion:Wino-bino(+),Obs_Up.csv','orig/AuxFig8aWZExcl.UpperLimitObs.Wino-bino(+)(Deltam).csv', 'orig/AuxFig8bWZExcl.UpperLimitExp.Wino-bino(+)(Deltam).csv'],
# 					dataFormats	= ['csv','csv','csv','csv','csv','csv','csv','csv'])

#+++++++ next mass plane block ++++++++++++++
TChiWZ_plus_off                  = TChiWZoff.addMassPlane(2*[[x, x-y]])
TChiWZ_plus_off.figure           = 'Aux. Fig. 6a'
TChiWZ_plus_off.figureUrl        = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/figaux_06a.png'
TChiWZ_plus_off.dataUrl          = 'https://www.hepdata.net/record/ins1866951?table=AuxFig%208c%20WZ%20Excl.%20Upper%20Limit%20Obs.%20Wino-bino(%2b)%20($\Delta%20m$)'
TChiWZ_plus_off.exclusionDataUrl = 'https://www.hepdata.net/record/ins1866951?version=1&table=Fig%2016b%20WZ%20Exclusion%3A%20Wino-bino(%2B)%20(%24%5CDelta%20m%24)%2C%20Obs'
TChiWZ_plus_off.setSources(dataLabels 	= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits', 'expectedUpperLimits' ],
					units 		= [None, None, None, None, None, None, 'pb','pb'],
					dataFiles 	= ['orig/newFig16aWZExclusion:Wino-bino(+),Exp.csv','orig/newFig16aWZExclusion:Wino-bino(+),Exp_Down.csv','orig/newFig16aWZExclusion:Wino-bino(+),Exp_Up.csv','orig/newFig16aWZExclusion:Wino-bino(+),Obs.csv','orig/newFig16aWZExclusion:Wino-bino(+),Obs_Down.csv','orig/newFig16aWZExclusion:Wino-bino(+),Obs_Up.csv','orig/AuxFig8cWZExcl.UpperLimitObs.Wino-bino(+)(Deltam).csv', 'orig/AuxFig8dWZExcl.UpperLimitExp.Wino-bino(+)(Deltam).csv'],
					dataFormats	= ['csv','csv','csv','csv','csv','csv','csv','csv'])

#+++++++ next mass plane block ++++++++++++++
TChiWZ_higgsino                  = TChiWZoff.addMassPlane([[x, x-y], [x-y/2,x-y]])
TChiWZ_higgsino.validationTarball = "TChiWZoffCd2.tar.gz"
TChiWZ_higgsino.figure           = 'Aux. Fig. 6c'
TChiWZ_higgsino.figureUrl        = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/figaux_06c.png'
TChiWZ_higgsino.dataUrl          = 'https://www.hepdata.net/record/ins1866951?table=AuxFig%208g%20WZ%20Excl.%20Upper%20Limit%20Obs.%20Higgsino%20($\Delta%20m$)'
TChiWZ_higgsino.exclusionDataUrl = 'https://www.hepdata.net/record/ins1866951?table=Fig%2016d%20WZ%20Exclusion:%20Higgsino%20($\Delta%20m$),%20Obs'
TChiWZ_higgsino.setSources(dataLabels 	= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits', 'expectedUpperLimits' ],
					units 		= [None, None, None, None, None, None, 'pb','pb'],
					dataFiles 	= ['orig/Fig16dWZExclusion:Higgsino(Deltam),Exp.csv','orig/Fig16dWZExclusion:Higgsino(Deltam),Exp_Down.csv','orig/Fig16dWZExclusion:Higgsino(Deltam),Exp_Up.csv','orig/Fig16dWZExclusion:Higgsino(Deltam),Obs.csv','orig/Fig16dWZExclusion:Higgsino(Deltam),Obs_Down.csv','orig/Fig16dWZExclusion:Higgsino(Deltam),Obs_Up.csv','orig/AuxFig8gWZExcl.UpperLimitObs.Higgsino(Deltam).csv', 'orig/AuxFig8hWZExcl.UpperLimitExp.Higgsino(Deltam).csv'],
					dataFormats	= ['csv','csv','csv','csv','csv','csv','csv','csv'])

#+++++++ next mass plane block ++++++++++++++
TChiWH                  = TChiWH.addMassPlane(2*[[x, y]])
TChiWH.figure           = 'Fig. 17'
TChiWH.figureUrl        = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/fig_17.png'
TChiWH.dataUrl          = 'https://www.hepdata.net/record/ins1866951?table=AuxFig%209a%20Wh%20Excl.%20Upper%20Limit%20Obs.'
TChiWH.exclusionDataUrl = 'https://www.hepdata.net/record/ins1866951?table=Fig%2017%20Wh%20Exclusion,%20Obs'
TChiWH.setSources(dataLabels 	= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits', 'expectedUpperLimits' ],
					units 		= [None, None, None, None, None, None, 'pb','pb'],
					dataFiles 	= ['orig/Fig17WhExclusion,Exp.csv','orig/Fig17WhExclusion,Exp_Down.csv','orig/Fig17WhExclusion,Exp_Up.csv','orig/Fig17WhExclusion,Obs.csv','orig/Fig17WhExclusion,Obs_Down.csv','orig/Fig17WhExclusion,Obs_Up.csv','orig/AuxFig9aWhExcl.UpperLimitObs..csv', 'orig/AuxFig9bWhExcl.UpperLimitExp..csv'],
					dataFormats	= ['csv','csv','csv','csv','csv','csv','csv','csv'])


databaseCreator.create()
