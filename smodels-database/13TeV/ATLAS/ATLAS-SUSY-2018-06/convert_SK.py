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
info 			= MetaInfoInput('ATLAS-SUSY-2018-06')
info.url 		= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-06/'
info.sqrts 		= 13
info.lumi 		= 139
info.prettyName = '3 leptons EW-ino'
info.private 	= False
info.arxiv 		= 'https://arxiv.org/abs/1912.08479'
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

#+++++++ next mass plane block ++++++++++++++
TChiWZ_1                  = TChiWZ.addMassPlane(2*[[x, y]])
TChiWZ_1.figure           = 'Aux. Fig. 12'
TChiWZ_1.figureUrl        = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-06/figaux_12a.png'
TChiWZ_1.dataUrl          = 'https://doi.org/10.17182/hepdata.91127.v1/t15'
TChiWZ_1.exclusionDataUrl = 'https://doi.org/10.17182/hepdata.91127.v1/t9'
TChiWZ_1.setSources(dataLabels 	= ['expExclusion', 'obsExclusion','upperLimits','expectedUpperLimits' ],
					units 		= [None, None,'pb','pb'],
					dataFiles 	= ['orig/ExpectedLimit3LeRJR.csv','orig/ObservedLimit3LeRJR.csv','orig/ObservedUpperLimits.csv', 'orig/ExpectedUpperLimits.csv'],
					dataFormats	= ['csv','csv','csv','csv'])

TChiWZoff.addMassPlane(TChiWZ_1)

databaseCreator.create()
