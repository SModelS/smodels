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
info = MetaInfoInput('ATLAS-CONF-2013-028')
info.sqrts = '8.0'
info.private = False
info.lumi = '20.7'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-028/'
info.prettyName = 'ATLAS hadronic stau'
info.comment = 'Will be supersededBy ATLAS-SUSY-2013-14 (at the moment cannto be implemented -Jan2017)'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiChipmStauL = dataset.addTxName('TChiChipmStauL')
TChiChipmStauL.checked ="AL"
TChiChipmStauL.constraint ="2.*([[['nu'],['ta']],[['ta+'],['ta-']]] + [[['ta'],['nu']],[['ta+'],['ta-']]]+[[['nu'],['ta']],[['ta-'],['ta+']]] + [[['ta'],['nu']],[['ta-'],['ta+']]])"
TChiChipmStauL.conditionDescription ="[[['nu'],['ta']],[['ta+'],['ta-']]] ~ [[['ta'],['nu']],[['ta+'],['ta-']]],[[['nu'],['ta']],[['ta+'],['ta-']]] ~ [[['nu'],['ta']],[['ta-'],['ta+']]], [[['nu'],['ta']],[['ta+'],['ta-']]] ~ [[['ta'],['nu']],[['ta-'],['ta+']]]"
TChiChipmStauL.condition ="Csim([[['nu'],['ta']],[['ta+'],['ta-']]],[[['ta'],['nu']],[['ta+'],['ta-']]],[[['nu'],['ta']],[['ta-'],['ta+']]],[[['ta'],['nu']],[['ta-'],['ta+']]])"
TChiChipmStauL.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChiChipmStauL050 = TChiChipmStauL.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChiChipmStauL050.figure = 'Fig.(aux) 8a'
TChiChipmStauL050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-028/figaux_08a.png'
TChiChipmStauL050.dataUrl = 'Not defined'
TChiChipmStauL050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/TChiChipmStauL_excl.txt', 'orig/TChiChipmStauL.txt'],
                 dataFormats= ['txt', 'txt'])

#+++++++ next txName block ++++++++++++++
TChipChimStauSnu = dataset.addTxName('TChipChimStauSnu')
TChipChimStauSnu.checked ="A"
TChipChimStauSnu.constraint ="[[['ta-'],['nu']],[['nu'],['ta+']]] + [[['ta+'],['nu']],[['nu'],['ta-']]] + [[['ta+'],['nu']],[['ta-'],['nu']]] + [[['nu'],['ta+']],[['nu'],['ta-']]]"
TChipChimStauSnu.conditionDescription ="[[['ta-'],['nu']],[['nu'],['ta+']]] ~ [[['ta+'],['nu']],[['nu'],['ta-']]], [[['ta-'],['nu']],[['nu'],['ta+']]] ~ [[['ta+'],['nu']],[['ta-'],['nu']]] ,[[['ta-'],['nu']],[['nu'],['ta+']]] ~ [[['nu'],['ta+']],[['nu'],['ta-']]]"
TChipChimStauSnu.condition ="Csim([[['ta-'],['nu']],[['nu'],['ta+']]],[[['ta+'],['nu']],[['nu'],['ta-']]],[[['ta+'],['nu']],[['ta-'],['nu']]],[[['nu'],['ta+']],[['nu'],['ta-']]])"
TChipChimStauSnu.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChipChimStauSnu050 = TChipChimStauSnu.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChipChimStauSnu050.figure = 'Fig.(aux) 9a'
TChipChimStauSnu050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-028/figaux_09a.png'
TChipChimStauSnu050.dataUrl = 'Not defined'
TChipChimStauSnu050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/TChipChimStauSnu_excl.txt', 'orig/TChipChimStauSnu.txt'],
                 dataFormats= ['txt', 'txt'])



databaseCreator.create()
