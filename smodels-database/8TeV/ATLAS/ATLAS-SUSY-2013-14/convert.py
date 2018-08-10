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
info = MetaInfoInput('ATLAS-SUSY-2013-14')
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.publication = 'link.springer.com/article/10.1007/JHEP10(2014)096'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-14/'
info.arxiv = 'http://arxiv.org/abs/1407.0350'
info.contact = '?'
info.prettyName = 'ATLAS di-tau'
info.supersedes = 'ATLAS-CONF-2013-028'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiChipmStauL = dataset.addTxName('TChiChipmStauL')
TChiChipmStauL.constraint ="2.*([[['ta'],['ta']],[['nu'],['ta']]]+[[['ta'],['ta']],[['ta'],['nu']]])"
TChiChipmStauL.conditionDescription ="[[['ta'],['ta']],[['nu'],['ta']]] ~ [[['ta'],['ta']],[['ta'],['nu']]]"
TChiChipmStauL.condition ="Csim([[['ta'],['ta']],[['nu'],['ta']]],[[['ta'],['ta']],[['ta'],['nu']]])"
TChiChipmStauL.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChiChipmStauL050 = TChiChipmStauL.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChiChipmStauL050.figure = 'Fig.(aux) 11b'
TChiChipmStauL050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-14/figaux_11b.png'
TChiChipmStauL050.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1304288/d29'
TChiChipmStauL050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiChipmStauL.txt', 'orig/limit_TChiChipmStauL.txt'],
                 dataFormats= ['txt', 'txt'])

#+++++++ next txName block ++++++++++++++
TChipChimStauSnu = dataset.addTxName('TChipChimStauSnu')
TChipChimStauSnu.constraint ="[[['ta-'],['nu']],[['nu'],['ta+']]] + [[['ta+'],['nu']],[['nu'],['ta-']]] + [[['ta+'],['nu']],[['ta-'],['nu']]] + [[['nu'],['ta+']],[['nu'],['ta-']]]"
TChipChimStauSnu.conditionDescription ="[[['ta-'],['nu']],[['nu'],['ta+']]] ~ [[['ta+'],['nu']],[['nu'],['ta-']]], [[['ta-'],['nu']],[['nu'],['ta+']]] ~ [[['ta+'],['nu']],[['ta-'],['nu']]], [[['ta-'],['nu']],[['nu'],['ta+']]] ~ [[['nu'],['ta+']],[['nu'],['ta-']]]"
TChipChimStauSnu.condition ="Csim([[['ta-'],['nu']],[['nu'],['ta+']]],[[['ta+'],['nu']],[['nu'],['ta-']]],[[['ta+'],['nu']],[['ta-'],['nu']]],[[['nu'],['ta+']],[['nu'],['ta-']]])"
TChipChimStauSnu.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChipChimStauSnu050 = TChipChimStauSnu.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChipChimStauSnu050.figure = 'Fig.(aux) 11a'
TChipChimStauSnu050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-14/figaux_11a.png'
TChipChimStauSnu050.dataUrl = 'Not defined'
TChipChimStauSnu050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChipChimStauSnu.txt', 'orig/limit_TChipChimStauSnu.txt'],
                 dataFormats= ['txt', 'txt'])



databaseCreator.create()
