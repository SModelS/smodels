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
info = MetaInfoInput('CMS-PAS-SUS-16-024')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/'
info.sqrts = 13
info.lumi = 12.9
info.private = False 
info.comment = 'https://cds.cern.ch/record/2205168. Superseded by CMS-SUS-16-039.'
info.supersededBy = 'CMS-SUS-16-039'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiChipmStauL = dataset.addTxName('TChiChipmStauL')
TChiChipmStauL.checked = '' 
TChiChipmStauL.constraint = "[[['ta'],['ta']],[['nu'],['L']]]"
TChiChipmStauL.conditionDescription = None
TChiChipmStauL.condition = None
TChiChipmStauL.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChiChipmStauL_1 = TChiChipmStauL.addMassPlane(2*[[x, .5*(x+y), y]])
TChiChipmStauL_1.figure = 'Figure 9'
TChiChipmStauL_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_009.png'
TChiChipmStauL_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_009.root'
TChiChipmStauL_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_009.root'
TChiChipmStauL_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_009.root'
TChiChipmStauL_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_009.root'
TChiChipmStauL_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/CMS-PAS-SUS-16-024_Figure_009.root', 'orig/CMS-PAS-SUS-16-024_Figure_009.root', 'orig/CMS-PAS-SUS-16-024_Figure_009.root', 'orig/CMS-PAS-SUS-16-024_Figure_009.root', 'orig/CMS-PAS-SUS-16-024_Figure_009.root', 'orig/CMS-PAS-SUS-16-024_Figure_009.root', 'orig/CMS-PAS-SUS-16-024_Figure_009.root', 'orig/CMS-PAS-SUS-16-024_Figure_009.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['gr_exp_smoothed', 'gr_em1s_smoothed', 'gr_ep1s_smoothed', 'exp_xs', 'gr_obs_smoothed', 'gr_om1s_smoothed', 'gr_op1s_smoothed', 'obs_xs'])

#+++++++ next txName block ++++++++++++++
TChiWH = dataset.addTxName('TChiWH')
TChiWH.checked = ''
TChiWH.constraint = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription = None
TChiWH.condition = None
TChiWH.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChiWH_1 = TChiWH.addMassPlane(2*[[x, y]])
TChiWH_1.figure = 'Figure 10-b'
TChiWH_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_010-b.png'
TChiWH_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_010-b.root'
TChiWH_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_010-b.root'
TChiWH_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_010-b.root'
TChiWH_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_010-b.root'
TChiWH_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/CMS-PAS-SUS-16-024_Figure_010-b.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-b.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-b.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-b.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-b.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-b.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-b.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-b.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['gr_exp', 'gr_em1s', 'gr_ep1s', 'exp_xs', 'gr_obs', 'gr_om1s', 'gr_op1s', 'obs_xs'])

#+++++++ next txName block ++++++++++++++
TChiWZ = dataset.addTxName('TChiWZ')
TChiWZ.checked = ''
TChiWZ.constraint = "[[['W']],[['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition = None
TChiWZ.source = "CMS"
TChiWZ.massConstraint = None
TChiWZoff = dataset.addTxName('TChiWZoff')
TChiWZoff.checked = ''
TChiWZoff.constraint = "[[['e+','e-']],[['l','nu']]]+[[['mu+','mu-']],[['l','nu']]]"
TChiWZoff.conditionDescription = None
TChiWZoff.condition = None
TChiWZoff.massConstraint = [['dm <= 76.0'], ['dm <= 86.0']]
TChiWZoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChiWZ_1 = TChiWZ.addMassPlane(2*[[x, y]])
TChiWZ_1.figure = 'Figure 10-a'
TChiWZ_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_010-a.png'
TChiWZ_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_010-a.root'
TChiWZ_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_010-a.root'
TChiWZ_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_010-a.root'
TChiWZ_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-024/CMS-PAS-SUS-16-024_Figure_010-a.root'
TChiWZ_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/CMS-PAS-SUS-16-024_Figure_010-a.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-a.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-a.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-a.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-a.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-a.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-a.root', 'orig/CMS-PAS-SUS-16-024_Figure_010-a.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['gr_exp_smoothed', 'gr_em1s_smoothed', 'gr_ep1s_smoothed', 'exp_xs', 'gr_obs_smoothed', 'gr_om1s_smoothed', 'gr_op1s_smoothed', 'obs_xs'])
TChiWZoff.addMassPlane(TChiWZ_1)



databaseCreator.create()
