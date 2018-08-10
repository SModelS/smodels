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
info = MetaInfoInput('CMS-SUS-13-019')
info.url ='https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13019'
info.sqrts = 8
info.lumi = 19.5
info.prettyName = '>= 2 jets + Etmiss, M_T2'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1502.04358'
info.contact =''
info.publication ='http://link.springer.com/article/10.1007%2FJHEP05%282015%29078'
info.comment =''
info.supersedes =''
info.supersededBy =''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked =''
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = 'None'
T1tttt.condition = 'None'
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked =''
T1ttttoff.constraint = "[[['b','W','b','W']],[['b','W','b','W']]]" 
T1ttttoff.conditionDescription ='None'
T1ttttoff.condition ='None'
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure = 'Fig. 13c'
T1tttt_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1tttt-SUS13019-final_XSEC.png'       
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13019'
T1tttt_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2tt-SUS13019-final_XSEC.root'
T1tttt_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2tt-SUS13019-final_XSEC.root'
T1tttt_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2tt-SUS13019-final_XSEC.root'
T1tttt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/Extracted_T1tttt-SUS13019-final_XSEC.root', 'orig/Extracted_T1tttt-SUS13019-final_XSEC.root', 'orig/Extracted_T1tttt-SUS13019-final_XSEC.root', 'orig/Extracted_T1tttt-SUS13019-final_XSEC.root', 'orig/Extracted_T1tttt-SUS13019-final_XSEC.root', 'orig/Extracted_T1tttt-SUS13019-final_XSEC.root', 'orig/Extracted_T1tttt-SUS13019-final_XSEC.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['gr_cTH_exp_combined', 'gr_cTH_psig_combined', 'gr_cTH_msig_combined', 'gr_cTH_obs_combined', 'gr_mTH_obs_combined', 'gr_pTH_obs_combined', 'XSec_limit_combined'])
T1ttttoff.addMassPlane(T1tttt_1)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked =''
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.conditionDescription = 'None'
T2tt.condition ='None'
T2tt.source = "CMS"
T2tt.massConstraint = None
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint ="[[['b','W']],[['b','W']]]"
T2ttoff.conditionDescription ='None'
T2ttoff.condition ='None'
T2ttoff.massConstraint = [['dm <= 169.0'], ['dm <= 169.0']]
T2ttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure = 'Fig. 12c'
T2tt_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2tt-SUS13019-final_XSEC.png'
T2tt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13019'
T2tt_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2tt-SUS13019-final_XSEC.root'
T2tt_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2tt-SUS13019-final_XSEC.root'
T2tt_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2tt-SUS13019-final_XSEC.root'
T2tt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/Extracted_T2tt-SUS13019-final_XSEC.root', 'orig/Extracted_T2tt-SUS13019-final_XSEC.root', 'orig/Extracted_T2tt-SUS13019-final_XSEC.root', 'orig/Extracted_T2tt-SUS13019-final_XSEC.root', 'orig/Extracted_T2tt-SUS13019-final_XSEC.root', 'orig/Extracted_T2tt-SUS13019-final_XSEC.root', 'orig/Extracted_T2tt-SUS13019-final_XSEC.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['gr_cTH_exp_combined', 'gr_cTH_psig_combined', 'gr_cTH_msig_combined', 'gr_cTH_obs_combined', 'gr_mTH_obs_combined', 'gr_pTH_obs_combined', 'XSec_limit_combined'])
T2ttoff.addMassPlane(T2tt_1)

#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked = ''
T1bbbb.constraint = "[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription = 'None'
T1bbbb.condition = 'None'
T1bbbb.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb_1 = T1bbbb.addMassPlane(2*[[x, y]])
T1bbbb_1.figure ='Fig. 13b'
T1bbbb_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1bbbb-SUS13019-final_XSEC.png'
T1bbbb_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1bbbb-SUS13019-final_XSEC.png'
T1bbbb_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13019'
T1bbbb_1.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1bbbb-SUS13019-final_XSEC.root'
T1bbbb_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1bbbb-SUS13019-final_XSEC.root'
T1bbbb_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1bbbb-SUS13019-final_XSEC.root'
T1bbbb_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/Extracted_T1bbbb-SUS13019-final_XSEC.root', 'orig/Extracted_T1bbbb-SUS13019-final_XSEC.root', 'orig/Extracted_T1bbbb-SUS13019-final_XSEC.root', 'orig/Extracted_T1bbbb-SUS13019-final_XSEC.root', 'orig/Extracted_T1bbbb-SUS13019-final_XSEC.root', 'orig/Extracted_T1bbbb-SUS13019-final_XSEC.root', 'orig/Extracted_T1bbbb-SUS13019-final_XSEC.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['gr_cTH_exp_combined', 'gr_cTH_msig_combined', 'gr_cTH_psig_combined', 'gr_cTH_obs_combined', 'gr_mTH_obs_combined', 'gr_pTH_obs_combined', 'XSec_limit_combined'])

#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.checked =''
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription = "None"
T2.condition = "None"
T2.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2_1 = T2.addMassPlane(2*[[x, y]])
T2_1.figure ='Fig. 13a'
T2_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2qq-SUS13019-final_XSEC.png'
T2_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13019'
T2_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1qqqq-SUS13019-final_XSEC.root'
T2_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1qqqq-SUS13019-final_XSEC.root'
T2_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1qqqq-SUS13019-final_XSEC.root'
T2_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/Extracted_T2qq-SUS13019-final_XSEC.root', 'orig/Extracted_T2qq-SUS13019-final_XSEC.root', 'orig/Extracted_T2qq-SUS13019-final_XSEC.root', 'orig/Extracted_T2qq-SUS13019-final_XSEC.root', 'orig/Extracted_T2qq-SUS13019-final_XSEC.root', 'orig/Extracted_T2qq-SUS13019-final_XSEC.root', 'orig/Extracted_T2qq-SUS13019-final_XSEC.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['gr_cTH_exp_combined', 'gr_cTH_msig_combined', 'gr_cTH_psig_combined', 'gr_cTH_obs_combined', 'gr_mTH_obs_combined', 'gr_pTH_obs_combined', 'XSec_limit_combined'])

#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1_1 = T1.addMassPlane(2*[[x, y]])
T1_1.figure ='Fig. 13a'
T1_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1qqqq-SUS13019-final_XSEC.png'
T1_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13019'
T1_1.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1qqqq-SUS13019-final_XSEC.root'
T1_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1qqqq-SUS13019-final_XSEC.root'
T1_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T1qqqq-SUS13019-final_XSEC.root'
T1_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/Extracted_T1qqqq-SUS13019-final_XSEC.root', 'orig/Extracted_T1qqqq-SUS13019-final_XSEC.root', 'orig/Extracted_T1qqqq-SUS13019-final_XSEC.root', 'orig/Extracted_T1qqqq-SUS13019-final_XSEC.root', 'orig/Extracted_T1qqqq-SUS13019-final_XSEC.root', 'orig/Extracted_T1qqqq-SUS13019-final_XSEC.root', 'orig/Extracted_T1qqqq-SUS13019-final_XSEC.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['gr_cTH_exp_combined', 'gr_cTH_msig_combined', 'gr_cTH_psig_combined', 'gr_cTH_obs_combined', 'gr_mTH_obs_combined', 'gr_pTH_obs_combined', 'XSec_limit_combined'])

#+++++++ next txName block ++++++++++++++
T2bb = dataset.addTxName('T2bb')
T2bb.checked =''
T2bb.constraint = "[[['b']],[['b']]]"
T2bb.conditionDescription ='None'
T2bb.condition ='None'
T2bb.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2bb_1 = T2bb.addMassPlane(2*[[x, y]])
T2bb_1.figure ='Fig. 12b'
T2bb_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2bb-SUS13019-final_XSEC.png'
T2bb_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13019'
T2bb_1.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2bb-SUS13019-final_XSEC.root'
T2bb_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2bb-SUS13019-final_XSEC.root'
T2bb_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13019/T2bb-SUS13019-final_XSEC.root'
T2bb_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/Extracted_T2bb-SUS13019-final_XSEC.root', 'orig/Extracted_T2bb-SUS13019-final_XSEC.root', 'orig/Extracted_T2bb-SUS13019-final_XSEC.root', 'orig/Extracted_T2bb-SUS13019-final_XSEC.root', 'orig/Extracted_T2bb-SUS13019-final_XSEC.root', 'orig/Extracted_T2bb-SUS13019-final_XSEC.root', 'orig/Extracted_T2bb-SUS13019-final_XSEC.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['gr_cTH_exp_combined', 'gr_cTH_msig_combined', 'gr_cTH_psig_combined', 'gr_cTH_obs_combined', 'gr_mTH_obs_combined', 'gr_pTH_obs_combined', 'XSec_limit_combined'])



databaseCreator.create()
