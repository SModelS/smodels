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
info = MetaInfoInput('CMS-PAS-SUS-13-023')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13023'
info.sqrts = 8
info.lumi = 18.9 
info.prettyName = 'hadronic stop'
info.private = False
info.comment = 'PAS:http://inspirehep.net/record/1387812/files/SUS-13-023-pas.pdf'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked =''
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription = None
T2tt.condition = None
T2tt.source = "CMS"
T2tt.massConstraint = None
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked =''
T2ttoff.constraint ="[[['b','W']],[['b','W']]]"
T2ttoff.conditionDescription = None
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0'], ['dm <= 169.0']]
T2ttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure = 'Figure 13'
T2tt_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2tt__observed_xsection_UL.png'
T2tt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13023'
T2tt_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2tt__observed_xsection_UL.root'
T2tt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2tt__observed_xsection_UL.root'
T2tt_1.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2tt__observed_xsection_UL.root'
T2tt_1.setSources(dataLabels= ['obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/T2tt_Obs.txt', 'orig/T2tt_ObsMinus.txt', 'orig/T2tt_ObsPlus.txt', 'orig/T2tt__observed_xsection_UL.root'],
                 dataFormats= ['txt', 'txt', 'txt', 'canvas'],objectNames= ['None', 'None', 'None', 'T2tt__observed_xsection_UL'],indices= [None, None, None, 2])
T2ttoff.addMassPlane(T2tt_1)

#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.checked = ''
T6bbWW.constraint ="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription = None
T6bbWW.condition = None
T6bbWW.source = "CMS"
T6bbWW.massConstraint = None
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.checked =''
T6bbWWoff.constraint ="2.3*([[['b'],['jet','jet']],[['b'],['jet','jet']]])"
T6bbWWoff.conditionDescription = None
T6bbWWoff.condition = None
T6bbWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6bbWWoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T6bbWW_1 = T6bbWW.addMassPlane(2*[[x, (0.75*x+0.25*y), y]])
T6bbWW_1.figure = 'Figure 13'
T6bbWW_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p25_observed_xsection_UL.png'
T6bbWW_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p25_observed_xsection_UL.root'
T6bbWW_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p25_observed_xsection_UL.root'
T6bbWW_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p25_observed_xsection_UL.root'
T6bbWW_1.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p25_observed_xsection_UL.root'
T6bbWW_1.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2bbWW_x025_Obs.txt', 'orig/T2bw_0p25_observed_xsection_UL.root'],
                 dataFormats= ['txt', 'canvas'],objectNames= ['None', 'T2bw_0p25_observed_xsection_UL'],indices= [None, 2])
T6bbWWoff.addMassPlane(T6bbWW_1)
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane(2*[[x, (0.25*x+0.75*y), y]])
T6bbWW_2.figure = 'Figure 13'
T6bbWW_2.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023//T2bw_0p75_observed_xsection_UL.pdf'
T6bbWW_2.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p75_observed_xsection_UL.root'
T6bbWW_2.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p75_observed_xsection_UL.root'
T6bbWW_2.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p75_observed_xsection_UL.root'
T6bbWW_2.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p75_observed_xsection_UL.root'
T6bbWW_2.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2bbWW_x075_Obs.txt', 'orig/T2bw_0p75_observed_xsection_UL.root'],
                 dataFormats= ['txt', 'canvas'],objectNames= ['None', 'T2bw_0p75_observed_xsection_UL'],indices= [None, 2])
T6bbWWoff.addMassPlane(T6bbWW_2)
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane(2*[[x, (0.5*x+0.5*y), y]])
T6bbWW_3.figure = 'Figure 13'
T6bbWW_3.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023//T2bw_0p50_observed_xsection_UL.pdf'
T6bbWW_3.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p50_observed_xsection_UL.root'
T6bbWW_3.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p50_observed_xsection_UL.root'
T6bbWW_3.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p50_observed_xsection_UL.root'
T6bbWW_3.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13023/T2bw_0p50_observed_xsection_UL.root'
T6bbWW_3.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2bbWW_x050_Obs.txt', 'orig/T2bw_0p50_observed_xsection_UL.root'],
                 dataFormats= ['txt', 'canvas'],objectNames= ['None', 'T2bw_0p50_observed_xsection_UL'],indices= [None, 2])
T6bbWWoff.addMassPlane(T6bbWW_3)



databaseCreator.create()
