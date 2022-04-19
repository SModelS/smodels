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
info = MetaInfoInput('CMS-SUS-12-024')
info.sqrts = '8.0'
info.private = False
info.lumi = '19.4'
info.publication = 'http://www.sciencedirect.com/science/article/pii/S0370269313005339'
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS12024'
info.arxiv = 'http://arxiv.org/abs/1305.2390'
info.contact = 'Keith Ulmer <keith.ulmer@cern.ch>, Josh Thompson <joshua.thompson@cern.ch>, Alessandro Gaz <alessandro.gaz@cern.ch>'
info.prettyName = '0 leptons + >= 3 (1b-)jets + Etmiss'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="A"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane(2*[[x, y]])
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.C"
T1bbbb.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/graph_smoothed_Obs_T1bbbb.root', 'orig/hXsec_exp_corr_T1bbbb.root'],
                 dataFormats= ['root', 'root'],objectNames= ['graph_smoothed_Obs_T1bbbb', 'hXsec_exp_corr'])

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="A"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="A"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.conditionDescription ="None"
T1ttttoff.condition ="None"
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane(2*[[x, y]])
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusion_corrected.C"
T1tttt.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/graph_smoothed_Obs_T1tttt.root', 'orig/hXsec_obs_final_T1tttt.root'],
                 dataFormats= ['root', 'root'],objectNames= ['graph_smoothed_Obs_T1tttt', 'hXsec_obs_final'])
T1ttttoff.addMassPlane(T1tttt)



databaseCreator.create()
