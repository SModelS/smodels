#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse
import types

argparser = argparse.ArgumentParser(description =  
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath', 
help = 'path to the package smodels_utils',\
type = str)
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
help = 'path to the package smodels_utils',\
type = str)
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
info = MetaInfoInput('CMS-PAS-SUS-13-016')
info.sqrts = '8.0'
info.private = False
info.lumi = '19.7'
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13016'
info.prettyName = '2 OS leptons + >= 4 (2b-)jets + Etmiss'




#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("sr0")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "sr0", observedN = 1, expectedBG = 1.2 , bgError = 1.04862, upperLimit = '2.110E-01*fb', expectedUpperLimit = '2.119E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="VM"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.condition ="None"
T1tttt.conditionDescription ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.addSource('efficiencyMap',"orig/Acceptance_OS_T1tttt.root", "root", objectName = "AccANN")
T1tttt_1.addSource('obsExclusion',"orig/Results_SUS13016.root", "root", objectName = "GR_ns")
T1tttt_1.addSource('obsExclusionM1',"orig/Results_SUS13016.root", "root", objectName = "GR_ns0")
T1tttt_1.addSource('obsExclusionP1',"orig/Results_SUS13016.root", "root", objectName = "GR_ns1")
T1tttt_1.addSource('expExclusion',"orig/Results_SUS13016.root", "root", objectName = "GR_nsE")
T1tttt_1.addSource('expExclusionM1',"orig/Results_SUS13016.root", "root", objectName = "GR_nsE0")
T1tttt_1.addSource('expExclusionP1',"orig/Results_SUS13016.root", "root", objectName = "GR_nsE1")
T1tttt_1.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13016/Results_SUS13016.root"
T1ttttoff.addMassPlane(T1tttt_1)

databaseCreator.create()
