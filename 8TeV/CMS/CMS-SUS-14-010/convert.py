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
info = MetaInfoInput('CMS-SUS-14-010')
info.url = "https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14010"
info.sqrts = 8
info.lumi = 19.5
info.prettyName = 'b-jets + 4Ws (combination of 5 analyses, 0,1,2,...leptons)'
info.private = False
info.arxiv ="http://arxiv.org/abs/arXiv:1412.4109"
info.contact =""
info.publication ="Phys. Lett. B 745 (2015) 5"
info.comment =""
info.supersedes =""
info.supersededBy =""


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked =""
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked =""
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription ="None"
T1ttttoff.condition ="None"
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure ='5a'
T1tttt_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14010/T1CombXSEC.pdf'
T1tttt_1.dataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14010/T1tttt_combination.root"
T1tttt_1.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14010/T1CombXSEC.pdf'
T1tttt_1.dataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14010/T1tttt_combination.root"
T1tttt_1.exclusionDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14010/T1tttt_combination.root"
T1tttt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/T1tttt_combination.root', 'orig/T1tttt_combination.root', 'orig/T1tttt_combination.root', 'orig/T1tttt_combination.root', 'orig/T1tttt_combination.root', 'orig/T1tttt_combination.root', 'orig/T1tttt_combination.root', 'orig/T1tttt_combination.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['expected_graph', 'expected_minus_graph', 'expected_plus_graph', 'expectedUL', 'observed_graph', 'observed_minus_graph', 'observed_plus_graph', 'observedUL'],units= [None, None, None, 'fb', None, None, None, 'fb'])
T1ttttoff.addMassPlane(T1tttt_1)



databaseCreator.create()
