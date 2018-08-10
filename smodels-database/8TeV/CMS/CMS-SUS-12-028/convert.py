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
info = MetaInfoInput('CMS-SUS-12-028')
info.url ="https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS12028"
info.sqrts = 8
info.lumi = 11.7
info.prettyName = 'jets + Etmiss, alpha_T'
info.private = False
info.arxiv ='http://arxiv.org/abs/1303.2985'
info.contact ='Edward Laird <edward.laird@cern.ch>, Rob Bainbridge <robert.bainbridge@cern.ch>'
info.publication ='Eur. Phys. J. C 73 (2013) 2568'
info.comment ="Only the expected upper limit exclusion line is present - no observed one"
info.supersedes ='CMS-PAS-SUS-12-016'
info.supersededBy =""


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked =""
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = "None"
T1tttt.condition ="None"
T1tttt.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure ="Figure 8"
T1tttt_1.figureUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/t1tttt.pdf"
T1tttt_1.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T1tttt.root"
T1tttt_1.histoDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T1tttt.pdf"
T1tttt_1.dataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T1tttt.root"
T1tttt_1.exclusionDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T1tttt.root"
T1tttt_1.setSources(dataLabels= ['expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T1tttt.root', 'orig/T1tttt.root', 'orig/T1tttt.root'],
                 dataFormats= ['root', 'root', 'root'],objectNames= ['ExpectedUpperLimit', 'UpperLimit_graph', 'UpperLimit'])

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked =""
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure ="t2tt"
T2tt_1.figureUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/t2tt.pdf"
T2tt_1.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T2tt.root"
T2tt_1.histoDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T2tt.pdf"
T2tt_1.dataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T2tt.root"
T2tt_1.exclusionDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T2tt.root"
T2tt_1.setSources(dataLabels= ['expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2tt.root', 'orig/T2tt_exc.dat', 'orig/T2tt.root'],
                 dataFormats= ['root', 'txt', 'root'],objectNames= ['ExpectedUpperLimit', 'None', 'UpperLimit'])

#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked =""
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb_1 = T1bbbb.addMassPlane(2*[[x, y]])
T1bbbb_1.figure ="Figure 8"
T1bbbb_1.figureUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/t1bbbb.pdf"
T1bbbb_1.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T1bbbb.root"
T1bbbb_1.histoDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/t1bbbb.pdf"
T1bbbb_1.dataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T1bbbb.root"
T1bbbb_1.exclusionDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T1bbbb.root"
T1bbbb_1.setSources(dataLabels= ['expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T1bbbb.root', 'orig/T1bbbb.root', 'orig/T1bbbb.root'],
                 dataFormats= ['root', 'root', 'root'],objectNames= ['ExpectedUpperLimit', 'UpperLimit_graph', 'UpperLimit'])

#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.checked =""
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2_1 = T2.addMassPlane(2*[[x, y]])
T2_1.figure ="Figure 8"
T2_1.figureUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/t2.pdf"
T2_1.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T2.root"
T2_1.histoDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/t2.pdf"
T2_1.dataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T2.root"
T2_1.exclusionDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T2.root"
T2_1.setSources(dataLabels= ['expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2.root', 'orig/T2.root', 'orig/T2.root'],
                 dataFormats= ['root', 'root', 'root'],objectNames= ['T2_ExpectedUpperLimit', 'UpperLimit_xs0p8', 'T2_UpperLimit'])

#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.checked =""
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1_1 = T1.addMassPlane(2*[[x, y]])
T1_1.figure ="Figure 8"
T1_1.figureUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/t1.pdf"
T1_1.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T1.root"
T1_1.histoDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/t1.pdf"
T1_1.dataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T1.root"
T1_1.exclusionDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T1.root"
T1_1.setSources(dataLabels= ['expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T1.root', 'orig/T1.root', 'orig/T1.root'],
                 dataFormats= ['root', 'root', 'root'],objectNames= ['ExpectedUpperLimit', 'UpperLimit_graph', 'UpperLimit'])

#+++++++ next txName block ++++++++++++++
T2bb = dataset.addTxName('T2bb')
T2bb.checked =""
T2bb.constraint ="[[['b']],[['b']]]"
T2bb.conditionDescription ="None"
T2bb.condition ="None"
T2bb.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2bb_1 = T2bb.addMassPlane(2*[[x, y]])
T2bb_1.figure ="Figure 8"
T2bb_1.figureUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/t2bb.pdf"
T2bb_1.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T2bb.root"
T2bb_1.histoDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/t2bb.pdf"
T2bb_1.dataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T2bb.root"
T2bb_1.exclusionDataUrl ="https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12028/T2bb.root"
T2bb_1.setSources(dataLabels= ['expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2bb.root', 'orig/T2bb.root', 'orig/T2bb.root'],
                 dataFormats= ['root', 'root', 'root'],objectNames= ['ExpectedUpperLimit', 'UpperLimit_graph', 'UpperLimit'])



databaseCreator.create()
