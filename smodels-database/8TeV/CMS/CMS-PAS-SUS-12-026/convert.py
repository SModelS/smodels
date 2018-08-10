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
info = MetaInfoInput('CMS-PAS-SUS-12-026')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS12026'
info.sqrts = 8
info.lumi = 9.2
info.prettyName = '>= 3 leptons (+jets) + Etmiss'
info.private = False
info.arxiv = ''
info.contact =''
info.publication = ''
info.comment = 'Only CDS entry http://cdsweb.cern.ch/record/1496879/files/SUS-12-026-pas.pdf'
info.supersededBy = 'CMS-SUS-13-002'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked =''
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked =''
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure = 'Fig. 7'
T1tttt_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12026/curve_T1tttt_overlay_observed.pdf'
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS12026'
T1tttt_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12026/Limit_T1tttt.root'
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12026/Limit_T1tttt.root'
T1tttt_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12026/Limit_T1tttt.root'
T1tttt_1.setSources(dataLabels= ['expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/Limit_T1tttt.root', 'orig/T1tttt_Excl.txt', 'orig/Limit_T1tttt.root'],
                 dataFormats= ['root', 'txt', 'root'],objectNames= ['LimitExp', 'None', 'LimitObs'],units= ['fb', None, 'fb'])
T1ttttoff.addMassPlane(T1tttt_1)



databaseCreator.create()
