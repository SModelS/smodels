#!/usr/bin/env python3

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
help = 'path to smodels',\
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
info = MetaInfoInput('CMS-PAS-SUS-15-002')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
info.sqrts = 13
info.lumi = 2.2
info.prettyName = '>= 4jets + Etmiss, HT, HTmiss'
info.private = False
info.arxiv = ''
info.contact = ''
info.publication = ''
info.comment = 'Only CDS entry https://cds.cern.ch/record/2114817'
info.supersededBy = 'CMS-SUS-15-002'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked = ''
T1bbbb.constraint = "[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription = None
T1bbbb.condition = None
T1bbbb.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb_1 = T1bbbb.addMassPlane(2*[[x, y]])
T1bbbb_1.figure = 'Fig.8-a'
T1bbbb_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/CMS-PAS-SUS-15-002_Figure_008-a.png'
T1bbbb_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
T1bbbb_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
T1bbbb_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1bbbbSmooth.root'
T1bbbb_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
T1bbbb_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/MassScanT1bbbbSmooth.root', 'orig/MassScanT1bbbbSmooth.root', 'orig/MassScanT1bbbbSmooth.root', 'orig/MassScanT1bbbbSmooth.root', 'orig/MassScanT1bbbbSmooth.root', 'orig/MassScanT1bbbbSmooth.root', 'orig/MassScanT1bbbbSmooth.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['ExpLim', 'ExpLimSdn', 'ExpLimSup', 'ObsLim', 'ObsLimSdn', 'ObsLimSup', 'MassScan2D'])

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked = ''
T1tttt.constraint =  "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked = ''
T1ttttoff.constraint = "[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure = 'Fig.8-b'
T1tttt_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/CMS-PAS-SUS-15-002_Figure_008-b.png'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
T1tttt_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1ttttSmooth.root'
T1tttt_1.exclusionDataUrl =  'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1ttttSmooth.root'
T1tttt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/MassScanT1ttttSmooth.root', 'orig/MassScanT1ttttSmooth.root', 'orig/MassScanT1ttttSmooth.root', 'orig/MassScanT1ttttSmooth.root', 'orig/MassScanT1ttttSmooth.root', 'orig/MassScanT1ttttSmooth.root', 'orig/MassScanT1ttttSmooth.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['ExpLim', 'ExpLimSdn', 'ExpLimSup', 'ObsLim', 'ObsLimSdn', 'ObsLimSup', 'MassScan2D'])
T1ttttoff.addMassPlane(T1tttt_1)

#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.checked =''
T1.constraint = "[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription = None
T1.condition = None
T1.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1_1 = T1.addMassPlane(2*[[x, y]])
T1_1.figure = 'Fig.8-c'
T1_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/CMS-PAS-SUS-15-002_Figure_008-c.png'
# T1_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
T1_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
T1_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1qqqqSmooth.root'
T1_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1qqqqSmooth.root'
T1_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/MassScanT1qqqqSmooth.root', 'orig/MassScanT1qqqqSmooth.root', 'orig/MassScanT1qqqqSmooth.root', 'orig/MassScanT1qqqqSmooth.root', 'orig/MassScanT1qqqqSmooth.root', 'orig/MassScanT1qqqqSmooth.root', 'orig/MassScanT1qqqqSmooth.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['ExpLim', 'ExpLimSdn', 'ExpLimSup', 'ObsLim', 'ObsLimSdn', 'ObsLimSup', 'MassScan2D'])



databaseCreator.create()
