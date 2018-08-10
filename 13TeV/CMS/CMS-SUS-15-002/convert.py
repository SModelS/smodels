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

databaseCreator.ncpus = 1


#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-15-002')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/'
info.sqrts = 13
info.lumi = 2.2
info.prettyName = 'multijets + Etmiss, HT'
info.private = False
info.arxiv = ''
info.contact = ''
info.publication = ''
info.comment = 'Do not know where the off constraint for T1tttt came from. Superseded by CMS-SUS-16-033 and CMS-SUS-16-036.'
info.supersedes = 'CMS-SUS-PAS-15-002'
info.supersededBy = 'CMS-SUS-16-033'


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
T1bbbb_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008-a.png'
T1bbbb_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1bbbb_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1bbbb_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1bbbb_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1bbbb_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['ExpLim_T1bbbb', 'ExpLimMinus_T1bbbb', 'ExpLimPlus_T1bbbb', 'T1bbbb_MassScan2DExp', 'ObsLim_T1bbbb', 'ObsLimMinus_T1bbbb', 'ObsLimPlus_T1bbbb', 'T1bbbb_MassScan2DObs'])

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
T1tttt_1.figure = 'Fig.9-a'
T1tttt_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008-a.png'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1tttt_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1tttt_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1tttt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['ExpLim_T1tttt', 'ExpLimMinus_T1tttt', 'ExpLimPlus_T1tttt', 'T1tttt_MassScan2DExp', 'ObsLim_T1tttt', 'ObsLimMinus_T1tttt', 'ObsLimPlus_T1tttt', 'T1tttt_MassScan2DObs'])
T1ttttoff.addMassPlane(T1tttt_1)

#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.checked =''
T1.constraint = "[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription = None
T1.condition = None
T1.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1qqqq_1 = T1.addMassPlane(2*[[x, y]])
T1qqqq_1.figure = 'Fig.8-c'
T1qqqq_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008-a.png'
T1qqqq_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1qqqq_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1qqqq_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1qqqq_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-15-002_Figure_008.root'
T1qqqq_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root', 'orig/CMS-SUS-15-002_Figure_008.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['ExpLim_T1qqqq', 'ExpLimMinus_T1qqqq', 'ExpLimPlus_T1qqqq', 'T1qqqq_MassScan2DExp', 'ObsLim_T1qqqq', 'ObsLimMinus_T1qqqq', 'ObsLimPlus_T1qqqq', 'T1qqqq_MassScan2DObs'])



databaseCreator.create()
