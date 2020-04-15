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
info = MetaInfoInput('CMS-SUS-17-005')
info.url = 'https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'multijets + Etmiss, top tagging'
info.private = False
info.arxiv = ''
info.contact = ''
info.publication = ''
info.comment = ''



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2bbWWoff = dataset.addTxName('T2bbWWoff')
T2bbWWoff.checked = ''
T2bbWWoff.constraint = "[[['b','q','q']],[['b','q','q']]]"
T2bbWWoff.conditionDescription = None
T2bbWWoff.condition = None
T2bbWWoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2bbWWoff_1 = T2bbWWoff.addMassPlane(2*[[x, x-y]])
T2bbWWoff_1.figure = 'Fig.7-b'
T2bbWWoff_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_007-b.png'
T2bbWWoff_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_007-b.root'
T2bbWWoff_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_007-b.root'
T2bbWWoff_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_007-b.root'
T2bbWWoff_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_007-b.root'
T2bbWWoff_1.setSources(dataLabels= ['obsExclusionM1','expExclusionM1', 'expExclusionP1','obsExclusion', 'expExclusion','upperLimits','expectedUpperLimits','obsExclusionP1'],
                 dataFiles= ['orig/CMS-SUS-17-005_Figure_007-b.root', 'orig/CMS-SUS-17-005_Figure_007-b.root','orig/CMS-SUS-17-005_Figure_007-b.root', 
                 'orig/CMS-SUS-17-005_Figure_007-b.root','orig/CMS-SUS-17-005_Figure_007-b.root', 'orig/CMS-SUS-17-005_Figure_007-b.root','orig/CMS-SUS-17-005_Figure_007-b.root', 
                 'orig/CMS-SUS-17-005_Figure_007-b.root'],
                 dataFormats= ['root','root','root','root','root','root','root','root'],objectNames= ['contour_ObsMinus1','contour_ExpMinus1','contour_ExpPlus1','contour_Obs','contour_Exp','xsecUL_Obs','xsecUL_Exp','contour_ObsPlus1'])
                 
                 
                 
	
  	
  	

                 
                 
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
T1tttt_1.figure = 'Fig.13-a'
T1tttt_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-a.png'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-a.root'
T1tttt_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-a.root'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-a.root'
T1tttt_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-a.root'
T1tttt_1.setSources(dataLabels= ['upperLimits','expectedUpperLimits'],
                 dataFiles= ['orig/CMS-SUS-17-005_Figure_013-a.root', 'orig/CMS-SUS-17-005_Figure_013-a.root'],
                 dataFormats= ['root', 'root'],objectNames= ['combined_obsLimit_BR100pct', 'combined_expLimit_BR100pct'])
T1ttttoff.addMassPlane(T1tttt_1)

#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T5tctc')
T1.checked =''
T1.constraint = "[[['t'],['c']],[['t'],['c']]]"
T1.conditionDescription = None
T1.condition = None
T1.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1qqqq_1 = T1.addMassPlane(2*[[x,y+20, y]])
T1qqqq_1.figure = 'Fig.13-b'
T1qqqq_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-b.png'
T1qqqq_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-b.root'
T1qqqq_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-b.root'
T1qqqq_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-16-009_Figure_013-b.root'
T1qqqq_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-b.root'
T1qqqq_1.setSources(dataLabels= ['upperLimits','expectedUpperLimits'],
                 dataFiles= ['orig/CMS-SUS-17-005_Figure_013-b.root', 'orig/CMS-SUS-17-005_Figure_013-b.root'],
                 dataFormats= ['root', 'root'],objectNames= ['combined_obsLimit_BR100pct', 'combined_expLimit_BR100pct'])



databaseCreator.create()
