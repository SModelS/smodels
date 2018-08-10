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

databaseCreator.ncpus = 1


#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-15-008')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-008/index.html'
info.sqrts = 13
info.lumi = 2.3
info.prettyName = 'same-sign dilepton'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1605.03171'
info.publication = 'http://dx.doi.org/10.1140/epjc/s10052-016-4261-z'
info.comment = 'Only T1tttt can be implemneted (other topologies have only up to 2 mass planes). Superseded by CMS-SUS-16-035.'
info.supersededBy = 'CMS-SUS-16-035'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure = 'Fig.4a'
T1tttt_1.figureUrl ='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-008/CMS-SUS-15-008_Figure_004-a.png'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-008/CMS-SUS-15-008_Figure_004-a.root'
T1tttt_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-008/CMS-SUS-15-008_Figure_004-a.root'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-008/CMS-SUS-15-008_Figure_004-a.root'
T1tttt_1.exclusionDataUrl ='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-008/CMS-SUS-15-008_Figure_004-a.root'
T1tttt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/CMS-SUS-15-008_Figure_004-a.root', 'orig/CMS-SUS-15-008_Figure_004-a.root', 'orig/CMS-SUS-15-008_Figure_004-a.root', 'orig/CMS-SUS-15-008_Figure_004-a.root', 'orig/CMS-SUS-15-008_Figure_004-a.root', 'orig/CMS-SUS-15-008_Figure_004-a.root', 'orig/CMS-SUS-15-008_Figure_004-a.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['ssexp', 'ssexp_m1s', 'ssexp_p1s', 'ssobs', 'ssobs_m1s', 'ssobs_p1s', 'xsec'])
T1ttttoff.addMassPlane(T1tttt_1)



databaseCreator.create()
