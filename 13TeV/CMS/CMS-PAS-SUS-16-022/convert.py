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
info = MetaInfoInput('CMS-PAS-SUS-16-022')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-022/'
info.sqrts = 13
info.lumi = 12.9
info.prettyName = '>= 3 leptons'
info.comment = 'https://cds.cern.ch/record/2205165'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked =''
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription =None
T1tttt.condition =None
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked =''
T1ttttoff.constraint = "[[['b','W','b','W']],[['b','W','b','W']]]" 
T1ttttoff.conditionDescription =None
T1ttttoff.condition =None
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure = 'Figure 5-a'
T1tttt_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-022/CMS-PAS-SUS-16-022_Figure_005-a.png'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-022/CMS-PAS-SUS-16-022_Figure_005-a.root'
T1tttt_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-022/CMS-PAS-SUS-16-022_Figure_005-a.root'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-022/CMS-PAS-SUS-16-022_Figure_005-a.root'
T1tttt_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-022/CMS-PAS-SUS-16-022_Figure_005-a.root' 
T1tttt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/CMS-PAS-SUS-16-022_Figure_005-a.root', 'orig/CMS-PAS-SUS-16-022_Figure_005-a.root', 'orig/CMS-PAS-SUS-16-022_Figure_005-a.root', 'orig/CMS-PAS-SUS-16-022_Figure_005-a.root', 'orig/CMS-PAS-SUS-16-022_Figure_005-a.root', 'orig/CMS-PAS-SUS-16-022_Figure_005-a.root', 'orig/CMS-PAS-SUS-16-022_Figure_005-a.root', 'orig/CMS-PAS-SUS-16-022_Figure_005-a.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['gr_exp_smoothed', 'gr_em1s_smoothed', 'gr_ep1s_smoothed', 'exp_xs0', 'gr_obs_smoothed', 'gr_om1s_smoothed', 'gr_op1s_smoothed', 'obs_xs0'])
T1ttttoff.addMassPlane(T1tttt_1)



databaseCreator.create()
