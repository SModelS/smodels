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
type = types.StringType)
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
help = 'path to the package smodels_utils',\
type = types.StringType)
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
info = MetaInfoInput('CMS-SUS-16-042')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-042/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = '1L + jets + Etmiss (with Delta Phi)'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1709.09814'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Submitted to Phys. Lett. B'
info.comment = 'Moriond 2017.'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++txName block++++++++++++++++++++
T1tttt=dataset.addTxName('T1tttt')
T1tttt.checked=''
T1tttt.constraint="[[['t','t']],[['t','t']]]"
T1tttt.condition=None
T1tttt.conditionDescription = None
T1tttt.source="CMS"
T1tttt.massConstraint=[['dm>=338.0'],['dm>=338.0']]

T1ttttoff=dataset.addTxName('T1ttttoff')
T1ttttoff.checked=''
T1ttttoff.constraint="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.condition=None
T1ttttoff.source="CMS"
T1ttttoff.massConstraint=[['160 <dm < 338.0'], ['160 < dm < 338.0']]

#++++++next mass plane block+++++++++

T1tttt_1 = T1tttt.addMassPlane(2*[[x,y]])
T1tttt_1.figure='Fig. 5-a'
T1tttt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-042/CMS-SUS-16-042_Figure_005-a.png'
T1tttt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-042/CMS-SUS-16-042_Figure_005-a.root'
T1tttt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-042/CMS-SUS-16-042_Figure_005-a.root'
T1tttt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-042/CMS-SUS-16-042_Figure_005-a.root'
T1tttt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-042_Figure_005-a.root','orig/CMS-SUS-16-042_Figure_005-a.root','orig/CMS-SUS-16-042_Figure_005-a.root','orig/CMS-SUS-16-042_Figure_005-a.root','orig/CMS-SUS-16-042_Figure_005-a.root','orig/CMS-SUS-16-042_Figure_005-a.root','orig/CMS-SUS-16-042_Figure_005-a.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['T1ttttExpectedLimit;1','T1ttttExpectedLimitDown;1','T1ttttExpectedLimitUp;1','T1ttttObservedLimit;1','T1ttttObservedLimitDown;1','T1ttttObservedLimitUp;1','T1ttttObservedExcludedXsec;1'],units=[None,None,None,None,None,None,'pb'])

T1ttttoff.addMassPlane(T1tttt_1)

#++++++++txName block+++++++++++++++++

T5WW=dataset.addTxName('T5WW')
T5WW.checked=''
T5WW.constraint="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.condition=None
T5WW.conditionDescription = None
T5WW.source="CMS"
T5WW.massConstraint=None


T5WWoff=dataset.addTxName('T5WWoff')
T5WWoff.checked=''
T5WWoff.constraint="7.1*([[['jet','jet'],['l','nu']],[['jet','jet'],['jet','jet']]])"
T5WWoff.condition=None
T5WWoff.conditionDescription = None
T5WWoff.source="CMS"
T5WWoff.massConstraint=[['dm>0','dm<76.0'],['dm>0.','dm<76.0']]


#++++++next mass plane block+++++++++

T5WW_1 = T5WW.addMassPlane(2*[[x,0.5*x+0.5*y,y]])
T5WW_1.figure='Fig. 5-b'
T5WW_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-042/CMS-SUS-16-042_Figure_005-b.png'
T5WW_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-042/CMS-SUS-16-042_Figure_005-b.root'
T5WW_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-042/CMS-SUS-16-042_Figure_005-b.root'
T5WW_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-042/CMS-SUS-16-042_Figure_005-b.root'
T5WW_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-042_Figure_005-b.root','orig/CMS-SUS-16-042_Figure_005-b.root','orig/CMS-SUS-16-042_Figure_005-b.root','orig/CMS-SUS-16-042_Figure_005-b.root','orig/CMS-SUS-16-042_Figure_005-b.root','orig/CMS-SUS-16-042_Figure_005-b.root','orig/CMS-SUS-16-042_Figure_005-b.root'],
                  dataFormats=['root','root','root','root','root','root','root'],objectNames=['T5qqqqWWExpectedLimit;1','T5qqqqWWExpectedLimitDown;1','T5qqqqWWExpectedLimitUp;1','T5qqqqWWObservedLimit;1','T5qqqqWWObservedLimitDown;1','T5qqqqWWObservedLimitUp;1','T5qqqqWWObservedExcludedXsec;1'],units=[None,None,None,None,None,None,'pb'])

T5WWoff.addMassPlane(T5WW_1)

databaseCreator.create()
