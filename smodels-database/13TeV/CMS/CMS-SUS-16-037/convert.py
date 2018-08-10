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
info = MetaInfoInput('CMS-SUS-16-037')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-037/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = '1L + jets + Etmiss with MJ (sum of masses of large radius jets)'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1705.04673'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Phys. Rev. Lett. 119 (2017) 151802, http://dx.doi.org/10.1103/PhysRevLett.119.151802'
info.comment = ' '


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++txName block +++++++++++++++++

T1tttt=dataset.addTxName('T1tttt')
T1tttt.checked=''
T1tttt.constraint="[[['t','t']],[['t','t']]]"
T1tttt.condition=None
T1tttt.conditionDescription = None
T1tttt.source="CMS"
T1tttt.massConstraint = None

T1ttttoff=dataset.addTxName('T1ttttoff')
T1ttttoff.checked=''
T1ttttoff.constraint="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.condition=None
T1ttttoff.source="CMS"
T1ttttoff.massConstraint=[['160 < dm < 338.0'], ['160 < dm < 338.0']] 

#++++++++++mass plane++++++++++++++++++

T1tttt_1 = T1tttt.addMassPlane(2*[[x,y]])
T1tttt_1.figure='Fig. 2'
T1tttt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-037/CMS-SUS-16-037_Figure_002.png'
T1tttt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-037/CMS-SUS-16-037_Figure_002.root'
T1tttt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-037/CMS-SUS-16-037_Figure_002.root'
T1tttt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-037/CMS-SUS-16-037_Figure_002.root'
T1tttt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-037_Figure_002.root','orig/CMS-SUS-16-037_Figure_002.root','orig/CMS-SUS-16-037_Figure_002.root','orig/CMS-SUS-16-037_Figure_002.root','orig/CMS-SUS-16-037_Figure_002.root','orig/CMS-SUS-16-037_Figure_002.root','orig/CMS-SUS-16-037_Figure_002.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['T1ttttExpectedLimit;1','T1ttttExpectedLimitDown;1','T1ttttExpectedLimitUp;1','T1ttttObservedLimit;1','T1ttttObservedLimitDown;1','T1ttttObservedLimitUp;1','T1ttttObservedExcludedXsec;1'],
                    units=[None,None,None,None,None,None,'pb'])
T1ttttoff.addMassPlane(T1tttt_1)
#++++++++txName block++++++++++++++++++

T5tttt=dataset.addTxName('T5tttt')
T5tttt.checked=''
T5tttt.constraint="[[['t'],['t']],[['t'],['t']]]"
T5tttt.condition=None
T5tttt.conditionDescription = None
T5tttt.source="CMS"
T5tttt.massConstraint=None

T5ttofftt=dataset.addTxName('T5ttofftt')
T5ttofftt.checked=''
T5ttofftt.constraint="[[['b','W'],['t']],[['b','W'],['t']]]"
T5ttofftt.condition=None
T5ttofftt.conditionDescription=None
T5ttofftt.source="CMS"
T5ttofftt.massConstraint=[['80<dm<169.0','dm>169.0'],['80<dm<169.0','dm>169.0']]

#++++++++++mass plane++++++++++++++++++

T5tttt_1 = T5tttt.addMassPlane(2*[[x,y+175.0,y]])
T5tttt_1.figure='aux. Fig. 5'
T5tttt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-037/CMS-SUS-16-037_Figure-aux_005.png'
T5tttt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-037/CMS-SUS-16-037_Figure-aux_005.root'
T5tttt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-037/CMS-SUS-16-037_Figure-aux_005.root'
T5tttt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-037/CMS-SUS-16-037_Figure-aux_005.root'
T5tttt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                  dataFiles=['orig/CMS-SUS-16-037_Figure-aux_005.root','orig/CMS-SUS-16-037_Figure-aux_005.root','orig/CMS-SUS-16-037_Figure-aux_005.root','orig/CMS-SUS-16-037_Figure-aux_005.root','orig/CMS-SUS-16-037_Figure-aux_005.root','orig/CMS-SUS-16-037_Figure-aux_005.root','orig/CMS-SUS-16-037_Figure-aux_005.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['T5ttttExpectedLimit;1','T5ttttExpectedLimitDown;1','T5ttttExpectedLimitUp;1','T5ttttObservedLimit;1',
'T5ttttObservedLimitDown;1','T5ttttObservedLimitUp;1','T5ttttObservedExcludedXsec;1'],
                     units=[None,None,None,None,None,None,'pb'])


T5ttofftt.addMassPlane(T5tttt_1)



databaseCreator.create()
