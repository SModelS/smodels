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
info = MetaInfoInput('CMS-SUS-16-041')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'Multileptons + jets + Etmiss'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1710.09154'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Submitted to J. High Energy Phys.'
info.comment = 'Moriond 2017. Omitted the gluino -> q qbar W/Z decay simplified model since the branching ratio information for the gluino to the chargino or neutralino2 was not provided.'


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
T1ttttoff.massConstraint=[['160 < dm < 338.0'], ['160 < dm < 338.0']]

#++++++next mass plane block+++++++++

T1tttt_1 = T1tttt.addMassPlane(2*[[x,y]])
T1tttt_1.figure='Fig. 5-a'
T1tttt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_005-a.png'
T1tttt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_005-a.root'
T1tttt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_005-a.root'
T1tttt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_005-a.root'
T1tttt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-041_Figure_005-a.root','orig/CMS-SUS-16-041_Figure_005-a.root','orig/CMS-SUS-16-041_Figure_005-a.root','orig/CMS-SUS-16-041_Figure_005-a.root','orig/CMS-SUS-16-041_Figure_005-a.root','orig/CMS-SUS-16-041_Figure_005-a.root','orig/CMS-SUS-16-041_Figure_005-a.root'],
                    dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1'],indices=[3,5, 4, 6, 8,7, 2],units=[None,None,None,None,None,None,'pb'])
T1ttttoff.addMassPlane(T1tttt_1)

#+++++txName block +++++++++++++++++

T6ttWW=dataset.addTxName('T6ttWW')
T6ttWW.checked=''
T6ttWW.constraint="[[['t'],['W']],[['t'],['W']]]"
T6ttWW.condition=None
T6ttWW.conditionDescription = None
T6ttWW.source="CMS"
T6ttWW.massConstraint=None

T6ttWWoff=dataset.addTxName('T6ttWWoff')
T6ttWWoff.checked=''
T6ttWWoff.constraint="20.66*[[['t'],['l','nu']],[['t'],['l','nu']]]"
T6ttWWoff.condition=None
T6ttWWoff.conditionDescription = None
T6ttWWoff.source="CMS"
T6ttWWoff.massConstraint=[['dm>169.0','dm<76.0'],['dm>169','dm<76.0']]

T6ttoffWW=dataset.addTxName('T6ttoffWW')
T6ttoffWW.checked=''
T6ttoffWW.constraint="[[['b','W'],['W']],[['b','W'],['W']]]"
T6ttoffWW.condition=None
T6ttoffWW.conditionDescription = None
T6ttoffWW.source="CMS"
T6ttoffWW.massConstraint=[['76<dm<169.0','dm>76.0'],['76<dm<169','dm>76.0']]

#++++++next mass plane block+++++++++

T6ttWW_1 = T6ttWW.addMassPlane(2*[[x,y,50.0]])
T6ttWW_1.figure='Fig. 6'
T6ttWW_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_006.png'
T6ttWW_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_006.root'
T6ttWW_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_006.root'
T6ttWW_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_006.root'
T6ttWW_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-041_Figure_006.root','orig/CMS-SUS-16-041_Figure_006.root','orig/CMS-SUS-16-041_Figure_006.root','orig/CMS-SUS-16-041_Figure_006.root','orig/CMS-SUS-16-041_Figure_006.root','orig/CMS-SUS-16-041_Figure_006.root','orig/CMS-SUS-16-041_Figure_006.root'],
                  dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1'],indices=[3,5, 4, 6, 8,7, 2],units=[None,None,None,None,None,None,'pb'])

T6ttWWoff.addMassPlane(T6ttWW_1)
T6ttoffWW.addMassPlane(T6ttWW_1)

#+++++txName block +++++++++++++++++

T6ZZtt=dataset.addTxName('T6ZZtt')
T6ZZtt.checked=''
T6ZZtt.constraint="[[['Z'],['t']],[['Z'],['t']]]"
T6ZZtt.condition=None
T6ZZtt.conditionDescription = None
T6ZZtt.source="CMS"
T6ZZtt.massConstraint=None

T6ZZofftt=dataset.addTxName('T6ZZofftt')
T6ZZofftt.checked=''
T6ZZofftt.constraint="221.44*([[['l+', 'l-'],['t']],[['l+', 'l-'],['t']]])"
T6ZZofftt.condition=None
T6ZZofftt.conditionDescription = None
T6ZZofftt.source="CMS"
T6ZZofftt.massConstraint=[['dm<86.0','dm>169.0'],['dm<86','dm>169.0']]
#++++++next mass plane block+++++++++

T6ZZtt_1 = T6ZZtt.addMassPlane(2*[[x,y,y-175.0]])
T6ZZtt_1.figure='Fig. 7-c'
T6ZZtt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_007-c.png'
T6ZZtt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_007-c.root'
T6ZZtt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_007-c.root'
T6ZZtt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_007-c.root'
T6ZZtt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-041_Figure_007-c.root','orig/CMS-SUS-16-041_Figure_007-c.root','orig/CMS-SUS-16-041_Figure_007-c.root','orig/CMS-SUS-16-041_Figure_007-c.root','orig/CMS-SUS-16-041_Figure_007-c.root','orig/CMS-SUS-16-041_Figure_007-c.root','orig/CMS-SUS-16-041_Figure_007-c.root'],
                  dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1'],indices=[3,5, 4, 6, 8,7, 2],units=[None,None,None,None,None,None,'pb'])

T6ZZofftt.addMassPlane(T6ZZtt_1)

#+++++txName block +++++++++++++++++

T6HHtt=dataset.addTxName('T6HHtt')
T6HHtt.checked=''
T6HHtt.constraint="[[['higgs'],['t']],[['higgs'],['t']]]"
T6HHtt.condition=None
T6HHtt.conditionDescription = None
T6HHtt.source="CMS"
T6HHtt.massConstraint=None

#++++++next mass plane block+++++++++

T6HHtt_1 = T6HHtt.addMassPlane(2*[[x,y,y-175.0]])
T6HHtt_1.figure='Fig. 7-a'
T6HHtt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_007-a.png'
T6HHtt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_007-a.root'
T6HHtt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_007-a.root'
T6HHtt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-041/CMS-SUS-16-041_Figure_007-a.root'
T6HHtt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-041_Figure_007-a.root','orig/CMS-SUS-16-041_Figure_007-a.root','orig/CMS-SUS-16-041_Figure_007-a.root','orig/CMS-SUS-16-041_Figure_007-a.root','orig/CMS-SUS-16-041_Figure_007-a.root','orig/CMS-SUS-16-041_Figure_007-a.root','orig/CMS-SUS-16-041_Figure_007-a.root'],
                  dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1'],indices=[3,5, 4, 6, 8,7, 2],units=[None,None,None,None,None,None,'pb'])


databaseCreator.create()
