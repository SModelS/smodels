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
info = MetaInfoInput('CMS-SUS-16-035')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = '2 SS leptons'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1704.07323'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Eur.Phys.J. C77 (2017) no.9, 578, http://dx.doi.org/10.1140/epjc/s10052-017-5079-z'
info.comment = 'Moriond 2017. Seven SMS interpretations in paper, all implemented.'


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
T1ttttoff.constraint="[[[b,W,b,W]],[[b,W,b,W]]]"
T1ttttoff.condition=None
T1ttttoff.source="CMS"
T1ttttoff.massConstraint=[['dm < 338.0'], ['dm < 338.0']]

#++++++next mass plane block+++++++++

T1tttt_1 = T1tttt.addMassPlane(2*[[x,y]])
T1tttt_1.figure='Fig. 5-a'
T1tttt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-a.png'
T1tttt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-a.root'
T1tttt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-a.root'
T1tttt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-a.root'
T1tttt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-035_Figure_005-a.root','orig/CMS-SUS-16-035_Figure_005-a.root','orig/CMS-SUS-16-035_Figure_005-a.root','orig/CMS-SUS-16-035_Figure_005-a.root','orig/CMS-SUS-16-035_Figure_005-a.root','orig/CMS-SUS-16-035_Figure_005-a.root','orig/CMS-SUS-16-035_Figure_005-a.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['ssexp;1','ssexp_m2s;1','ssexp_p2s;1','ssobs;1','ssobs_m1s;1','ssobs_p1s;1','xsec;1'],units=[None,None,None,None,None,None,'pb'])
T1ttttoff.addMassPlane(T1tttt_1)

#+++++++txName block++++++++++++++++++++
T5tttt=dataset.addTxName('T5tttt')
T5tttt.checked=''
T5tttt.constraint="[[['t'],['t']],[['t'],['t']]]"
T5tttt.condition=None
T5tttt.conditionDescription = None
T5tttt.source="CMS"
T5tttt.massConstraint=None

#++++++next mass plane block+++++++++

T5tttt_1 = T5tttt.addMassPlane(2*[[x,y+175.,y]])
T5tttt_1.figure='Fig. 5-c'
T5tttt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-c.png'
T5tttt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-c.root'
T5tttt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-c.root'
T5tttt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-c.root'
T5tttt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root'],
                       dataFormats=['root','root','root','root','root','root','root'],objectNames=['ssexp;1','ssexp_m1s;1','ssexp_p1s;1','ssobs;1','ssobs_m1s;1','ssobs_p1s;1','xsec;1'],units=[None,None,None,None,None,None,'pb'])

#++++++next mass plane block+++++++++

T5tttt_2 = T5tttt.addMassPlane(2*[[x,y+171.,y]])
T5tttt_2.figure='Fig. 5-c'
T5tttt_2.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-c.png'
T5tttt_2.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-c.root'
T5tttt_2.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-c.root'
T5tttt_2.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-c.root'
T5tttt_2.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root','orig/CMS-SUS-16-035_Figure_005-c.root'],
                       dataFormats=['root','root','root','root','root','root','root'],objectNames=['ssexp;1','ssexp_m1s;1','ssexp_p1s;1','ssobs;1','ssobs_m1s;1','ssobs_p1s;1','xsec;1'],units=[None,None,None,None,None,None,'pb'])

#++++++++txName block+++++++++++++++++

T5WW=dataset.addTxName('T5WW')
T5WW.checked=''
T5WW.constraint="2*([[['jet','jet'],['W+']],[['jet','jet'],['W+']]]+[[['jet','jet'],['W-']],[['jet','jet'],['W-']]])"
T5WW.condition=None
T5WW.conditionDescription = None
T5WW.source="CMS"
T5WW.massConstraint=None

#++++++ mass plane block+++++++++

T5WW_1 = T5WW.addMassPlane(2*[[x,0.5*x+0.5*y,y]])
T5WW_1.figure='Fig. 6-a'
T5WW_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_006-a.png'
T5WW_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_006-a.root'
T5WW_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_006-a.root'
T5WW_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_006-a.root'
T5WW_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-035_Figure_006-a.root','orig/CMS-SUS-16-035_Figure_006-a.root','orig/CMS-SUS-16-035_Figure_006-a.root','orig/CMS-SUS-16-035_Figure_006-a.root','orig/CMS-SUS-16-035_Figure_006-a.root','orig/CMS-SUS-16-035_Figure_006-a.root','orig/CMS-SUS-16-035_Figure_006-a.root'],
                  dataFormats=['root','root','root','root','root','root','root'],objectNames=['ssexp;1','ssexp_m2s;1','ssexp_p2s;1','ssobs;1','ssobs_m1s;1','ssobs_p1s;1','xsec;1'],units=[None,None,None,None,None,None,'pb'])

#+++++++++txName block+++++++++++++++

T5WWoff=dataset.addTxName('T5WWoff')
T5WWoff.checked=''
T5WWoff.constraint="45.35*([[['jet','jet'],['l+','nu']],[['jet','jet'],['l+','nu']]]+[[['jet','jet'],['l-','nu']],[['jet','jet'],[''l-','nu'']]])"
T5WWoff.condition=None
T5WWoff.conditionDescription = None
T5WWoff.source="CMS"
T5WWoff.massConstraint=[['dm>0.0','dm<76.0'],['dm>0','dm<76.0']]

#++++++next mass plane block+++++++++

T5WWoff_1 = T5WWoff.addMassPlane(2*[[x,y+20.0,y]])
T5WWoff_1.figure='Fig. 6-b'
T5WWoff_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_006-b.png'
T5WWoff_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_006-b.root'
T5WWoff_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_006-b.root'
T5WWoff_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_006-b.root'
T5WWoff_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-035_Figure_006-b.root','orig/CMS-SUS-16-035_Figure_006-b.root','orig/CMS-SUS-16-035_Figure_006-b.root','orig/CMS-SUS-16-035_Figure_006-b.root','orig/CMS-SUS-16-035_Figure_006-b.root','orig/CMS-SUS-16-035_Figure_006-b.root','orig/CMS-SUS-16-035_Figure_006-b.root'],
                  dataFormats=['root','root','root','root','root','root','root'],objectNames=['ssexp;1','ssexp_m2s;1','ssexp_p2s;1','ssobs;1','ssobs_m1s;1','ssobs_p1s;1','xsec;1'],units=[None,None,None,None,None,None,'pb'])

#+++++++++txName block+++++++++++++++

T5tctc=dataset.addTxName('T5tctc')
T5tctc.checked=''
T5tctc.constraint="[[['t'],['jet']],[['t'],['jet']]]"
T5tctc.condition=None
T5tctc.conditionDescription = None
T5tctc.source="CMS"
T5tctc.massConstraint=[['dm>169.0','dm>.0'],['dm>169','dm>0.0']]

#++++++next mass plane block+++++++++

T5tctc_1 = T5tctc.addMassPlane(2*[[x,y+20.0,y]])
T5tctc_1.figure='Fig. 5-d'
T5tctc_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-d.png'
T5tctc_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-d.root'
T5tctc_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-d.root'
T5tctc_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-d.root'
T5tctc_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-035_Figure_005-d.root','orig/CMS-SUS-16-035_Figure_005-d.root','orig/CMS-SUS-16-035_Figure_005-d.root','orig/CMS-SUS-16-035_Figure_005-d.root','orig/CMS-SUS-16-035_Figure_005-d.root','orig/CMS-SUS-16-035_Figure_005-d.root','orig/CMS-SUS-16-035_Figure_005-d.root'],
                  dataFormats=['root','root','root','root','root','root','root'],objectNames=['ssexp;1','ssexp_m1s;1','ssexp_p1s;1','ssobs;1','ssobs_m1s;1','ssobs_p1s;1','xsec;1'],units=[None,None,None,None,None,None,'pb'])

#+++++++++txName block+++++++++++++++

T5ttbbWWoff=dataset.addTxName('T5ttbbWWoff')
T5ttbbWWoff.checked=''
T5ttbbWWoff.constraint="[[['t','b'],['jet','jet']],[['t','b'],['jet','jet']]]"
T5ttbbWWoff.condition=None
T5ttbbWWoff.conditionDescription = None
T5ttbbWWoff.source="CMS"
T5ttbbWWoff.massConstraint=[['dm>169.0','dm<76.0'],['dm>169','dm<76.0']]

#++++++next mass plane block+++++++++

T5ttbbWWoff_1 = T5ttbbWWoff.addMassPlane(2*[[x,y+5.0,y]])
T5ttbbWWoff_1.figure='Fig. 5-b'
T5ttbbWWoff_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-b.png'
T5ttbbWWoff_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-b.root'
T5ttbbWWoff_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-b.root'
T5ttbbWWoff_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_005-b.root'
T5ttbbWWoff_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-035_Figure_005-b.root','orig/CMS-SUS-16-035_Figure_005-b.root','orig/CMS-SUS-16-035_Figure_005-b.root','orig/CMS-SUS-16-035_Figure_005-b.root','orig/CMS-SUS-16-035_Figure_005-b.root','orig/CMS-SUS-16-035_Figure_005-b.root','orig/CMS-SUS-16-035_Figure_005-b.root'],
                  dataFormats=['root','root','root','root','root','root','root'],objectNames=['ssexp;1','ssexp_m2s;1','ssexp_p2s;1','ssobs;1','ssobs_m1s;1','ssobs_p1s;1','xsec;1'],units=[None,None,None,None,None,None,'pb'])

#+++++++++txName block+++++++++++++++

T6ttWW=dataset.addTxName('T6ttWW')
T6ttWW.checked=''
T6ttWW.constraint="[[['t'],['W']],[['t'],['W']]]"
T6ttWW.condition=None
T6ttWW.conditionDescription = None
T6ttWW.source="CMS"
T6ttWW.massConstraint=[['dm>169.0','dm>76.0'],['dm>169','dm>76.0']]

T6ttoffWW=dataset.addTxName('T6ttoffWW')
T6ttoffWW.checked=''
T6ttoffWW.constraint="[[['b','W'],['W']],[['b','W'],['W']]]"
T6ttoffWW.condition=None
T6ttoffWW.conditionDescription = None
T6ttoffWW.source="CMS"
T6ttoffWW.massConstraint=[['76<dm<169.0','dm>76.0'],['76<dm<169','dm>76.0']]

T6ttWWoff=dataset.addTxName('T6ttWWoff')
T6ttWWoff.checked=''
T6ttWWoff.constraint="2.228*([[['t'],['jet','jet']],[['t'],['jet','jet']]])"
T6ttWWoff.condition=None
T6ttWWoff.conditionDescription = None
T6ttWWoff.source="CMS"
T6ttWWoff.massConstraint=[['dm>169.0','dm<76.0'],['dm>169','dm<76.0']]

#++++++next mass plane block+++++++++

T6ttWW_1 = T6ttWW.addMassPlane(2*[[x,y,50]])
T6ttWW_1.figure='Fig. 7'
T6ttWW_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_007.png'
T6ttWW_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_007.root'
T6ttWW_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_007.root'
T6ttWW_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/CMS-SUS-16-035_Figure_007.root'
T6ttWW_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-035_Figure_007.root','orig/CMS-SUS-16-035_Figure_007.root','orig/CMS-SUS-16-035_Figure_007.root','orig/CMS-SUS-16-035_Figure_007.root','orig/CMS-SUS-16-035_Figure_007.root','orig/CMS-SUS-16-035_Figure_007.root','orig/CMS-SUS-16-035_Figure_007.root'],
                  dataFormats=['root','root','root','root','root','root','root'],objectNames=['ssexp;1','ssexp_m1s;1','ssexp_p1s;1','ssobs;1','ssobs_m1s;1','ssobs_p1s;1','xsec;1'],units=[None,None,None,None,None,None,'pb'])

T6ttoffWW.addMassPlane(T6ttWW_1);
T6ttWWoff.addMassPlane(T6ttWW_1);

databaseCreator.create()
