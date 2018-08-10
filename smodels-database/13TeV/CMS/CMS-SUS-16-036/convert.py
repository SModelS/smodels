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
info = MetaInfoInput('CMS-SUS-16-036')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = '0L + jets + Etmiss (using MT2)'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1705.04650'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Eur. Phys. J. C 77 (2017) 710, http://dx.doi.org/10.1140/epjc/s10052-017-5267-x'
info.comment = 'Moriond 2017. The mixed decay topology for stop pair production with stop->t \chi1_0 (50% BR) and stop->b \chi^{\pm} -> b W^{\pm} \chi1_0 (50%), has not been implemented since the branching ratio is fixed and the mass difference between the chargino and neutralino is fixed to 5 GeV. '


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++txName block +++++++++++++++++

T1bbbb=dataset.addTxName('T1bbbb')
T1bbbb.checked=''
T1bbbb.constraint="[[['b','b']],[['b','b']]]"
T1bbbb.condition=None
T1bbbb.conditionDescription = None
T1bbbb.source="CMS"

#++++++next mass plane block+++++++++

T1bbbb_1 = T1bbbb.addMassPlane(2*[[x,y]])
T1bbbb_1.figure='Fig. 6-a'
T1bbbb_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-a.png'
T1bbbb_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-a.root'
T1bbbb_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-a.root'
T1bbbb_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-a.root'
T1bbbb_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-036_Figure_006-a.root','orig/CMS-SUS-16-036_Figure_006-a.root','orig/CMS-SUS-16-036_Figure_006-a.root','orig/CMS-SUS-16-036_Figure_006-a.root','orig/CMS-SUS-16-036_Figure_006-a.root','orig/CMS-SUS-16-036_Figure_006-a.root','orig/CMS-SUS-16-036_Figure_006-a.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['ExpLim;1','ExpLim_1Sdn;1','ExpLim_1Sup;1','ObsLim;1','ObsLim_1Sdn;1','ObsLim_1Sup;1','ObsXS_2D;1'],units=[None,None,None,None,None,None,'pb'])

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
T1tttt_1.figure='Fig. 6-b'
T1tttt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-b.png'
T1tttt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-b.root'
T1tttt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-b.root'
T1tttt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-b.root'
T1tttt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-036_Figure_006-b.root','orig/CMS-SUS-16-036_Figure_006-b.root','orig/CMS-SUS-16-036_Figure_006-b.root','orig/CMS-SUS-16-036_Figure_006-b.root','orig/CMS-SUS-16-036_Figure_006-b.root','orig/CMS-SUS-16-036_Figure_006-b.root','orig/CMS-SUS-16-036_Figure_006-b.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['ExpLim;1','ExpLim_1Sdn;1','ExpLim_1Sup;1','ObsLim;1','ObsLim_1Sdn;1','ObsLim_1Sup;1','ObsXS_2D;1'],units=[None,None,None,None,None,None,'pb'])
T1ttttoff.addMassPlane(T1tttt_1)

#+++++txName block +++++++++++++++++

T1=dataset.addTxName('T1')
T1.checked=''
T1.constraint="[[['jet','jet']],[['jet','jet']]]"
T1.condition=None
T1.conditionDescription = None
T1.source="CMS"

#++++++next mass plane block+++++++++

T1_1 = T1.addMassPlane(2*[[x,y]])
T1_1.figure='Fig. 6-c'
T1_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-c.png'
T1_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-c.root'
T1_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-c.root'
T1_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_006-c.root'
T1_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-036_Figure_006-c.root','orig/CMS-SUS-16-036_Figure_006-c.root','orig/CMS-SUS-16-036_Figure_006-c.root','orig/CMS-SUS-16-036_Figure_006-c.root','orig/CMS-SUS-16-036_Figure_006-c.root','orig/CMS-SUS-16-036_Figure_006-c.root','orig/CMS-SUS-16-036_Figure_006-c.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['ExpLim;1','ExpLim_1Sdn;1','ExpLim_1Sup;1','ObsLim;1','ObsLim_1Sdn;1','ObsLim_1Sup;1','ObsXS_2D;1'],units=[None,None,None,None,None,None,'pb'])

#+++++txName block +++++++++++++++++

T2=dataset.addTxName('T2')
T2.checked=''
T2.constraint="[[['jet']],[['jet']]]"
T2.condition=None
T2.conditionDescription = None
T2.source="CMS"

#++++++next mass plane block+++++++++

T2_1 = T2.addMassPlane(2*[[x,y]])
T2_1.figure='Fig. 7-c'
T2_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-c.png'
T2_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-c.root'
T2_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-c.root'
T2_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-c.root'
T2_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-036_Figure_007-c.root','orig/CMS-SUS-16-036_Figure_007-c.root','orig/CMS-SUS-16-036_Figure_007-c.root','orig/CMS-SUS-16-036_Figure_007-c.root','orig/CMS-SUS-16-036_Figure_007-c.root','orig/CMS-SUS-16-036_Figure_007-c.root','orig/CMS-SUS-16-036_Figure_007-c.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['ExpLim;1','ExpLim_1Sdn;1','ExpLim_1Sup;1','ObsLim;1','ObsLim_1Sdn;1','ObsLim_1Sup;1','ObsXS_2D;1'],units=[None,None,None,None,None,None,'pb'])

#+++++txName block +++++++++++++++++

T2bb=dataset.addTxName('T2bb')
T2bb.checked=''
T2bb.constraint="[[['b']],[['b']]]"
T2bb.condition=None
T2bb.conditionDescription = None
T2bb.source="CMS"

#++++++next mass plane block+++++++++

T2bb_1 = T2bb.addMassPlane(2*[[x,y]])
T2bb_1.figure='Fig. 7-a'
T2bb_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-a.png'
T2bb_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-a.root'
T2bb_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-a.root'
T2bb_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-a.root'
T2bb_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-036_Figure_007-a.root','orig/CMS-SUS-16-036_Figure_007-a.root','orig/CMS-SUS-16-036_Figure_007-a.root','orig/CMS-SUS-16-036_Figure_007-a.root','orig/CMS-SUS-16-036_Figure_007-a.root','orig/CMS-SUS-16-036_Figure_007-a.root','orig/CMS-SUS-16-036_Figure_007-a.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['ExpLim;1','ExpLim_1Sdn;1','ExpLim_1Sup;1','ObsLim;1','ObsLim_1Sdn;1','ObsLim_1Sup;1','ObsXS_2D;1'],units=[None,None,None,None,None,None,'pb'])


#+++++txName block +++++++++++++++++

T2tt=dataset.addTxName('T2tt')
T2tt.checked=''
T2tt.constraint="[[['t']],[['t']]]"
T2tt.condition=None
T2tt.conditionDescription = None
T2tt.massConstraint=[['dm>169'],['dm>169']]
T2tt.source="CMS"

T2ttoff=dataset.addTxName('T2ttoff')
T2ttoff.checked=''
T2ttoff.constraint="[[['b','W']],[['b','W']]]"
T2ttoff.condition=None
T2ttoff.conditionDescription=None
T2ttoff.massConstraint=[['80<dm<169'],['80<dm<169']]
T2ttoff.source="CMS"
#++++++next mass plane block+++++++++

T2tt_1 = T2tt.addMassPlane(2*[[x,y]])
T2tt_1.figure='Fig. 7-b'
T2tt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-b.png'
T2tt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-b.root'
T2tt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-b.root'
T2tt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_007-b.root'
T2tt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-036_Figure_007-b.root','orig/CMS-SUS-16-036_Figure_007-b.root','orig/CMS-SUS-16-036_Figure_007-b.root','orig/CMS-SUS-16-036_Figure_007-b.root','orig/CMS-SUS-16-036_Figure_007-b.root','orig/CMS-SUS-16-036_Figure_007-b.root','orig/CMS-SUS-16-036_Figure_007-b.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['ExpLim;1','ExpLim_1Sdn;1','ExpLim_1Sup;1','ObsLim;1','ObsLim_1Sdn;1','ObsLim_1Sup;1','ObsXS_2D;1'],units=[None,None,None,None,None,None,'pb'])

T2ttoff.addMassPlane(T2tt_1)
#+++++txName block +++++++++++++++++

T6bbWW=dataset.addTxName('T6bbWW')
T6bbWW.checked=''
T6bbWW.constraint="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.condition=None
T6bbWW.conditionDescription = None
T6bbWW.source="CMS"
T6bbWW.massConstraint=None
#++++++next mass plane block+++++++++

T6bbWW_1 = T6bbWW.addMassPlane(2*[[x,0.5*x+0.5*y,y]])
T6bbWW_1.figure='Fig. 8-a'
T6bbWW_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_008-a.png'
T6bbWW_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_008-a.root'
T6bbWW_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_008-a.root'
T6bbWW_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_008-a.root'
T6bbWW_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-036_Figure_008-a.root','orig/CMS-SUS-16-036_Figure_008-a.root','orig/CMS-SUS-16-036_Figure_008-a.root','orig/CMS-SUS-16-036_Figure_008-a.root','orig/CMS-SUS-16-036_Figure_008-a.root','orig/CMS-SUS-16-036_Figure_008-a.root','orig/CMS-SUS-16-036_Figure_008-a.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['ExpLim;1','ExpLim_1Sdn;1','ExpLim_1Sup;1','ObsLim;1','ObsLim_1Sdn;1','ObsLim_1Sup;1','ObsXS_2D;1'],units=[None,None,None,None,None,None,'pb'])


#+++++txName block +++++++++++++++++

T2cc=dataset.addTxName('T2cc')
T2cc.checked=''
T2cc.constraint="[[['jet']],[['jet']]]"
T2cc.condition=None
T2cc.conditionDescription = None
T2cc.source="CMS"

#++++++next mass plane block+++++++++

T2cc_1 = T2cc.addMassPlane(2*[[x,y]])
T2cc_1.figure='Fig. 8-c'
T2cc_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_008-b.png'
T2cc_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_008-b.root'
T2cc_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_008-b.root'
T2cc_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-036/CMS-SUS-16-036_Figure_008-b.root'
T2cc_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-036_Figure_008-c.root','orig/CMS-SUS-16-036_Figure_008-c.root','orig/CMS-SUS-16-036_Figure_008-c.root','orig/CMS-SUS-16-036_Figure_008-c.root','orig/CMS-SUS-16-036_Figure_008-c.root','orig/CMS-SUS-16-036_Figure_008-c.root','orig/CMS-SUS-16-036_Figure_008-c.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['ExpLim;1','ExpLim_1Sdn;1','ExpLim_1Sup;1','ObsLim;1','ObsLim_1Sdn;1','ObsLim_1Sup;1','ObsXS_2D;1'],units=[None,None,None,None,None,None,'pb'])

databaseCreator.create()
