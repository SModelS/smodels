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
argparser.add_argument ('-utilsPath', '-utilsPath', 
help = 'path to the package smodels_utils',\
type = types.StringType)
argparser.add_argument ('-smodelsPath', '-smodelsPath', 
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
info = MetaInfoInput('CMS-SUS-16-049')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'All hadronic stop'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1707.03316'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'J. High Energy Phys. 10 (2017) 005, http://dx.doi.org/10.1007/JHEP10(2017)005'
info.comment = 'Moriond 2017. Paper contains 6 SMS interpretations, 5 of which are implemented here. Not implemented: Fig. 9 with stop -> b chi+1 on one branch and stop -> t chi01 on the other branch, because both the BRs and the chargino mass are fixed.'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++txName block +++++++++++++++++

T2tt=dataset.addTxName('T2tt')
T2tt.checked=''
T2tt.constraint="[[['t']],[['t']]]"
T2tt.condition=None
T2tt.conditionDescription = None
T2tt.massConstraint=[['dm > 169.']]*2  #Use only on-shell region to avoid interpolating in the excluded band
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
T2tt_1.figure='Fig. 7'
T2tt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_007.png'
T2tt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_007.root'
T2tt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_007.root'
T2tt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_007.root'
T2tt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-049_Figure_007.root','orig/CMS-SUS-16-049_Figure_007.root','orig/CMS-SUS-16-049_Figure_007.root','orig/CMS-SUS-16-049_Figure_007.root','orig/CMS-SUS-16-049_Figure_007.root','orig/CMS-SUS-16-049_Figure_007.root','orig/CMS-SUS-16-049_Figure_007.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['graph_smoothed_Exp;1','graph_smoothed_ExpM;1','graph_smoothed_ExpP;1','graph_smoothed_Obs;1','graph_smoothed_ObsM;1','graph_smoothed_ObsP;1','hXsec_obs_corrRemoved;1'],units=[None,None,None,None,None,None,'pb'])
T2ttoff.addMassPlane(T2tt_1)
#+++++txName block +++++++++++++++++

T2cc=dataset.addTxName('T2cc')
T2cc.checked=''
T2cc.constraint="[[['jet']],[['jet']]]"
T2cc.condition=None
T2cc.conditionDescription = None
T2cc.massConstraint=None
T2cc.source="CMS"

#++++++next mass plane block+++++++++

T2cc_1 = T2cc.addMassPlane(2*[[x,x-y]])
T2cc_1.figure='Fig. 12'
T2cc_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_012.png'
T2cc_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_012.root'
T2cc_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_012.root'
T2cc_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_012.root'
T2cc_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-049_Figure_012.root','orig/CMS-SUS-16-049_Figure_012.root','orig/CMS-SUS-16-049_Figure_012.root','orig/CMS-SUS-16-049_Figure_012.root','orig/CMS-SUS-16-049_Figure_012.root','orig/CMS-SUS-16-049_Figure_012.root','orig/CMS-SUS-16-049_Figure_012.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['graph_smoothed_Exp;1','graph_smoothed_ExpM;1','graph_smoothed_ExpP;1','graph_smoothed_Obs;1','graph_smoothed_ObsM;1','graph_smoothed_ObsP_1;1','hXsec_obs_corr;1'],units=[None,None,None,None,None,None,'pb'])

#+++++txName block +++++++++++++++++

T6bbWW=dataset.addTxName('T6bbWW')
T6bbWW.checked=''
T6bbWW.constraint="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.condition=None
T6bbWW.conditionDescription = None
T6bbWW.massConstraint=None
T6bbWW.source="CMS"

#++++++next mass plane block+++++++++

T6bbWW_1 = T6bbWW.addMassPlane(2*[[x,0.5*x+0.5*y,y]])
T6bbWW_1.figure='Fig. 8'
T6bbWW_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_008.png'
T6bbWW_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_008.root'
T6bbWW_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_008.root'
T6bbWW_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_008.root'
T6bbWW_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-049_Figure_008.root','orig/CMS-SUS-16-049_Figure_008.root','orig/CMS-SUS-16-049_Figure_008.root','orig/CMS-SUS-16-049_Figure_008.root','orig/CMS-SUS-16-049_Figure_008.root','orig/CMS-SUS-16-049_Figure_008.root','orig/CMS-SUS-16-049_Figure_008.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['graph_smoothed_Exp;1','graph_smoothed_ExpM;1','graph_smoothed_ExpP;1','graph_smoothed_Obs;1','graph_smoothed_ObsM;1','graph_smoothed_ObsP;1','hXsec_obs_corr;1'],units=[None,None,None,None,None,None,'pb'])

#+++++txName block +++++++++++++++++

T2ttC=dataset.addTxName('T2ttC')
T2ttC.checked=''
T2ttC.constraint="2.23*([[['b','jet','jet']],[['b','jet','jet']]])"
T2ttC.condition=None
T2ttC.conditionDescription = None
T2ttC.massConstraint=[['dm<80.0'],['dm<80.0']]
T2ttC.source="CMS"

#++++++next mass plane block+++++++++


T2ttC_1 = T2ttC.addMassPlane(2*[[x,x-y]])
T2ttC_1.figure='Fig. 10'
T2ttC_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_010.png'
T2ttC_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_010.root'
T2ttC_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_010.root'
T2ttC_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_010.root'
T2ttC_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-049_Figure_010.root','orig/CMS-SUS-16-049_Figure_010.root','orig/CMS-SUS-16-049_Figure_010.root','orig/CMS-SUS-16-049_Figure_010.root','orig/CMS-SUS-16-049_Figure_010.root','orig/CMS-SUS-16-049_Figure_010.root','orig/CMS-SUS-16-049_Figure_010.root'],
                   dataFormats=['root','root','root','root','root','root','root'],objectNames=['graph_smoothed_Exp;1','graph_smoothed_ExpM;1','graph_smoothed_ExpP;1','graph_smoothed_Obs;1','graph_smoothed_ObsM;1',
'graph_smoothed_ObsP;1','hXsec_obs_corr;1'],units=[None,None,None,None,None,None,'pb'])


#+++++txName block +++++++++++++++++

T6bbWWoff=dataset.addTxName('T6bbWWoff')
T6bbWWoff.checked=''
T6bbWWoff.constraint="2.23*([[['b'],['jet','jet']],[['b'],['jet','jet']]])"
T6bbWWoff.condition = None
T6bbWWoff.conditionDescription = None
T6bbWWoff.massConstraint=[['dm>0','dm<76.0'],['dm>0','dm<76.0']]
T6bbWWoff.source="CMS"

#++++++next mass plane block+++++++++


T6bbWWoff_1 = T6bbWWoff.addMassPlane(2*[[x,x-0.5*y,x-y]])
T6bbWWoff_1.figure='Fig. 11'
T6bbWWoff_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_011.png'
T6bbWWoff_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_011.root'
T6bbWWoff_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_011.root'
T6bbWWoff_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-049/CMS-SUS-16-049_Figure_011.root'
T6bbWWoff_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-049_Figure_011.root','orig/CMS-SUS-16-049_Figure_011.root','orig/CMS-SUS-16-049_Figure_011.root','orig/CMS-SUS-16-049_Figure_011.root','orig/CMS-SUS-16-049_Figure_011.root','orig/CMS-SUS-16-049_Figure_011.root','orig/CMS-SUS-16-049_Figure_011.root'],
                   dataFormats=['root','root','root','root','root','root','root'],objectNames=['graph_smoothed_Exp;1','graph_smoothed_ExpM;1','graph_smoothed_ExpP;1','graph_smoothed_Obs;1','graph_smoothed_ObsM;1',
'graph_smoothed_ObsP;1','hXsec_obs_corr;1'],units=[None,None,None,None,None,None,'pb'])

databaseCreator.create()
