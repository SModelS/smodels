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
info = MetaInfoInput('CMS-SUS-16-051')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-051/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = '1L stop'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1706.04402'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'JHEP 10 (2017) 019, http://dx.doi.org/10.1007/JHEP10(2017)019'
info.comment = 'Moriond 2017. Fig. 7 with asymmetric decay of stop to b+chargino and top+neutralino not implemented because both, BRs and chargino mass are fixed.'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++txName block +++++++++++++++++

T2tt=dataset.addTxName('T2tt')
T2tt.checked=''
T2tt.constraint="[[['t']],[['t']]]"
T2tt.condition=None
T2tt.conditionDescription = None
T2tt.massConstraint=None
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
T2tt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-051/CMS-SUS-16-051_Figure_005.png'
T2tt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-051/CMS-SUS-16-051_Figure_006.root'
T2tt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-051/CMS-SUS-16-051_Figure_006.root'
T2tt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-051/CMS-SUS-16-051_Figure_006.root'
T2tt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-051_Figure_005.root','orig/CMS-SUS-16-051_Figure_005.root','orig/CMS-SUS-16-051_Figure_005.root','orig/CMS-SUS-16-051_Figure_005.root','orig/CMS-SUS-16-051_Figure_005.root','orig/CMS-SUS-16-051_Figure_005.root','orig/CMS-SUS-16-051_Figure_005.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['gExp;1','gExp1m;1','gExp1p;1','gObs;1','gObs1m;1','gObs1p;1','hObsXsec;1'],units=[None,None,None,None,None,None,'pb'])

T2ttoff.addMassPlane(T2tt_1)

#+++++txName block +++++++++++++++++

T6bbWW=dataset.addTxName('T6bbWW')
T6bbWW.round_to = 6
T6bbWW.checked=''
T6bbWW.constraint="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.condition=None
T6bbWW.conditionDescription = None
T6bbWW.source="CMS"
T6bbWW.massConstraint=None
#T6bbWW.round_to = 6

#++++++next mass plane block+++++++++

T6bbWW_1 = T6bbWW.addMassPlane(2*[[x,0.5*x+0.5*y,y]])
T6bbWW_1.figure='Fig. 6'
T6bbWW_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-051/CMS-SUS-16-051_Figure_006.png.png'
T6bbWW_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-051/CMS-SUS-16-051_Figure_006.root'
T6bbWW_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-051/CMS-SUS-16-051_Figure_006.root'
T6bbWW_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-051/CMS-SUS-16-051_Figure_006.root'
T6bbWW_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-051_Figure_006.root','orig/CMS-SUS-16-051_Figure_006.root','orig/CMS-SUS-16-051_Figure_006.root','orig/CMS-SUS-16-051_Figure_006.root','orig/CMS-SUS-16-051_Figure_006.root','orig/CMS-SUS-16-051_Figure_006.root','orig/CMS-SUS-16-051_Figure_006.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['gExp;1','gExp1m;1','gExp1p;1','gObs;1','gObs1m;1','gObs1p;1','hObsXsec;1'],units=[None,None,None,None,None,None,'pb'])



databaseCreator.create()
