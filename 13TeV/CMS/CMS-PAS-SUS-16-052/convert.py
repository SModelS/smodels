#!/usr/bin/env python3

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
type = str)
argparser.add_argument ('-smodelsPath', '--smodelsPath',
help = 'path to the package smodels_utils',\
type = str)
argparser.add_argument ('-no', '--noUpdate',
help = 'do not update the lastUpdate field.',\
action= "store_true" )
args = argparser.parse_args()

if args.noUpdate:
    os.environ["SMODELS_NOUPDATE"]="1"

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
from smodels_utils.dataPreparation.datasetCreation import DatasetsFromLatex
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z

# DataSetInput.ntoys = 200

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-PAS-SUS-16-052')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/index.html'
info.sqrts = '13.0*TeV'
info.lumi = 35.9
info.prettyName = 'soft lepton, <= 2 jets'
info.private = False
# info.arxiv = 'https://arxiv.org/abs/....'
info.comment = 'EPS 2017. https://cds.cern.ch/record/2273394, http://inspirehep.net/record/1609006'
info.implementedBy = 'WW'
info.contact = 'CMS'

dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2bbWWoff = dataset.addTxName('T2bbWWoff')
T2bbWWoff.checked =''
T2bbWWoff.constraint = "3.47*([[['b','l','nu']],[['b','jet','jet']]])"
#T2bbWWoff.constraint = "3.47*([[['b','mu','nu']],[['b','jet','jet']]]+[[['b','e','nu']],[['b','jet','jet']]])"
T2bbWWoff.conditionDescription = None
# T2bbWWoff.massConstraint = [['dm <= 76.']]*2
T2bbWWoff.condition =None
T2bbWWoff.source ='CMS'                                                                #+++++++ next mass plane block ++++++++++++++
T2bbWWoff_1 = T2bbWWoff.addMassPlane([[x,x-y]]*2)
#----exclusion source----
T2bbWWoff_1.addSource( 'upperLimits', 'orig/CMS-PAS-SUS-16-052_Figure_004.root', 'canvas', objectName = 'cCONT_', index=2 )
T2bbWWoff_1.addSource( 'obsExclusion', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                    'canvas', objectName = "cCONT_", index=6 )
T2bbWWoff_1.addSource( 'obsExclusionP1', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                    'canvas', objectName = "cCONT_", index=7 )
T2bbWWoff_1.addSource( 'obsExclusionM1', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                    'canvas', objectName = "cCONT_", index=8 )
T2bbWWoff_1.addSource( 'expExclusion', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                    'canvas', objectName = "cCONT_", index=3 )
T2bbWWoff_1.addSource( 'expExclusionP1', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                    'canvas', objectName = "cCONT_", index=4 )
T2bbWWoff_1.addSource( 'expExclusionM1', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                    'canvas', objectName = "cCONT_", index=5 )
T2bbWWoff_1.dataUrl = "https://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/CMS-PAS-SUS-16-052_Figure_004.root"
T2bbWWoff_1.figureUrl = "https://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/CMS-PAS-SUS-16-052_Figure_004.png"

#+++++++ next txName block ++++++++++++++
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.checked =''
T6bbWWoff.constraint = "3.47*([[['b'],['mu','nu']],[['b'],['jet','jet']]]+[[['b'],['e','nu']],[['b'],['jet','jet']]])"
T6bbWWoff.conditionDescription = None
#T6bbWWoff.massConstraint = [['dm <= 76.']]*2
T6bbWWoff.condition =None
T6bbWWoff.source ='CMS'                                                                #+++++++ next mass plane block ++++++++++++++
T6bbWWoff_1 = T6bbWWoff.addMassPlane([[x,x-.5*y,x-y]]*2)
#----exclusion source----
T6bbWWoff_1.addSource( 'upperLimits', 'orig/CMS-PAS-SUS-16-052_Figure_005.root', 'canvas', objectName = 'cCONT_', index=2 )
T6bbWWoff_1.addSource( 'obsExclusion', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                    'canvas', objectName = "cCONT_", index=6 )
T6bbWWoff_1.addSource( 'obsExclusionP1', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                    'canvas', objectName = "cCONT_", index=7 )
T6bbWWoff_1.addSource( 'obsExclusionM1', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                    'canvas', objectName = "cCONT_", index=8 )
T6bbWWoff_1.addSource( 'expExclusion', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                    'canvas', objectName = "cCONT_", index=3 )
T6bbWWoff_1.addSource( 'expExclusionP1', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                    'canvas', objectName = "cCONT_", index=4 )
T6bbWWoff_1.addSource( 'expExclusionM1', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                    'canvas', objectName = "cCONT_", index=5 )
T6bbWWoff_1.dataUrl = "https://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/CMS-PAS-SUS-16-052_Figure_005.root"
T6bbWWoff_1.figureUrl = "https://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/CMS-PAS-SUS-16-052_Figure_005.png"

databaseCreator.create()
