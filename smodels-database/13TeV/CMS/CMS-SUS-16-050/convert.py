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

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-16-050')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/index.html'
info.sqrts = '13.0*TeV'
info.lumi = 35.9
info.prettyName = '0L + top tag'
info.publication = "Phys. Rev. D 97 (2018) 012007, http://dx.doi.org/10.1103/PhysRevD.97.012007"
info.private = False
info.arxiv = 'https://arxiv.org/abs/1710.11188'
info.comment = 'no charm tagging is performed in analyses, therefore we apply the T5tctc result also to up quarks' ## LHCP 2017. root files do not contain limit curves. root files Figure_009-c.root and Figure_009-d.root are swapped!'
info.implementedBy = 'WW'
info.supersedes = "CMS-PAS-SUS-16-050"
info.contact = 'CMS collaboration, cms-phys-conveners-sus@cern.ch'

dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

T2tt = dataset.addTxName('T2tt')
T2tt.checked=''
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription = None
T2tt.condition =None
T2tt.source ='CMS'

T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource( 'upperLimits', 'orig/CMS-SUS-16-050_Figure_008.root', 'root',
		objectName = 'combined_obsLimit_BR100pct' )
T2tt_1.addSource( 'expectedUpperLimits', 'orig/CMS-SUS-16-050_Figure_008.root', 'root',
		objectName = 'combined_expLimit_BR100pct' )
#T2tt_1.removeArea ( [ [ 150,0 ], [200,0 ], [300,100], [260,110] ] )
T2tt_1.addSource( 'obsExclusion', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ObsLim' )
T2tt_1.addSource( 'obsExclusionP1', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ObsLimSup' )
T2tt_1.addSource( 'obsExclusionM1', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ObsLimSdn' )
T2tt_1.addSource( 'expExclusion', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ExpLim' )
T2tt_1.addSource( 'expExclusionP1', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ExpLimSup' )
T2tt_1.addSource( 'expExclusionM1', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ExpLimSdn' )
## fixme add expected
T2tt_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/CMS-SUS-16-050_Figure_008.root"
T2tt_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/CMS-SUS-16-050_Figure_008.png"

T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked =''
T2ttoff.constraint = "[[['b','W']],[['b','W']]]"
T2ttoff.conditionDescription =None
T2ttoff.condition =None
T2ttoff.massConstraint = [['dm <= 169.0'], ['dm <= 169.0']]
T2ttoff.source = "CMS"
T2ttoff.addMassPlane(T2tt_1)

T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked=''
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition =None
T1tttt.source ='CMS'

T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.addSource( 'upperLimits', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root',
		objectName = 'combined_obsLimit_BR100pct' )
T1tttt_1.addSource( 'expectedUpperLimits', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root',
		objectName = 'combined_expLimit_BR100pct' )
T1tttt_1.addSource( 'obsExclusion', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ObsLim' )
T1tttt_1.addSource( 'obsExclusionP1', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ObsLimSup' )
T1tttt_1.addSource( 'obsExclusionM1', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ObsLimSdn' )
T1tttt_1.addSource( 'expExclusion', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ExpLim' )
T1tttt_1.addSource( 'expExclusionP1', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ExpLimSup' )
T1tttt_1.addSource( 'expExclusionM1', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ExpLimSdn' )
T1tttt_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/CMS-SUS-16-050_Figure_009-a.root"
T1tttt_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/CMS-SUS-16-050_Figure_009-a.png"

T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
T1ttttoff.addMassPlane(T1tttt_1)

T5tctc = dataset.addTxName('T5tctc')
T5tctc.checked=''
T5tctc.constraint = "[[['t'],['jet']],[['t'],['jet']]]"
T5tctc.conditionDescription = None
T5tctc.condition =None
T5tctc.source ='CMS'

T5tctc_1 = T5tctc.addMassPlane([[x,y+20,y]]*2)
T5tctc_1.addSource( 'upperLimits', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root',
		objectName = 'combined_obsLimit_BR100pct' )
T5tctc_1.addSource( 'expectedUpperLimits', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root',
		objectName = 'combined_expLimit_BR100pct' )
T5tctc_1.addSource( 'obsExclusion', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ObsLim' )
T5tctc_1.addSource( 'obsExclusionP1', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ObsLimSup' )
T5tctc_1.addSource( 'obsExclusionM1', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ObsLimSdn' )
T5tctc_1.addSource( 'expExclusion', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ExpLim' )
T5tctc_1.addSource( 'expExclusionP1', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ExpLimSup' )
T5tctc_1.addSource( 'expExclusionM1', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ExpLimSdn' )
T5tctc_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/CMS-SUS-16-050_Figure_009-d.root"
T5tctc_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/CMS-SUS-16-050_Figure_009-d.png"

T5tttt = dataset.addTxName('T5tttt')
T5tttt.checked=''
T5tttt.constraint = "[[['t'],['t']],[['t'],['t']]]"
T5tttt.conditionDescription = None
T5tttt.condition =None
T5tttt.source ='CMS'

T5tttt_1 = T5tttt.addMassPlane([[x,y+175,y]]*2)
T5tttt_1.addSource( 'upperLimits', 'orig/CMS-SUS-16-050_Figure_009-c.root', 'root',
		objectName = 'combined_obsLimit_BR100pct' )
T5tttt_1.addSource( 'expectedUpperLimits', 'orig/CMS-SUS-16-050_Figure_009-c.root', 'root',
		objectName = 'combined_expLimit_BR100pct' )
T5tttt_1.addSource( 'obsExclusion', 'orig/CMS-SUS-16-050_Figure_009-c.root', 'root', objectName = 'ObsLim' )
T5tttt_1.addSource( 'obsExclusionP1', 'orig/CMS-SUS-16-050_Figure_009-c.root', 'root', objectName = 'ObsLimSup' )
T5tttt_1.addSource( 'obsExclusionM1', 'orig/CMS-SUS-16-050_Figure_009-c.root', 'root', objectName = 'ObsLimSdn' )
T5tttt_1.addSource( 'expExclusion', 'orig/CMS-SUS-16-050_Figure_009-c.root', 'root', objectName = 'ExpLim' )
T5tttt_1.addSource( 'expExclusionP1', 'orig/CMS-SUS-16-050_Figure_009-c.root', 'root', objectName = 'ExpLimSup' )
T5tttt_1.addSource( 'expExclusionM1', 'orig/CMS-SUS-16-050_Figure_009-c.root', 'root', objectName = 'ExpLimSdn' )
T5tttt_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/CMS-SUS-16-050_Figure_009-c.root"
T5tttt_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/CMS-SUS-16-050_Figure_009-c.png"

databaseCreator.create()
