#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse
import types

argparser = argparse.ArgumentParser(description = \
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath', 
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

sys.path.append(os.path.abspath(utilsPath))
from smodels_utils.dataPreparation.inputObjects import TxNameInput, MetaInfoInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.origPlotObjects import x, y

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-PAS-15-002')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
info.sqrts = 13
info.lumi = 2.2
info.prettyName = ''
info.private = False
info.arxiv = ''
info.contact = ''
info.publication = ''
info.comment = 'Only CDS entry https://cds.cern.ch/record/2114817'
info.supersedes = 'None'
#info.supersededBy = 'SUS-15-002'

#+++++++ next txName block ++++++++++++++
T1bbbb = TxNameInput('T1bbbb')
T1bbbb.on.checked = ''
T1bbbb.off.checked = ''
T1bbbb.on.constraint = "[[['b','b']],[['b','b']]]"
#T1bbbb.off.constraint =
T1bbbb.on.conditionDescription = None
#T1bbbb.off.conditionDescription =
T1bbbb.on.condition = None
#T1bbbb.off.condition =
#T1bbbb.branchingRatio =

#+++++++ next mass plane block ++++++++++++++
T1bbbb_1 = T1bbbb.addMassPlane(motherMass = x, lspMass = y)
#----figure----
T1bbbb_1.figure = 'Fig.8-a'
T1bbbb_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/CMS-PAS-SUS-15-002_Figure_008-a.png'
#----limit source----
T1bbbb_1.obsUpperLimit.setSource( 'orig/MassScanT1bbbbSmooth.root', 'root', objectName = 'MassScan2D', index = None )
#T1bbbb_1.expUpperLimit.setSource( 'orig/MassScanT1bbbbSmooth.root', 'root', objectName = '', index = None )
#----exclusion source----
T1bbbb_1.obsExclusion.setSource( 'orig/MassScanT1bbbbSmooth.root', 'root', objectName = 'ObsLim', index = None )
T1bbbb_1.obsExclusionM1.setSource( 'orig/MassScanT1bbbbSmooth.root', 'root', objectName = 'ObsLimSdn', index = None )
T1bbbb_1.obsExclusionP1.setSource( 'orig/MassScanT1bbbbSmooth.root', 'root', objectName = 'ObsLimSup', index = None )
T1bbbb_1.expExclusion.setSource( 'orig/MassScanT1bbbbSmooth.root', 'root', objectName = 'ExpLim', index = None )
T1bbbb_1.expExclusionM1.setSource( 'orig/MassScanT1bbbbSmooth.root', 'root', objectName = 'ExpLimSdn', index = None )
T1bbbb_1.expExclusionP1.setSource( 'orig/MassScanT1bbbbSmooth.root', 'root', objectName = 'ExpLimSup', index = None  )
#----global url settings ----
T1bbbb_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
#----limit url settings ----
T1bbbb_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
T1bbbb_1.obsUpperLimit.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1bbbbSmooth.root'
#T1bbbb_1.expUpperlimit.dataUrl =
#----exclusion url settings ----
T1bbbb_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
T1bbbb_1.obsExclusion.dataUrl =  'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1bbbbSmooth.root'
T1bbbb_1.obsExclusionM1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1bbbbSmooth.root'
T1bbbb_1.obsExclusionP1.dataUrl ='http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1bbbbSmooth.root'
T1bbbb_1.expExclusion.dataUrl ='http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1bbbbSmooth.root'
T1bbbb_1.expExclusionM1.dataUrl ='http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1bbbbSmooth.root'
T1bbbb_1.expExclusionP1.dataUrl ='http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1bbbbSmooth.root'

#+++++++ next txName block ++++++++++++++
T1tttt = TxNameInput('T1tttt')
T1tttt.on.checked = ''
T1tttt.off.checked = ''
T1tttt.on.constraint =  "[[['t','t']],[['t','t']]]"
T1tttt.off.constraint = "[[['b','W','b','W']],[['b','W','b','W']]]"
T1tttt.on.conditionDescription = None
T1tttt.off.conditionDescription = None
T1tttt.on.condition = None
T1tttt.off.condition = None
#T1tttt.branchingRatio =

#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(motherMass = x, lspMass = y)
#----figure----
T1tttt_1.figure = 'Fig.8-b'
T1tttt_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/CMS-PAS-SUS-15-002_Figure_008-b.png'
#----limit source----
T1tttt_1.obsUpperLimit.setSource( 'orig/MassScanT1ttttSmooth.root', 'root', objectName = 'MassScan2D', index = None )
#T1tttt_1.expUpperLimit.setSource( path, type, objectName = None, index = None )
#----exclusion source----
T1tttt_1.obsExclusion.setSource( 'orig/MassScanT1ttttSmooth.root', 'root', objectName = 'ObsLim', index = None )
T1tttt_1.obsExclusionM1.setSource(  'orig/MassScanT1ttttSmooth.root', 'root', objectName = 'ObsLimSdn', index = None )
T1tttt_1.obsExclusionP1.setSource( 'orig/MassScanT1ttttSmooth.root', 'root', objectName = 'ObsLimSup', index = None )
T1tttt_1.expExclusion.setSource( 'orig/MassScanT1ttttSmooth.root', 'root', objectName = 'ExpLim', index = None )
T1tttt_1.expExclusionM1.setSource( 'orig/MassScanT1ttttSmooth.root', 'root', objectName = 'ExpLimSdn', index = None )
T1tttt_1.expExclusionP1.setSource( 'orig/MassScanT1ttttSmooth.root', 'root', objectName = 'ExpLimSup', index = None )
#----global url settings ----
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
#----limit url settings ----
T1tttt_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
T1tttt_1.obsUpperLimit.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1ttttSmooth.root'
#T1tttt_1.expUpperlimit.dataUrl =
#----exclusion url settings ----
T1tttt_1.exclusionDataUrl =  'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1ttttSmooth.root'
T1tttt_1.obsExclusion.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1ttttSmooth.root'
T1tttt_1.obsExclusionM1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1ttttSmooth.root'
T1tttt_1.obsExclusionP1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1ttttSmooth.root'
T1tttt_1.expExclusion.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1ttttSmooth.root'
T1tttt_1.expExclusionM1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1ttttSmooth.root'
T1tttt_1.expExclusionP1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1ttttSmooth.root'

#+++++++ next txName block ++++++++++++++
T1 = TxNameInput('T1')
T1.on.checked =''
T1.off.checked =''
T1.on.constraint = "[[['jet','jet']],[['jet','jet']]]"
#T1.off.constraint =
T1.on.conditionDescription = None
#T1.off.conditionDescription =
T1.on.condition = None
#T1.off.condition =
#T1.branchingRatio =

#+++++++ next mass plane block ++++++++++++++
T1_1 = T1.addMassPlane(motherMass = x, lspMass = y)
#----figure----
T1_1.figure = 'Fig.8-c'
T1_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/CMS-PAS-SUS-15-002_Figure_008-c.png'
#----limit source----
T1_1.obsUpperLimit.setSource( 'orig/MassScanT1qqqqSmooth.root', 'root', objectName = 'MassScan2D', index = None)
#T1_1.expUpperLimit.setSource( path, type, objectName = None, index = None )
#----exclusion source----
T1_1.obsExclusion.setSource( 'orig/MassScanT1qqqqSmooth.root', 'root', objectName = 'ObsLim', index = None  )
T1_1.obsExclusionM1.setSource( 'orig/MassScanT1qqqqSmooth.root', 'root', objectName = 'ObsLimSdn', index = None  )
T1_1.obsExclusionP1.setSource( 'orig/MassScanT1qqqqSmooth.root', 'root', objectName = 'ObsLimSup', index = None  )
T1_1.expExclusion.setSource( 'orig/MassScanT1qqqqSmooth.root', 'root', objectName = 'ExpLim', index = None )
T1_1.expExclusionM1.setSource( 'orig/MassScanT1qqqqSmooth.root', 'root', objectName = 'ExpLimSdn', index = None  )
T1_1.expExclusionP1.setSource( 'orig/MassScanT1qqqqSmooth.root', 'root', objectName = 'ExpLimSup', index = None)
#----global url settings ----
T1_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
#----limit url settings ----
T1_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/index.html'
T1_1.obsUpperLimit.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1Smooth.root'
#T1_1.expUpperlimit.dataUrl =
#----exclusion url settings ----
T1_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1qqqqSmooth.root'
T1_1.obsExclusion.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1qqqqSmooth.root'
T1_1.obsExclusionM1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1qqqqSmooth.root'
T1_1.obsExclusionP1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1qqqqSmooth.root'
T1_1.expExclusion.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1qqqqSmooth.root'
T1_1.expExclusionM1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1qqqqSmooth.root'
T1_1.expExclusionP1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/superseded/SUS-15-002/MassScanT1qqqqSmooth.root'

databaseCreator.create()
