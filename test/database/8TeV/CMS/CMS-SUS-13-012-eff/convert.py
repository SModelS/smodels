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
info = MetaInfoInput('CMS-SUS-13-012')
#info.comment =
info.sqrts = '8.0'
info.private = False
info.lumi = '19.5'
info.publication = 'JHEP06(2014)055'
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13012'
#info.superseded_by =
info.arxiv = 'http://arxiv.org/abs/1402.4770'
info.contact = 'cms-pag-conveners-sus@NOSPAMcernSPAMNOT.ch'
info.prettyName = ''
info.implementedBy = 'Federico A.'
# info.supersedes =


#+++++++ T2 Block  ++++++++++++++
T2 = TxNameInput('T2')
#T2.on.checked =" "
#T1bbbb.off.checked =
T2.on.constraint ="[[['jet']],[['jet']]]"
#T1bbbb.off.constraint =
T2.on.conditionDescription ="None"
#T1bbbb.off.condition =
T2.on.condition ="None"
T2 = T2.addMassPlane(motherMass = x, lspMass = y)

#----exclusion source----
T2.obsExclusion.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_obsExclOneTimesProspino", index = None )
T2.obsExclusionM1.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_obsExclMinusSysErrProspino", index = None )
T2.obsExclusionP1.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_obsExclPlusSysErrProspino", index = None )
T2.expExclusion.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_expExclOneTimesProspino", index = None )
T2.expExclusionM1.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_expExclMinusOneSigmaProspino", index = None )
T2.expExclusionP1.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_expExclPlusOneSigmaProspino", index = None )

T2.figure = "Fig_7a"
T2.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7a.pdf"

# +++++++ T5WZ block ++++++++++++++
T5WZ = TxNameInput('T5WZ')
#T5WZ.on.checked =" "
#T5WZ.off.checked =
T5WZ.on.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['Z']]] "
#T5WZ.off.constraint =
T5WZ.on.conditionDescription ="None"
#T5WZ.off.condition =
T5WZ.on.condition ="None"
T5WZ = T5WZ.addMassPlane(motherMass = x, interMass0 = 0.5*(x+y), lspMass = y)

#----exclusion source----
T5WZ.obsExclusion.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "combined_obsExclOneTimesProspino", index = None )
T5WZ.obsExclusionM1.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "combined_obsExclMinusSysErrProspino", index = None )
T5WZ.obsExclusionP1.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "combined_obsExclPlusSysErrProspino", index = None )
T5WZ.expExclusion.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "combined_expExclOneTimesProspino", index = None )
T5WZ.expExclusionM1.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "combined_expExclMinusOneSigmaProspino", index = None )
T5WZ.expExclusionP1.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "combined_expExclPlusOneSigmaProspino", index = None )

T5WZ.figure = "Fig_7d"
T5WZ.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7d.pdf"

#+++++++ T1tttt Block  ++++++++++++
T1tttt = TxNameInput('T1tttt')
T1tttt.on.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.off.constraint = "[[['b','W','b','W']],[['b','W','b','W']]]"
T1tttt.on.conditionDescription ="None"
T1tttt.off.conditionDescription ="None"
T1tttt.off.condition = "None"
T1tttt.on.condition ="None"
T1tttt = T1tttt.addMassPlane(motherMass = x, lspMass = y)

T1tttt.obsExclusion.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_obsExclOneTimesProspino", index = None )
T1tttt.obsExclusionM1.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_obsExclMinusSysErrProspino", index = None )
T1tttt.obsExclusionP1.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_obsExclPlusSysErrProspino", index = None )
T1tttt.expExclusion.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_expExclOneTimesProspino", index = None )
T1tttt.expExclusionM1.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_expExclMinusOneSigmaProspino", index = None )
T1tttt.expExclusionP1.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_expExclPlusOneSigmaProspino", index = None )

T1tttt.figure = "Fig_7c"
T1tttt.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7c.pdf"

# +++++++ T1 block ++++++++++++++
T1 = TxNameInput('T1')
#T1.on.checked =" "
#T1bbbb.off.checked =
T1.on.constraint ="[[['jet','jet']],[['jet','jet']]]"
#T1bbbb.off.constraint =
T1.on.conditionDescription ="None"
#T1bbbb.off.condition =
T1.on.condition ="None"
#T1bbbb.off.fuzzycondition =
T1 = T1.addMassPlane(motherMass = x, lspMass = y)




#----limit source----

## T1.expUpperLimit.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_expLimit", index = None )
## T1.obsUpperLimit.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_obsLimit", index = None )
#----exclusion source----
T1.obsExclusion.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_obsExclOneTimesProspino", index = None )
T1.obsExclusionM1.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_obsExclMinusSysErrProspino", index = None )
T1.obsExclusionP1.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_obsExclPlusSysErrProspino", index = None )
T1.expExclusion.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_expExclOneTimesProspino", index = None )
T1.expExclusionM1.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_expExclMinusOneSigmaProspino", index = None )
T1.expExclusionP1.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_expExclPlusOneSigmaProspino", index = None )

#----global url settings ----
#T1.dataUrl =
#T1.histoDataUrl =
#T1.exclusionDataUrl =
#----figure----
T1.figure = "Fig_7b"
T1.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7b.pdf"
#----limit url settings ----
T1.efficiencyMap.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/c_AccEffMap_T1qqqq.tar"
## T1.expUpperLimit.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/SUS13012_XsecLimits_T1qqqq.root"
#----exclusion url settings ----
T1.obsExclusion.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/SUS13012_XsecLimits_T1qqqq.root"
#T1bbbb.obsExclusionM1.dataUrl =
#T1bbbb.obsExclusionP1.dataUrl =
#T1bbbb.expExclusion.dataUrl =
#T1bbbb.expExclusionM1.dataUrl =
#T1bbbb.expExclusionP1.dataUrl =



#EFFICIENCIES
T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_200MHT300", index = None, dataset="3NJet6_500HT800_200MHT300" )
T2.efficiencyMap.setStatistics ( observedN=6159, expectedBG=6088, bgError=665 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_200MHT300", index = None, dataset="3NJet6_500HT800_200MHT300" )
T1tttt.efficiencyMap.setStatistics ( observedN=6159, expectedBG=6088, bgError=665 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_200MHT300", index = None, dataset="3NJet6_500HT800_200MHT300" )
T1.efficiencyMap.setStatistics ( observedN=6159, expectedBG=6088, bgError=665 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_200MHT300", index = None, dataset="3NJet6_500HT800_200MHT300" )
T5WZ.efficiencyMap.setStatistics ( observedN=6159, expectedBG=6088, bgError=665 )

databaseCreator.create()





T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_300MHT450", index = None, dataset="3NJet6_500HT800_300MHT450" )
T2.efficiencyMap.setStatistics ( observedN=2305, expectedBG=2278, bgError=266 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_300MHT450", index = None, dataset="3NJet6_500HT800_300MHT450" )
T1tttt.efficiencyMap.setStatistics ( observedN=2305, expectedBG=2278, bgError=266 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_300MHT450", index = None, dataset="3NJet6_500HT800_300MHT450" )
T1.efficiencyMap.setStatistics ( observedN=2305, expectedBG=2278, bgError=266 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_300MHT450", index = None, dataset="3NJet6_500HT800_300MHT450" )
T5WZ.efficiencyMap.setStatistics ( observedN=2305, expectedBG=2278, bgError=266 )

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_450MHT600", index = None, dataset="3NJet6_500HT800_450MHT600" )
T2.efficiencyMap.setStatistics ( observedN=454, expectedBG=418, bgError=66 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_450MHT600", index = None, dataset="3NJet6_500HT800_450MHT600" )
T1tttt.efficiencyMap.setStatistics ( observedN=454, expectedBG=418, bgError=66 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_450MHT600", index = None, dataset="3NJet6_500HT800_450MHT600" )
T1.efficiencyMap.setStatistics ( observedN=454, expectedBG=418, bgError=66 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_450MHT600", index = None, dataset="3NJet6_500HT800_450MHT600" )
T5WZ.efficiencyMap.setStatistics ( observedN=454, expectedBG=418, bgError=66 )

databaseCreator.create( True )



T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_600MHTinf", index = None, dataset="3NJet6_500HT800_600MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=62, expectedBG=57.4, bgError=11.2 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_600MHTinf", index = None, dataset="3NJet6_500HT800_600MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=62, expectedBG=57.4, bgError=11.2 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_600MHTinf", index = None, dataset="3NJet6_500HT800_600MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=62, expectedBG=57.4, bgError=11.2 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_500HT800_600MHTinf", index = None, dataset="3NJet6_500HT800_600MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=62, expectedBG=57.4, bgError=11.2 )

databaseCreator.create( True )








T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_200MHT300", index = None, dataset="3NJet6_800HT1000_200MHT300" )
T2.efficiencyMap.setStatistics ( observedN=808, expectedBG=777, bgError=107 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_200MHT300", index = None, dataset="3NJet6_800HT1000_200MHT300" )
T1tttt.efficiencyMap.setStatistics ( observedN=808, expectedBG=777, bgError=107 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_200MHT300", index = None, dataset="3NJet6_800HT1000_200MHT300" )
T1.efficiencyMap.setStatistics ( observedN=808, expectedBG=777, bgError=107 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_200MHT300", index = None, dataset="3NJet6_800HT1000_200MHT300" )
T5WZ.efficiencyMap.setStatistics ( observedN=808, expectedBG=777, bgError=107 )

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_300MHT450", index = None, dataset="3NJet6_800HT1000_300MHT450" )
T2.efficiencyMap.setStatistics ( observedN=305, expectedBG=330, bgError=40 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_300MHT450", index = None, dataset="3NJet6_800HT1000_300MHT450" )
T1tttt.efficiencyMap.setStatistics ( observedN=305, expectedBG=330, bgError=40 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_300MHT450", index = None, dataset="3NJet6_800HT1000_300MHT450" )
T1.efficiencyMap.setStatistics ( observedN=305, expectedBG=330, bgError=40 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_300MHT450", index = None, dataset="3NJet6_800HT1000_300MHT450" )
T5WZ.efficiencyMap.setStatistics ( observedN=305, expectedBG=330, bgError=40 )

databaseCreator.create( True )


T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_450MHT600", index = None, dataset="3NJet6_800HT1000_450MHT600" )
T2.efficiencyMap.setStatistics ( observedN=124, expectedBG=108, bgError=15 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_450MHT600", index = None, dataset="3NJet6_800HT1000_450MHT600" )
T1tttt.efficiencyMap.setStatistics ( observedN=124, expectedBG=108, bgError=15 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_450MHT600", index = None, dataset="3NJet6_800HT1000_450MHT600" )
T1.efficiencyMap.setStatistics ( observedN=124, expectedBG=108, bgError=15 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_450MHT600", index = None, dataset="3NJet6_800HT1000_450MHT600" )
T5WZ.efficiencyMap.setStatistics ( observedN=124, expectedBG=108, bgError=15 )

databaseCreator.create( True )



T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_600MHTinf", index = None, dataset="3NJet6_800HT1000_600MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=52, expectedBG=54.8, bgError=9.7 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_450MHT600", index = None, dataset="3NJet6_800HT1000_450MHT600" )
T1tttt.efficiencyMap.setStatistics ( observedN=52, expectedBG=54.8, bgError=9.7 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_600MHTinf", index = None, dataset="3NJet6_800HT1000_600MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=52, expectedBG=54.8, bgError=9.7 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_800HT1000_600MHTinf", index = None, dataset="3NJet6_800HT1000_600MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=52, expectedBG=54.8, bgError=9.7 )

databaseCreator.create( True )



T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_200MHT300", index = None, dataset="3NJet6_1000HT1250_200MHT300" )
T2.efficiencyMap.setStatistics ( observedN=335, expectedBG=305, bgError=41 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_200MHT300", index = None, dataset="3NJet6_1000HT1250_200MHT300" )
T1tttt.efficiencyMap.setStatistics ( observedN=335, expectedBG=305, bgError=41 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_200MHT300", index = None, dataset="3NJet6_1000HT1250_200MHT300" )
T1.efficiencyMap.setStatistics ( observedN=335, expectedBG=305, bgError=41 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_200MHT300", index = None, dataset="3NJet6_1000HT1250_200MHT300" )
T5WZ.efficiencyMap.setStatistics ( observedN=335, expectedBG=305, bgError=41 )

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_300MHT450", index = None, dataset="3NJet6_1000HT1250_300MHT450" )
T2.efficiencyMap.setStatistics ( observedN=129, expectedBG=137, bgError=20 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_300MHT450", index = None, dataset="3NJet6_1000HT1250_300MHT450" )
T1tttt.efficiencyMap.setStatistics ( observedN=129, expectedBG=137, bgError=20 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_300MHT450", index = None, dataset="3NJet6_1000HT1250_300MHT450" )
T1.efficiencyMap.setStatistics ( observedN=129, expectedBG=137, bgError=20 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_300MHT450", index = None, dataset="3NJet6_1000HT1250_300MHT450" )
T5WZ.efficiencyMap.setStatistics ( observedN=129, expectedBG=137, bgError=20 )

databaseCreator.create( True )



T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_450MHT600", index = None, dataset="3NJet6_1000HT1250_450MHT600" )
T2.efficiencyMap.setStatistics ( observedN=34, expectedBG=32.3, bgError=34 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_450MHT600", index = None, dataset="3NJet6_1000HT1250_450MHT600" )
T1tttt.efficiencyMap.setStatistics ( observedN=34, expectedBG=32.3, bgError=34 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_450MHT600", index = None, dataset="3NJet6_1000HT1250_450MHT600" )
T1.efficiencyMap.setStatistics ( observedN=34, expectedBG=32.3, bgError=34 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_450MHT600", index = None, dataset="3NJet6_1000HT1250_450MHT600" )
T5WZ.efficiencyMap.setStatistics ( observedN=34, expectedBG=32.3, bgError=34 )

databaseCreator.create( True )






T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_600MHTinf", index = None, dataset="3NJet6_1000HT1250_600MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=32, expectedBG=22.8, bgError=5.2 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_600MHTinf", index = None, dataset="3NJet6_1000HT1250_600MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=32, expectedBG=22.8, bgError=5.2 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_600MHTinf", index = None, dataset="3NJet6_1000HT1250_600MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=32, expectedBG=22.8, bgError=5.2 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_1000HT1250_600MHTinf", index = None, dataset="3NJet6_1000HT1250_600MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=32, expectedBG=22.8, bgError=5.2 )

databaseCreator.create( True )






T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_200MHT300", index = None, dataset="3NJet6_1250HT1500_200MHT300" )
T2.efficiencyMap.setStatistics ( observedN=98, expectedBG=109, bgError=18 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_200MHT300", index = None, dataset="3NJet6_1250HT1500_200MHT300" )
T1tttt.efficiencyMap.setStatistics ( observedN=98, expectedBG=109, bgError=18 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_200MHT300", index = None, dataset="3NJet6_1250HT1500_200MHT300" )
T1.efficiencyMap.setStatistics ( observedN=98, expectedBG=109, bgError=18 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_200MHT300", index = None, dataset="3NJet6_1250HT1500_200MHT300" )
T5WZ.efficiencyMap.setStatistics ( observedN=98, expectedBG=109, bgError=18 )

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_300MHT450", index = None, dataset="3NJet6_1250HT1500_300MHT450" )
T2.efficiencyMap.setStatistics ( observedN=38, expectedBG=42.8, bgError=9.5 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_300MHT450", index = None, dataset="3NJet6_1250HT1500_300MHT450" )
T1tttt.efficiencyMap.setStatistics ( observedN=38, expectedBG=42.8, bgError=9.5 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_300MHT450", index = None, dataset="3NJet6_1250HT1500_300MHT450" )
T1.efficiencyMap.setStatistics ( observedN=38, expectedBG=42.8, bgError=9.5 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_300MHT450", index = None, dataset="3NJet6_1250HT1500_300MHT450" )
T5WZ.efficiencyMap.setStatistics ( observedN=38, expectedBG=42.8, bgError=9.5 )

databaseCreator.create( True )







T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_450MHTinf", index = None, dataset="3NJet6_1250HT1500_450MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=23, expectedBG=17.6, bgError=4.1 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_450MHTinf", index = None, dataset="3NJet6_1250HT1500_450MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=23, expectedBG=17.6, bgError=4.1 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_450MHTinf", index = None, dataset="3NJet6_1250HT1500_450MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=23, expectedBG=17.6, bgError=4.1 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_450MHTinf", index = None, dataset="3NJet6_1250HT1500_450MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=23, expectedBG=17.6, bgError=4.1 )

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_1500HTinf_200MHT300", index = None, dataset="3NJet6_1500HTinf_200MHT300" )
T2.efficiencyMap.setStatistics ( observedN=94, expectedBG=86, bgError=17 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_1500HTinf_200MHT300", index = None, dataset="3NJet6_1500HTinf_200MHT300" )
T1tttt.efficiencyMap.setStatistics ( observedN=94, expectedBG=86, bgError=17 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_1500HTinf_200MHT300", index = None, dataset="3NJet6_1500HTinf_200MHT300" )
T1.efficiencyMap.setStatistics ( observedN=94, expectedBG=86, bgError=17 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_1500HTinf_200MHT300", index = None, dataset="3NJet6_1500HTinf_200MHT300" )
T5WZ.efficiencyMap.setStatistics ( observedN=94, expectedBG=86, bgError=17 )

databaseCreator.create( True )





T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_3NJet6_1500HTinf_300MHTinf", index = None, dataset="3NJet6_1500HTinf_300MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=39, expectedBG=29.7, bgError=5.8 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_3NJet6_1500HTinf_300MHTinf", index = None, dataset="3NJet6_1500HTinf_300MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=39, expectedBG=29.7, bgError=5.8 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_3NJet6_1500HTinf_300MHTinf", index = None, dataset="3NJet6_1500HTinf_300MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=39, expectedBG=29.7, bgError=5.8 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_1500HTinf_300MHTinf", index = None, dataset="3NJet6_1500HTinf_300MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=39, expectedBG=29.7, bgError=5.8 )

databaseCreator.create( True )





T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_200MHT300", index = None, dataset="6NJet8_500HT800_200MHT300" )
T2.efficiencyMap.setStatistics ( observedN=266, expectedBG=290, bgError=65 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_200MHT300", index = None, dataset="6NJet8_500HT800_200MHT300" )
T1tttt.efficiencyMap.setStatistics ( observedN=266, expectedBG=290, bgError=65 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_200MHT300", index = None, dataset="6NJet8_500HT800_200MHT300" )
T1.efficiencyMap.setStatistics ( observedN=266, expectedBG=290, bgError=65 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_200MHT300", index = None, dataset="6NJet8_500HT800_200MHT300" )
T5WZ.efficiencyMap.setStatistics ( observedN=266, expectedBG=290, bgError=65 )

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_300MHT450", index = None, dataset="6NJet8_500HT800_300MHT450" )
T2.efficiencyMap.setStatistics ( observedN=62, expectedBG=52, bgError=12 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_300MHT450", index = None, dataset="6NJet8_500HT800_300MHT450" )
T1tttt.efficiencyMap.setStatistics ( observedN=62, expectedBG=52, bgError=12 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_300MHT450", index = None, dataset="6NJet8_500HT800_300MHT450" )
T1.efficiencyMap.setStatistics ( observedN=62, expectedBG=52, bgError=12 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_300MHT450", index = None, dataset="6NJet8_500HT800_300MHT450" )
T5WZ.efficiencyMap.setStatistics ( observedN=62, expectedBG=52, bgError=12 )

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_450MHTinf", index = None, dataset="6NJet8_500HT800_450MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=9, expectedBG=0.8, bgError=0.6 ) ### check error
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_450MHTinf", index = None, dataset="6NJet8_500HT800_450MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=9, expectedBG=0.8, bgError=0.6 ) ### check error
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_450MHTinf", index = None, dataset="6NJet8_500HT800_450MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=9, expectedBG=0.8, bgError=0.6 ) ### check error
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_500HT800_450MHTinf", index = None, dataset="6NJet8_500HT800_450MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=9, expectedBG=0.8, bgError=0.6 ) ### check error

databaseCreator.create( True )






T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_200MHT300", index = None, dataset="6NJet8_800HT1000_200MHT300" )
T2.efficiencyMap.setStatistics ( observedN=111, expectedBG=124, bgError=29 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_200MHT300", index = None, dataset="6NJet8_800HT1000_200MHT300" )
T1tttt.efficiencyMap.setStatistics ( observedN=111, expectedBG=124, bgError=29 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_200MHT300", index = None, dataset="6NJet8_800HT1000_200MHT300" )
T1.efficiencyMap.setStatistics ( observedN=111, expectedBG=124, bgError=29 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_200MHT300", index = None, dataset="6NJet8_800HT1000_200MHT300" )
T5WZ.efficiencyMap.setStatistics ( observedN=111, expectedBG=124, bgError=29 )

databaseCreator.create( True )





T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_300MHT450", index = None, dataset="6NJet8_800HT1000_300MHT450" )
T2.efficiencyMap.setStatistics ( observedN=35, expectedBG=28.6, bgError=6.9 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_300MHT450", index = None, dataset="6NJet8_800HT1000_300MHT450" )
T1tttt.efficiencyMap.setStatistics ( observedN=35, expectedBG=28.6, bgError=6.9 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_300MHT450", index = None, dataset="6NJet8_800HT1000_300MHT450" )
T1.efficiencyMap.setStatistics ( observedN=35, expectedBG=28.6, bgError=6.9 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_300MHT450", index = None, dataset="6NJet8_800HT1000_300MHT450" )
T5WZ.efficiencyMap.setStatistics ( observedN=35, expectedBG=28.6, bgError=6.9 )

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_450MHTinf", index = None, dataset="6NJet8_800HT1000_450MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=4, expectedBG=6.0, bgError=2.8 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_450MHTinf", index = None, dataset="6NJet8_800HT1000_450MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=4, expectedBG=6.0, bgError=2.8 )
T1.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_450MHTinf", index = None, dataset="6NJet8_800HT1000_450MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=4, expectedBG=6.0, bgError=2.8 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_800HT1000_450MHTinf", index = None, dataset="6NJet8_800HT1000_450MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=4, expectedBG=6.0, bgError=2.8 )

databaseCreator.create( True )






T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_200MHT300", index = None, dataset="6NJet8_1000HT1250_200MHT300" )
T2.efficiencyMap.setStatistics ( observedN=67, expectedBG=70, bgError=67 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_200MHT300", index = None, dataset="6NJet8_1000HT1250_200MHT300" )
T1tttt.efficiencyMap.setStatistics ( observedN=67, expectedBG=70, bgError=67 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_200MHT300", index = None, dataset="6NJet8_1000HT1250_200MHT300" )
T1.efficiencyMap.setStatistics ( observedN=67, expectedBG=70, bgError=67 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_200MHT300", index = None, dataset="6NJet8_1000HT1250_200MHT300" )
T5WZ.efficiencyMap.setStatistics ( observedN=67, expectedBG=70, bgError=67 )

databaseCreator.create( True )





T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_300MHT450", index = None, dataset="6NJet8_1000HT1250_300MHT450" )
T2.efficiencyMap.setStatistics ( observedN=20, expectedBG=21.6, bgError=5.8 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_300MHT450", index = None, dataset="6NJet8_1000HT1250_300MHT450" )
T1tttt.efficiencyMap.setStatistics ( observedN=20, expectedBG=21.6, bgError=5.8 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_300MHT450", index = None, dataset="6NJet8_1000HT1250_300MHT450" )
T1.efficiencyMap.setStatistics ( observedN=20, expectedBG=21.6, bgError=5.8 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_300MHT450", index = None, dataset="6NJet8_1000HT1250_300MHT450" )
T5WZ.efficiencyMap.setStatistics ( observedN=20, expectedBG=21.6, bgError=5.8 )

databaseCreator.create( True )







T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_450MHTinf", index = None, dataset="6NJet8_1000HT1250_450MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=4, expectedBG=2.2, bgError=1.1 ) ### check error
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_450MHTinf", index = None, dataset="6NJet8_1000HT1250_450MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=4, expectedBG=2.2, bgError=1.1 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_450MHTinf", index = None, dataset="6NJet8_1000HT1250_450MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=4, expectedBG=2.2, bgError=1.1 ) ### check error
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_1000HT1250_450MHTinf", index = None, dataset="6NJet8_1000HT1250_450MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=4, expectedBG=2.2, bgError=1.1 ) ### check error

databaseCreator.create( True )






T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_200MHT300", index = None, dataset="6NJet8_1250HT1500_200MHT300" )
T2.efficiencyMap.setStatistics ( observedN=24, expectedBG=28.0, bgError=8.2 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_200MHT300", index = None, dataset="6NJet8_1250HT1500_200MHT300" )
T1tttt.efficiencyMap.setStatistics ( observedN=24, expectedBG=28.0, bgError=8.2 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_200MHT300", index = None, dataset="6NJet8_1250HT1500_200MHT300" )
T1.efficiencyMap.setStatistics ( observedN=24, expectedBG=28.0, bgError=8.2 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_200MHT300", index = None, dataset="6NJet8_1250HT1500_200MHT300" )
T5WZ.efficiencyMap.setStatistics ( observedN=24, expectedBG=28.0, bgError=8.2 )

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_300MHT450", index = None, dataset="6NJet8_1250HT1500_300MHT450" )
T2.efficiencyMap.setStatistics ( observedN=5, expectedBG=9.4, bgError=3.6 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_300MHT450", index = None, dataset="6NJet8_1250HT1500_300MHT450" )
T1tttt.efficiencyMap.setStatistics ( observedN=5, expectedBG=9.4, bgError=3.6 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_300MHT450", index = None, dataset="6NJet8_1250HT1500_300MHT450" )
T1.efficiencyMap.setStatistics ( observedN=5, expectedBG=9.4, bgError=3.6 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_300MHT450", index = None, dataset="6NJet8_1250HT1500_300MHT450" )
T5WZ.efficiencyMap.setStatistics ( observedN=5, expectedBG=9.4, bgError=3.6 )

databaseCreator.create( True )





T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_450MHTinf", index = None, dataset="6NJet8_1250HT1500_450MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=2, expectedBG=0.5, bgError=0.4 ) #check error
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_450MHTinf", index = None, dataset="6NJet8_1250HT1500_450MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=2, expectedBG=0.5, bgError=0.4 ) #check error
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_450MHTinf", index = None, dataset="6NJet8_1250HT1500_450MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=2, expectedBG=0.5, bgError=0.4 ) #check error
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_1250HT1500_450MHTinf", index = None, dataset="6NJet8_1250HT1500_450MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=2, expectedBG=0.5, bgError=0.4 ) #check error

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_1500HTinf_200MHT300", index = None, dataset="6NJet8_1500HTinf_200MHT300" )
T2.efficiencyMap.setStatistics ( observedN=18, expectedBG=21.1, bgError=8.1 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_1500HTinf_200MHT300", index = None, dataset="6NJet8_1500HTinf_200MHT300" )
T1tttt.efficiencyMap.setStatistics ( observedN=18, expectedBG=21.1, bgError=8.1 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_1500HTinf_200MHT300", index = None, dataset="6NJet8_1500HTinf_200MHT300" )
T1.efficiencyMap.setStatistics ( observedN=18, expectedBG=21.1, bgError=8.1 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_1500HTinf_200MHT300", index = None, dataset="6NJet8_1500HTinf_200MHT300" )
T5WZ.efficiencyMap.setStatistics ( observedN=18, expectedBG=21.1, bgError=8.1 )

databaseCreator.create( True )







T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_6NJet8_1500HTinf_300MHTinf", index = None, dataset="6NJet8_1500HTinf_300MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=3, expectedBG=7.9, bgError=3.6 )
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_6NJet8_1500HTinf_300MHTinf", index = None, dataset="6NJet8_1500HTinf_300MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=3, expectedBG=7.9, bgError=3.6 )
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_6NJet8_1500HTinf_300MHTinf", index = None, dataset="6NJet8_1500HTinf_300MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=3, expectedBG=7.9, bgError=3.6 )
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_6NJet8_1500HTinf_300MHTinf", index = None, dataset="6NJet8_1500HTinf_300MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=3, expectedBG=7.9, bgError=3.6 )

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_8NJetinf_500HT800_200MHTinf", index = None, dataset="8NJetinf_500HT800_200MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=8, expectedBG=4.8, bgError=2.1 )#### check error
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_8NJetinf_500HT800_200MHTinf", index = None, dataset="8NJetinf_500HT800_200MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=8, expectedBG=4.8, bgError=2.1 )#### check error
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_8NJetinf_500HT800_200MHTinf", index = None, dataset="8NJetinf_500HT800_200MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=8, expectedBG=4.8, bgError=2.1 )#### check error
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_8NJetinf_500HT800_200MHTinf", index = None, dataset="8NJetinf_500HT800_200MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=8, expectedBG=4.8, bgError=2.1 )#### check error

databaseCreator.create( True )



T2.efficiencyMap.setSource("orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_8NJetinf_800HT1000_200MHTinf", index = None, dataset="8NJetinf_800HT1000_200MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=9, expectedBG=8.3, bgError=3.3 )#### check error
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_8NJetinf_800HT1000_200MHTinf", index = None, dataset="8NJetinf_800HT1000_200MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=9, expectedBG=8.3, bgError=3.3 )#### check error
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_8NJetinf_800HT1000_200MHTinf", index = None, dataset="8NJetinf_800HT1000_200MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=9, expectedBG=8.3, bgError=3.3 )#### check error
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_8NJetinf_800HT1000_200MHTinf", index = None, dataset="8NJetinf_800HT1000_200MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=9, expectedBG=8.3, bgError=3.3 )#### check error

databaseCreator.create( True )





T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_8NJetinf_1000HT1250_200MHTinf", index = None, dataset="8NJetinf_1000HT1250_200MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=8, expectedBG=5.6, bgError=2.1 )#### check error
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_8NJetinf_1000HT1250_200MHTinf", index = None, dataset="8NJetinf_1000HT1250_200MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=8, expectedBG=5.6, bgError=2.1 )#### check error
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_8NJetinf_1000HT1250_200MHTinf", index = None, dataset="8NJetinf_1000HT1250_200MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=8, expectedBG=5.6, bgError=2.1 )#### check error
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_8NJetinf_1000HT1250_200MHTinf", index = None, dataset="8NJetinf_1000HT1250_200MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=8, expectedBG=5.6, bgError=2.1 )#### check error

databaseCreator.create( True )




T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_8NJetinf_1250HT1500_200MHTinf", index = None, dataset="8NJetinf_1250HT1500_200MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=5, expectedBG=3.3, bgError=1.7 )#### check error
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_8NJetinf_1250HT1500_200MHTinf", index = None, dataset="8NJetinf_1250HT1500_200MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=5, expectedBG=3.3, bgError=1.7 )#### check error
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_8NJetinf_1250HT1500_200MHTinf", index = None, dataset="8NJetinf_1250HT1500_200MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=5, expectedBG=3.3, bgError=1.7 )#### check error
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_8NJetinf_1250HT1500_200MHTinf", index = None, dataset="8NJetinf_1250HT1500_200MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=5, expectedBG=3.3, bgError=1.7 )#### check error

databaseCreator.create( True )





T2.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "h_EffAcc_8NJetinf_1500HTinf_200MHTinf", index = None, dataset="8NJetinf_1500HTinf_200MHTinf" )
T2.efficiencyMap.setStatistics ( observedN=2, expectedBG=3.3, bgError=1.7 )#### check error
T1tttt.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "h_EffAcc_8NJetinf_1500HTinf_200MHTinf", index = None, dataset="8NJetinf_1500HTinf_200MHTinf" )
T1tttt.efficiencyMap.setStatistics ( observedN=2, expectedBG=3.3, bgError=1.7 )#### check error
T1.efficiencyMap.setSource(  "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "h_EffAcc_8NJetinf_1500HTinf_200MHTinf", index = None, dataset="8NJetinf_1500HTinf_200MHTinf" )
T1.efficiencyMap.setStatistics ( observedN=2, expectedBG=3.3, bgError=1.7 )#### check error
T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_8NJetinf_1500HTinf_200MHTinf", index = None, dataset="8NJetinf_1500HTinf_200MHTinf" )
T5WZ.efficiencyMap.setStatistics ( observedN=2, expectedBG=3.3, bgError=1.7 )#### check error

databaseCreator.create( True )













