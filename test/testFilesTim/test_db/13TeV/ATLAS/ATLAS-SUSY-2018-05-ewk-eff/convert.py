#!/usr/bin/env python3

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse

argparser = argparse.ArgumentParser(description =
    'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath',
    help = 'path to the package smodels_utils',\
    type = str )
argparser.add_argument ('-smodelsPath', '--smodelsPath',
    help = 'path to the package smodels_utils',\
    type = str )
argparser.add_argument ('-no', '--noUpdate',
    help = 'do not update the lastUpdate field.',\
    action= "store_true" )
argparser.add_argument ('-r', '--resetValidation',
    help = 'reset the validation flag',\
    action= "store_true" )

args = argparser.parse_args()

if args.noUpdate:
    os.environ["SMODELS_NOUPDATE"]="1"
if args.resetValidation:
    os.environ["SMODELS_RESETVALIDATION"]="1"

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
info                 = MetaInfoInput('ATLAS-SUSY-2018-05-ewk-eff')
info.url             = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/'
info.sqrts             = 13
info.lumi             = 139.
info.prettyName     = '2L + jets + MET'
info.private         = False
info.arxiv             = 'arXiv:2204.13072'
info.contact         = 'atlas-phys-susy-conveners@cern.ch'
info.publication     = ''
info.implementedBy   = "Sahana Narasimha"
info.datasetOrder    = ' "SRHigh4_cuts", "SRHigh8_1_cuts", "SRHigh8_2_cuts", "SRHigh16_1_cuts", "SRHigh16_2_cuts", "SRllbb_cuts", "SRInt_1_cuts", "SRInt_2_cuts", "SRLow_1_cuts", "SRLow_2_cuts", "SRLow2_cuts", "SROffShell_1_cuts", "SROffShell_2_cuts"'
info.jsonFiles       =  '{ "llh_jsons/ewk_signal_bkgonly.json" : ["SRHigh4_cuts", "SRHigh8_1_cuts", "SRHigh8_2_cuts", "SRHigh16_1_cuts", "SRHigh16_2_cuts", "SRllbb_cuts", "SRInt_1_cuts", "SRInt_2_cuts", "SRLow_1_cuts", "SRLow_2_cuts", "SRLow2_cuts", "SROffShell_1_cuts", "SROffShell_2_cuts"]}'



# +++++++ exclusive EWK-SR ++++++++++++++

# +++++++ dataset block ++++++++++++++
SROffShell_1_cuts = DataSetInput('SROffShell_1_cuts')
SROffShell_1_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SROffShell_1_cuts',
                 observedN = 6, expectedBG = 9.2, bgError = 1.7)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SROffShell_1_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SROffShell_1_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_09b"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_09b.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t130"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SROffShell_1_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_09b"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_09b.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t130"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SROffShell_1_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"

# +++++++ dataset block ++++++++++++++
SROffShell_2_cuts = DataSetInput('SROffShell_2_cuts')
SROffShell_2_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SROffShell_2_cuts',
                 observedN = 15, expectedBG = 12.5, bgError = 1.9)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SROffShell_2_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SROffShell_2_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_09d"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_09d.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t132"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SROffShell_2_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_09d"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_09d.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t132"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SROffShell_2_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"


# +++++++ dataset block ++++++++++++++
SRLow_1_cuts = DataSetInput('SRLow_1_cuts')
SRLow_1_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SRLow_1_cuts',
                 observedN = 10, expectedBG = 12.8, bgError = 3.4)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SRLow_1_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SRLow_1_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'


# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_09f"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_09f.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t134"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRLow_1_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_09f"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_09f.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t134"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRLow_1_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"


# +++++++ dataset block ++++++++++++++
SRLow_2_cuts = DataSetInput('SRLow_2_cuts')
SRLow_2_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SRLow_2_cuts',
                 observedN = 8, expectedBG = 10.5, bgError = 2.5)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SRLow_2_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SRLow_2_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'


# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_09h"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_09h.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t136"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRLow_2_cuts_TChiWZ.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_09h"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_09h.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t136"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRLow_2_cuts_TChiWZ.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"

# +++++++ dataset block ++++++++++++++
SRLow2_cuts = DataSetInput('SRLow2_cuts')
SRLow2_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SRLow2_cuts',
                 observedN = 8, expectedBG = 9.0, bgError = 4.0)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SRLow2_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SRLow2_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_10b"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_10b.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t138"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRLow2_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_10b"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_10b.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t138"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRLow2_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"

# +++++++ dataset block ++++++++++++++
SRInt_1_cuts = DataSetInput('SRInt_1_cuts')
SRInt_1_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SRInt_1_cuts',
                 observedN = 24, expectedBG = 22.8, bgError = 3.5)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SRInt_1_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SRInt_1_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'


# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_10d"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_10d.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t140"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRInt_1_cuts_TChiWZ.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_10d"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_10d.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t140"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRInt_1_cuts_TChiWZ.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"

# +++++++ dataset block ++++++++++++++
SRInt_2_cuts = DataSetInput('SRInt_2_cuts')
SRInt_2_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SRInt_2_cuts',
                 observedN = 14, expectedBG = 10.1, bgError = 1.0)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SRInt_2_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SRInt_2_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_10f"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_10f.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t142"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRInt_2_cuts_TChiWZ.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_10f"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_10f.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t142"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRInt_2_cuts_TChiWZ.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"

# +++++++ dataset block ++++++++++++++
SRHigh16_1_cuts = DataSetInput('SRHigh16_1_cuts')
SRHigh16_1_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SRHigh16_1_cuts',
                 observedN = 4, expectedBG = 3.9, bgError = 0.7)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SRHigh16_1_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SRHigh16_1_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_11b"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_11b.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t144"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRHigh16_1_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_11b"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_11b.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t144"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRHigh16_1_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"

# +++++++ dataset block ++++++++++++++
SRHigh16_2_cuts = DataSetInput('SRHigh16_2_cuts')
SRHigh16_2_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SRHigh16_2_cuts',
                 observedN = 3, expectedBG = 3.4, bgError = 0.9)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SRHigh16_2_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SRHigh16_2_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_11d"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_11d.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t146"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRHigh16_2_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_11d"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_11d.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t146"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRHigh16_2_cuts_TChiWZ_Bin_1.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"


# +++++++ dataset block ++++++++++++++
SRHigh8_1_cuts = DataSetInput('SRHigh8_1_cuts')
SRHigh8_1_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SRHigh8_1_cuts',
                 observedN = 0, expectedBG = 2.0, bgError = 0.23)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SRHigh8_1_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SRHigh8_1_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_11f"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_11f.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t148"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRHigh8_1_cuts_TChiWZ.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_11f"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_11f.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t148"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRHigh8_1_cuts_TChiWZ.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"


# +++++++ dataset block ++++++++++++++
SRHigh8_2_cuts = DataSetInput('SRHigh8_2_cuts')
SRHigh8_2_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SRHigh8_2_cuts',
                 observedN = 0, expectedBG = 2.0, bgError = 0.33)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SRHigh8_2_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SRHigh8_2_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'


# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_11h"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_11h.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t150"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRHigh8_2_cuts_TChiWZ.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_11h"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_11h.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t150"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRHigh8_2_cuts_TChiWZ.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"

# +++++++ dataset block ++++++++++++++
SRHigh4_cuts = DataSetInput('SRHigh4_cuts')
SRHigh4_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SRHigh4_cuts',
                 observedN = 1, expectedBG = 0.85, bgError = 0.34)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SRHigh4_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*[[['l+','l-']],[['jet','jet']]]"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SRHigh4_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'


# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_12b"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_12b.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t152"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRHigh4_cuts_TChiWZ.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_12b"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_12b.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t152"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRHigh4_cuts_TChiWZ.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"

# +++++++ dataset block ++++++++++++++
SRllbb_cuts = DataSetInput('SRllbb_cuts')
SRllbb_cuts.setInfo( dataType = 'efficiencyMap', dataId = 'SRllbb_cuts',
                 observedN = 0, expectedBG = 0.58, bgError = 0.2)
                 
# +++++++ next txName block ++++++++++++++
TChiWZoff = SRllbb_cuts.addTxName("TChiWZoff" )
TChiWZoff.checked              = 'no'
TChiWZoff.constraint           = "22.08*([[['l+','l-']],[['jet','jet']]])"
TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription = None
TChiWZoff.condition            = None
TChiWZoff.source               = 'ATLAS'

TChiWZ = SRllbb_cuts.addTxName("TChiWZ" )
TChiWZ.checked              = 'no'
TChiWZ.constraint           = "[[['W']], [['Z']]]"
TChiWZ.conditionDescription = None
TChiWZ.condition            = None
TChiWZ.source               = 'ATLAS'


# +++++++ next mass plane block ++++++++++++++
TChiWZoff_1           = TChiWZoff.addMassPlane( 2*[[ x, y ]] )
TChiWZoff_1.figure    = "Fig_12d"
TChiWZoff_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/figaux_12d.png"
TChiWZoff_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t154"
TChiWZoff_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1' ],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZoff_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRllbb_cuts_TChiWZ.csv"),dataFormat = "csv")
TChiWZoff_1.efficiencyMap._unit = "22.08"

TChiWZ_1           = TChiWZ.addMassPlane( 2*[[ x, y ]] )
TChiWZ_1.figure    = "Fig_015a"
TChiWZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-05/fig_15a.png"
TChiWZ_1.dataUrl   = "https://doi.org/10.17182/hepdata.116034.v1/t154"
TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
    dataFiles = [ "orig/TChiWZoff_Observed_Limit.csv",
                  "orig/TChiWZoff_Observed_LimitP1.csv",
                  "orig/TChiWZoff_Observed_LimitM1.csv",
                  "orig/TChiWZoff_Expected_Limit.csv",
                  "orig/TChiWZoff_Expected_LimitP1.csv",
                  "orig/TChiWZoff_Expected_LimitM1.csv"],
    units = [ None, None, None, None, None, None], dataFormats = [ 'csv' ]*6 )
TChiWZ_1.addSource("efficiencyMap", ("orig/EffFromPatch_SRllbb_cuts_TChiWZ.csv"),dataFormat = "csv")
#TChiWZ_1.efficiencyMap._unit = "1.0/22.08"



databaseCreator.create()
