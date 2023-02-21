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
info                 = MetaInfoInput('ATLAS-SUSY-2019-08')
info.url             = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/'
info.sqrts             = 13
info.lumi             = 139.
info.prettyName     = '1L + higgs + Etmiss (EWino)'
info.private         = False
info.arxiv             = 'arXiv:1909.09226'
info.contact         = 'atlas-phys-susy-conveners@cern.ch'
info.publication     = ''
info.implementedBy   = "Gael Alguero"
info.datasetOrder    = '"SR_HM_Low_MCT", "SR_HM_Med_MCT", "SR_HM_High_MCT", "SR_MM_Low_MCT", "SR_MM_Med_MCT", "SR_MM_High_MCT", "SR_LM_Low_MCT", "SR_LM_Med_MCT", "SR_LM_High_MCT"'
info.jsonFiles       = '{"BkgOnly.json" : ["SR_HM_Low_MCT", "SR_HM_Med_MCT", "SR_HM_High_MCT", "SR_MM_Low_MCT", "SR_MM_Med_MCT", "SR_MM_High_MCT", "SR_LM_Low_MCT", "SR_LM_Med_MCT", "SR_LM_High_MCT"]}'

# +++++++ exclusive SR_LM ++++++++++++++
# +++++++ dataset block ++++++++++++++
SR_LM_Low_MCT = DataSetInput('SR_LM_Low_MCT')
SR_LM_Low_MCT.setInfo( dataType = 'efficiencyMap', dataId = 'SR_LM_Low_MCT',
                 observedN = 16, expectedBG = 8.8, bgError = 2.8)

# +++++++ next txName block ++++++++++++++
TChiWH = SR_LM_Low_MCT.addTxName("TChiWH" )
TChiWH.checked              = 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription = None
TChiWH.condition            = None
TChiWH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWH1           = TChiWH.addMassPlane( 2*[[ x, y ]] )
TChiWH1.figure    = "fig_07b"
TChiWH1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/figaux_07b.png"
TChiWH1.dataUrl   = "https://www.hepdata.net/record/ins1755298?version=2&table=Eff_table_SR_LM_Low_MCT"

TChiWH1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion'],
    dataFiles = [ "orig/Observedlimit1lbb.csv",
                  "orig/Observedlimit1lbb(Up).csv",
                  "orig/Observedlimit1lbb(Down).csv",
                  "orig/Expectedlimit1lbb.csv"],
    units = [ None, None, None, None ], dataFormats = [ 'csv' ]*4 )
TChiWH1.addSource("efficiencyMap", ("orig/Acc_table_SR_LM_Low_MCT.csv", "orig/Eff_table_SR_LM_Low_MCT.csv"),
                  unit = "/10000", dataFormat = "mcsv")

# +++++++ dataset block ++++++++++++++
SR_LM_Med_MCT = DataSetInput('SR_LM_Med_MCT')
SR_LM_Med_MCT.setInfo( dataType = 'efficiencyMap', dataId = 'SR_LM_Med_MCT',
                 observedN = 11, expectedBG = 11.3, bgError = 3.1)

# +++++++ next txName block ++++++++++++++
TChiWH = SR_LM_Med_MCT.addTxName("TChiWH" )
TChiWH.checked              = 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription = None
TChiWH.condition            = None
TChiWH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWH1           = TChiWH.addMassPlane( 2*[[ x, y ]] )
TChiWH1.figure    = "fig_07c"
TChiWH1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/figaux_07c.png"
TChiWH1.dataUrl   = "https://www.hepdata.net/record/ins1755298?version=2&table=Eff_table_SR_LM_Med_MCT"

TChiWH1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion'],
    dataFiles = [ "orig/Observedlimit1lbb.csv",
                  "orig/Observedlimit1lbb(Up).csv",
                  "orig/Observedlimit1lbb(Down).csv",
                  "orig/Expectedlimit1lbb.csv"],
    units = [ None, None, None, None ], dataFormats = [ 'csv' ]*4 )
TChiWH1.addSource("efficiencyMap", ("orig/Acc_table_SR_LM_Med_MCT.csv", "orig/Eff_table_SR_LM_Med_MCT.csv"),
                  unit = "/10000", dataFormat = "mcsv")

# +++++++ dataset block ++++++++++++++
SR_LM_High_MCT = DataSetInput('SR_LM_High_MCT')
SR_LM_High_MCT.setInfo( dataType = 'efficiencyMap', dataId = 'SR_LM_High_MCT',
                 observedN = 7, expectedBG = 7.3, bgError = 1.5)

# +++++++ next txName block ++++++++++++++
TChiWH = SR_LM_High_MCT.addTxName("TChiWH" )
TChiWH.checked              = 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription = None
TChiWH.condition            = None
TChiWH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWH1           = TChiWH.addMassPlane( 2*[[ x, y ]] )
TChiWH1.figure    = "fig_07d"
TChiWH1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/figaux_07d.png"
TChiWH1.dataUrl   = "https://www.hepdata.net/record/ins1755298?version=2&table=Eff_table_SR_LM_High_MCT"

TChiWH1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion'],
    dataFiles = [ "orig/Observedlimit1lbb.csv",
                  "orig/Observedlimit1lbb(Up).csv",
                  "orig/Observedlimit1lbb(Down).csv",
                  "orig/Expectedlimit1lbb.csv"],
    units = [ None, None, None, None ], dataFormats = [ 'csv' ]*4 )
TChiWH1.addSource("efficiencyMap", ("orig/Acc_table_SR_LM_High_MCT.csv", "orig/Eff_table_SR_LM_High_MCT.csv"),
                  unit = "/10000", dataFormat = "mcsv")

# +++++++ exclusive SR_MM ++++++++++++++
# +++++++ dataset block ++++++++++++++
SR_MM_Low_MCT = DataSetInput('SR_MM_Low_MCT')
SR_MM_Low_MCT.setInfo( dataType = 'efficiencyMap', dataId = 'SR_MM_Low_MCT',
                 observedN = 4, expectedBG = 4.6, bgError = 1.7)

# +++++++ next txName block ++++++++++++++
TChiWH = SR_MM_Low_MCT.addTxName("TChiWH" )
TChiWH.checked              = 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription = None
TChiWH.condition            = None
TChiWH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWH1           = TChiWH.addMassPlane( 2*[[ x, y ]] )
TChiWH1.figure    = "fig_08b"
TChiWH1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/figaux_08b.png"
TChiWH1.dataUrl   = "https://www.hepdata.net/record/ins1755298?version=2&table=Eff_table_SR_MM_Low_MCT"

TChiWH1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion'],
    dataFiles = [ "orig/Observedlimit1lbb.csv",
                  "orig/Observedlimit1lbb(Up).csv",
                  "orig/Observedlimit1lbb(Down).csv",
                  "orig/Expectedlimit1lbb.csv"],
    units = [ None, None, None, None ], dataFormats = [ 'csv' ]*4 )
TChiWH1.addSource("efficiencyMap", ("orig/Acc_table_SR_MM_Low_MCT.csv", "orig/Eff_table_SR_MM_Low_MCT.csv"),
                  unit = "/10000", dataFormat = "mcsv")

# +++++++ dataset block ++++++++++++++
SR_MM_Med_MCT = DataSetInput('SR_MM_Med_MCT')
SR_MM_Med_MCT.setInfo( dataType = 'efficiencyMap', dataId = 'SR_MM_Med_MCT',
                 observedN = 7, expectedBG = 2.6, bgError = 1.3)

# +++++++ next txName block ++++++++++++++
TChiWH = SR_MM_Med_MCT.addTxName("TChiWH" )
TChiWH.checked              = 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription = None
TChiWH.condition            = None
TChiWH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWH1           = TChiWH.addMassPlane( 2*[[ x, y ]] )
TChiWH1.figure    = "fig_08c"
TChiWH1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/figaux_08c.png"
TChiWH1.dataUrl   = "https://www.hepdata.net/record/ins1755298?version=2&table=Eff_table_SR_MM_Med_MCT"

TChiWH1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion'],
    dataFiles = [ "orig/Observedlimit1lbb.csv",
                  "orig/Observedlimit1lbb(Up).csv",
                  "orig/Observedlimit1lbb(Down).csv",
                  "orig/Expectedlimit1lbb.csv"],
    units = [ None, None, None, None ], dataFormats = [ 'csv' ]*4 )
TChiWH1.addSource("efficiencyMap", ("orig/Acc_table_SR_MM_Med_MCT.csv", "orig/Eff_table_SR_MM_Med_MCT.csv"),
                  unit = "/10000", dataFormat = "mcsv")

# +++++++ dataset block ++++++++++++++
SR_MM_High_MCT = DataSetInput('SR_MM_High_MCT')
SR_MM_High_MCT.setInfo( dataType = 'efficiencyMap', dataId = 'SR_MM_High_MCT',
                 observedN = 2, expectedBG = 1.4, bgError = 0.6)

# +++++++ next txName block ++++++++++++++
TChiWH = SR_MM_High_MCT.addTxName("TChiWH" )
TChiWH.checked              = 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription = None
TChiWH.condition            = None
TChiWH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWH1           = TChiWH.addMassPlane( 2*[[ x, y ]] )
TChiWH1.figure    = "fig_08d"
TChiWH1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/figaux_08d.png"
TChiWH1.dataUrl   = "https://www.hepdata.net/record/ins1755298?version=2&table=Eff_table_SR_MM_High_MCT"

TChiWH1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion'],
    dataFiles = [ "orig/Observedlimit1lbb.csv",
                  "orig/Observedlimit1lbb(Up).csv",
                  "orig/Observedlimit1lbb(Down).csv",
                  "orig/Expectedlimit1lbb.csv"],
    units = [ None, None, None, None ], dataFormats = [ 'csv' ]*4 )
TChiWH1.addSource("efficiencyMap", ("orig/Acc_table_SR_MM_High_MCT.csv", "orig/Eff_table_SR_MM_High_MCT.csv"),
                  unit = "/10000", dataFormat = "mcsv")

# +++++++ exclusive SR_HM ++++++++++++++
# +++++++ dataset block ++++++++++++++
SR_HM_Low_MCT = DataSetInput('SR_HM_Low_MCT')
SR_HM_Low_MCT.setInfo( dataType = 'efficiencyMap', dataId = 'SR_HM_Low_MCT',
                 observedN = 6, expectedBG = 4.1, bgError = 1.9)

# +++++++ next txName block ++++++++++++++
TChiWH = SR_HM_Low_MCT.addTxName("TChiWH" )
TChiWH.checked              = 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription = None
TChiWH.condition            = None
TChiWH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWH1           = TChiWH.addMassPlane( 2*[[ x, y ]] )
TChiWH1.figure    = "fig_09b"
TChiWH1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/figaux_09b.png"
TChiWH1.dataUrl   = "https://www.hepdata.net/record/ins1755298?version=2&table=Eff_table_SR_HM_Low_MCT"

TChiWH1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion'],
    dataFiles = [ "orig/Observedlimit1lbb.csv",
                  "orig/Observedlimit1lbb(Up).csv",
                  "orig/Observedlimit1lbb(Down).csv",
                  "orig/Expectedlimit1lbb.csv"],
    units = [ None, None, None, None ], dataFormats = [ 'csv' ]*4 )
TChiWH1.addSource("efficiencyMap", ("orig/Acc_table_SR_HM_Low_MCT.csv", "orig/Eff_table_SR_HM_Low_MCT.csv"),
                  unit = "/10000", dataFormat = "mcsv")

# +++++++ dataset block ++++++++++++++
SR_HM_Med_MCT = DataSetInput('SR_HM_Med_MCT')
SR_HM_Med_MCT.setInfo( dataType = 'efficiencyMap', dataId = 'SR_HM_Med_MCT',
                 observedN = 5, expectedBG = 2.9, bgError = 1.3)

# +++++++ next txName block ++++++++++++++
TChiWH = SR_HM_Med_MCT.addTxName("TChiWH" )
TChiWH.checked              = 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription = None
TChiWH.condition            = None
TChiWH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWH1           = TChiWH.addMassPlane( 2*[[ x, y ]] )
TChiWH1.figure    = "fig_09c"
TChiWH1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/figaux_09c.png"
TChiWH1.dataUrl   = "https://www.hepdata.net/record/ins1755298?version=2&table=Eff_table_SR_HM_Med_MCT"

TChiWH1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion'],
    dataFiles = [ "orig/Observedlimit1lbb.csv",
                  "orig/Observedlimit1lbb(Up).csv",
                  "orig/Observedlimit1lbb(Down).csv",
                  "orig/Expectedlimit1lbb.csv"],
    units = [ None, None, None, None ], dataFormats = [ 'csv' ]*4 )
TChiWH1.addSource("efficiencyMap", ("orig/Acc_table_SR_HM_Med_MCT.csv", "orig/Eff_table_SR_HM_Med_MCT.csv"),
                  unit = "/10000", dataFormat = "mcsv")

# +++++++ dataset block ++++++++++++++
SR_HM_High_MCT = DataSetInput('SR_HM_High_MCT')
SR_HM_High_MCT.setInfo( dataType = 'efficiencyMap', dataId = 'SR_HM_High_MCT',
                 observedN = 3, expectedBG = 1.1, bgError = 0.5)

# +++++++ next txName block ++++++++++++++
TChiWH = SR_HM_High_MCT.addTxName("TChiWH" )
TChiWH.checked              = 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription = None
TChiWH.condition            = None
TChiWH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWH1           = TChiWH.addMassPlane( 2*[[ x, y ]] )
TChiWH1.figure    = "fig_09d"
TChiWH1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/figaux_09d.png"
TChiWH1.dataUrl   = "https://www.hepdata.net/record/ins1755298?version=2&table=Eff_table_SR_HM_High_MCT"

TChiWH1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion'],
    dataFiles = [ "orig/Observedlimit1lbb.csv",
                  "orig/Observedlimit1lbb(Up).csv",
                  "orig/Observedlimit1lbb(Down).csv",
                  "orig/Expectedlimit1lbb.csv"],
    units = [ None, None, None, None ], dataFormats = [ 'csv' ]*4 )
TChiWH1.addSource("efficiencyMap", ("orig/Acc_table_SR_HM_High_MCT.csv", "orig/Eff_table_SR_HM_High_MCT.csv"),
                  unit = "/10000", dataFormat = "mcsv")

databaseCreator.create()
