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
info                 = MetaInfoInput('ATLAS-SUSY-2018-04')
info.url             = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-04/'
info.sqrts             = 13
info.lumi             = 139.
info.prettyName     = '2 hadronic taus'
info.private         = False
info.arxiv             = 'https://arxiv.org/abs/1911.06660'
info.publication = 'Phys. Rev. D 101 (2020) 032009'
info.publicationDOI = 'https://doi.org/10.1103/PhysRevD.101.032009'
info.contact         = 'atlas-phys-susy-conveners@cern.ch'
info.jsonFiles = "{'SRcombined.json': ['SRlow', 'SRhigh']}"

# +++++++ SRlow block ++++++++++++++
SRlow = DataSetInput('SRlow')
SRlow.setInfo( dataType = 'efficiencyMap', dataId = "SRlow", 
               observedN = 10, expectedBG=6.0, bgError=1.7,
               upperLimit = '0.08*fb', expectedUpperLimit = '0.055*fb',
)

# +++++++ next txName block ++++++++++++++
TStauStau = SRlow.addTxName("TStauStau" )
TStauStau.checked              = 'no'
TStauStau.constraint           = "[[['ta']],[['ta']]]"
TStauStau.conditionDescription = None
TStauStau.condition            = None
TStauStau.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TStauStau1           = TStauStau.addMassPlane( 2*[[ x, y ]] )
TStauStau1.figure    = "figaux_01a"
TStauStau1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-04/figaux_01a.png"
TStauStau1.dataUrl   = "https://www.hepdata.net/record/ins1765529?version=1&table=X-section U.L 1"
TStauStau1.setSources(dataLabels = [ 'obsExclusion', 'expExclusion', 'efficiencyMap' ],
    dataFiles = [ "orig/HEPData-ins1765529-v1-Exclusion_contour_1_(Obs.).csv", \
                  "orig/HEPData-ins1765529-v1-Exclusion_contour_1_(Exp.).csv",
                  "orig/HEPData-ins1765529-v1-Eff_Acc_table_SRlow.csv" ],
    units = [ None, None, "%" ], dataFormats = [ 'csv', 'csv', 'csv'] )

# +++++++ SRhigh block ++++++++++++++
SRhigh = DataSetInput('SRhigh')
SRhigh.setInfo( dataType = 'efficiencyMap', dataId = "SRhigh", 
               observedN = 7, expectedBG=10.2, bgError=3.3,
               upperLimit = '0.05*fb', expectedUpperLimit = '0.065*fb',
)

# +++++++ next txName block ++++++++++++++
TStauStau = SRhigh.addTxName("TStauStau" )
TStauStau.checked              = 'no'
TStauStau.constraint           = "[[['ta']],[['ta']]]"
TStauStau.conditionDescription = None
TStauStau.condition            = None
TStauStau.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TStauStau1           = TStauStau.addMassPlane( 2*[[ x, y ]] )
TStauStau1.figure    = "fig_07a"
TStauStau1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-04/fig_07a.png"
TStauStau1.dataUrl   = "https://www.hepdata.net/record/ins1765529?version=1&table=X-section U.L 1"
TStauStau1.setSources(dataLabels = [ 'obsExclusion', 'expExclusion', 'efficiencyMap' ],
    dataFiles = [ "orig/HEPData-ins1765529-v1-Exclusion_contour_1_(Obs.).csv", \
                  "orig/HEPData-ins1765529-v1-Exclusion_contour_1_(Exp.).csv",
                  "orig/HEPData-ins1765529-v1-Eff_Acc_table_SRhigh.csv" ],
    units = [ None, None, '%' ], dataFormats = [ 'csv', 'csv', 'csv'] )

databaseCreator.create()
