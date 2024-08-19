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
info                 = MetaInfoInput('ATLAS-SUSY-2018-31')
info.url             = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-31/'
info.sqrts             = 13
info.lumi             = 139.
info.prettyName     = 'higgs + b-jets + MET'
info.private         = False
info.arxiv             = 'arXiv:1908.03122'
info.contact         = 'atlas-phys-susy-conveners@cern.ch'
info.publication     = ''

# +++++++ dataset block ++++++++++++++
SRA_incl = DataSetInput('SRA_incl')
SRA_incl.setInfo( dataType = 'efficiencyMap', dataId = "SRA_incl",
                  observedN = 17, expectedBG=17.1, bgError = 2.8,
#                 upperLimit = '2.637E-01*fb', expectedUpperLimit = '2.639E-01*fb', ## FIXME this line is nonsense
                  jsonfile = "orig/BkgOnlyA.json" )

# +++++++ next txName block ++++++++++++++
T6bbHH = SRA_incl.addTxName("T6bbHH" )
T6bbHH.checked              = 'NO'
T6bbHH.constraint           = "[[['b'],['higgs']],[['b'],['higgs']]]"
T6bbHH.conditionDescription = None
T6bbHH.condition            = None
T6bbHH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
T6bbHHM60           = T6bbHH.addMassPlane( 2*[[ x, y, 60. ]] )
T6bbHHM60.figure    = "figaux_05b (acceptance) and figaux_9b (efficiency)"
T6bbHHM60.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-31/figaux_05b.png"
T6bbHHM60.dataUrl   = "https://www.hepdata.net/record/ins1748602?version=1&table=Efficiency_SRB and others"
T6bbHHM60.addSource ( "obsExclusion", "orig/HEPData-ins1748602-v1-M60_Obs.csv",
                      unit = None, dataFormat = "csv" )
T6bbHHM60.addSource ( "efficiencyMap", 
        ( "orig/Efficiency_SRA_incl_m60.csv", "orig/Acceptance_SRA_incl_m60.csv" ), 
        unit = "/10000", dataFormat = "mcsv" )

T6bbHHDM130           = T6bbHH.addMassPlane( 2*[[ x, y, y- 130. ]] )
T6bbHHDM130.figure    = "fig_08b"
T6bbHHDM130.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-31/fig_08b.png"
T6bbHHDM130.dataUrl   = "https://www.hepdata.net/record/ins1748602?version=1&table=DM130_bestCombined_UpperLimit"
T6bbHHDM130.addSource ( "obsExclusion", "orig/HEPData-ins1748602-v1-DM130_Obs.csv",
                        unit = None, dataFormat = "csv" )
T6bbHHDM130.addSource ( "efficiencyMap", ( "orig/Efficiency_SRA_incl_dm130.csv", "orig/Acceptance_SRA_incl_dm130.csv" ), unit = "/10000", 
                     dataFormat = "mcsv" )

# +++++++ dataset block ++++++++++++++
SRB = DataSetInput('SRB')
SRB.setInfo( dataType = 'efficiencyMap', dataId = "SRB", 
             observedN = 3, expectedBG=3.3, bgError = 0.9,
            # upperLimit = '2.637E-01*fb', expectedUpperLimit = '2.639E-01*fb', # FIXME this line is nonsense
             jsonfile = "orig/BkgOnlyB.json" )

# +++++++ next txName block ++++++++++++++
T6bbHH = SRB.addTxName("T6bbHH" )
T6bbHH.checked              = 'NO'
T6bbHH.constraint           = "[[['b'],['higgs']],[['b'],['higgs']]]"
T6bbHH.conditionDescription = None
T6bbHH.condition            = None
T6bbHH.source               = 'ATLAS'

"""
# +++++++ next mass plane block ++++++++++++++
T6bbHHM60           = T6bbHH.addMassPlane( 2*[[ x, y, 60. ]] )
T6bbHHM60.figure    = "figaux_05b (acceptance) and figaux_9b (efficiency)"
T6bbHHM60.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-31/figaux_05b.png"
T6bbHHM60.dataUrl   = "https://www.hepdata.net/record/ins1748602?version=1&table=Efficiency_SRB and others"
T6bbHHM60.addSource ( "obsExclusion", "orig/HEPData-ins1748602-v1-M60_Obs.csv",
                      unit = None, dataFormat = "csv" )
T6bbHHM60.addSource ( "efficiencyMap", 
        ( "orig/Efficiency_SRB.csv", "orig/Acceptance_SRB.csv" ), 
        unit = "/10000", dataFormat = "mcsv" )
"""

T6bbHHDM130           = T6bbHH.addMassPlane( 2*[[ x, y, y-130. ]] )
T6bbHHDM130.figure    = "fig_08b"
T6bbHHDM130.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-31/fig_08b.png"
T6bbHHDM130.dataUrl   = "https://www.hepdata.net/record/ins1748602?version=1&table=DM130_bestCombined_UpperLimit"
T6bbHHDM130.addSource ( "obsExclusion", "orig/HEPData-ins1748602-v1-DM130_Obs.csv",
                        unit = None, dataFormat = "csv" )
T6bbHHDM130.addSource ( "efficiencyMap", ( "orig/Efficiency_SRB.csv", "orig/Acceptance_SRB.csv" ), unit = "/10000", 
                     dataFormat = "mcsv" )

# +++++++ dataset block ++++++++++++++
SRC_incl = DataSetInput('SRC_incl')
SRC_incl.setInfo( dataType = 'efficiencyMap', dataId = "SRC_incl", 
                  observedN = 47, expectedBG=37.9, bgError = 6.2,
                 # upperLimit = '2.637E-01*fb', expectedUpperLimit = '2.639E-01*fb', # FIXME this line is nonsense
                  jsonfile = "orig/BkgOnlyC.json" )

# +++++++ next txName block ++++++++++++++
T6bbHH = SRC_incl.addTxName("T6bbHH" )
T6bbHH.checked              = 'NO'
T6bbHH.constraint           = "[[['b'],['higgs']],[['b'],['higgs']]]"
T6bbHH.conditionDescription = None
T6bbHH.condition            = None
T6bbHH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
T6bbHHM60           = T6bbHH.addMassPlane( 2*[[ x, y, 60. ]] )
T6bbHHM60.figure    = "figaux_05b (acceptance) and figaux_9b (efficiency)"
T6bbHHM60.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-31/figaux_05b.png"
T6bbHHM60.dataUrl   = "https://www.hepdata.net/record/ins1748602?version=1&table=Efficiency_SRB and others"
T6bbHHM60.addSource ( "obsExclusion", "orig/HEPData-ins1748602-v1-M60_Obs.csv",
                      unit = None, dataFormat = "csv" )
T6bbHHM60.addSource ( "efficiencyMap", 
        ( "orig/Efficiency_SRC_incl.csv", "orig/Acceptance_SRC_incl.csv" ), 
        unit = "/10000", dataFormat = "mcsv" )

"""
T6bbHHDM130           = T6bbHH.addMassPlane( 2*[[ x, y, y - 130. ]] )
T6bbHHDM130.figure    = "fig_08b"
T6bbHHDM130.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-31/fig_08b.png"
T6bbHHDM130.dataUrl   = "https://www.hepdata.net/record/ins1748602?version=1&table=DM130_bestCombined_UpperLimit"
T6bbHHDM130.addSource ( "obsExclusion", "orig/HEPData-ins1748602-v1-DM130_Obs.csv",
                        unit = None, dataFormat = "csv" )
T6bbHHDM130.addSource ( "efficiencyMap", ( "orig/Efficiency_SRC_incl.csv", "orig/Acceptance_SRC_incl.csv" ), unit = "%", 
                     dataFormat = "mcsv" )
"""

databaseCreator.create()
