#!/usr/bin/env python3

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""

from smodels_utils.dataPreparation.argParser import getParserArgs
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z

getParserArgs()

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2018-12')
# info.comment = ''
info.sqrts = '13.0'
info.private = False
info.lumi = '139.'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-12/'
info.prettyName = '0 leptons + jets + Etmiss'
info.implementedBy = 'WW'

# https://www.hepdata.net/record/ins1793461-v2

SRs = [ "SRAT0", "SRATT", "SRATW", "SRBT0", "SRBTT", "SRBTW", "SRC1", "SRC2", "SRC3", "SRC4",
        "SRC5", "SRD0", "SRD1", "SRD2" ]

yields = { "SRATT": { "obsN": 4, "bgExp": 3.2, "expErr": 0.5, "ul": 0.04 },
           "SRATW": { "obsN": 8, "bgExp": 5.6, "expErr": 0.7, "ul": 0.06 },
           "SRAT0": { "obsN": 11, "bgExp": 17.3, "expErr": 1.7, "ul": 0.05 },
           "SRBTT": { "obsN": 67, "bgExp": 46, "expErr": 7, "ul": 0.28 },
           "SRBTW": { "obsN": 84, "bgExp": 81, "expErr": 7, "ul": 0.21 },
           "SRBT0": { "obsN": 292, "bgExp": 276, "expErr": 24, "ul": 0.51 },
           "SRC1": { "obsN": 53, "bgExp": 46, "expErr": 12, "ul": 0.19 },
           "SRC2": { "obsN": 57, "bgExp": 52, "expErr":  9, "ul": 0.24 },
           "SRC3": { "obsN": 38, "bgExp": 32, "expErr":  7, "ul": 0.17 },
           "SRC4": { "obsN":  9, "bgExp": 11.8, "expErr":  3.1, "ul": 0.06 },
           "SRC5": { "obsN":  4, "bgExp":  2.5, "expErr":  0.7, "ul": 0.05 },
           "SRD0": { "obsN": 5, "bgExp": 6.9, "expErr": 1.3, "ul": 0.04 }, 
           "SRD1": { "obsN": 4, "bgExp": 3.1, "expErr": 1.0, "ul": 0.04 }, 
           "SRD2": { "obsN": 10, "bgExp": 12.2, "expErr": 1.5, "ul": 0.05 } }

figures = { "SRATT": ( "03a", "03b" ), "SRATW": ( "04a", "04b" ), 
            "SRAT0": ( "05a", "05b" ), "SRBTT": ( "06a", "06b" ),
            "SRBTW": ( "07a", "07b" ), "SRBT0": ( "08a", "08b" ),
            "SRC1": ( "09a", "09b" ), "SRC2": ( "10a", "10b" ),
            "SRC3": ( "11a", "11b" ), "SRC4": ( "12a", "12b" ), 
            "SRC5": ( "13a", "13b" ), "SRD0": ( "14a", "14b" ),
            "SRD1": ( "15a", "15b" ), "SRD2": ( "16a", "16b" )
}

base = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-12/"

for SR,data in yields.items():
    #+++++++ dataset block ++++++++++++++
    dataset = DataSetInput( SR )
    dataset.setInfo( dataType = 'efficiencyMap', dataId = SR, observedN = data["obsN"],
                     expectedBG=data["bgExp"], bgError = data["expErr"],
                     upperLimit = '%s*fb' % data["ul"] )

    T2tt = dataset.addTxName('T2tt')
    # T2tt.checked ="VM"
    T2tt.constraint ="[[['t']],[['t']]]"
    T2tt.conditionDescription ="None"
    T2tt.condition ="None"
    T2tt.source = 'ATLAS'
    #+++++++ next mass plane block ++++++++++++++
    T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
    T2tt_1.figure = ";".join ( [ f'Figaux. {fignr}' for fignr in figures[SR] ] )
    T2tt_1.figureUrl = ";".join ( [ f'{base}figaux_{fignr}.png' for fignr in figures[SR] ] )
    # T2tt_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-12/figaux_03a.png; https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-12/figaux_03b.png'
    T2tt_1.dataUrl = 'https://www.hepdata.net/record/ins1793461?version=2'
    multiplier = "/100000"
    if "SRC" in SR or "SRD" in SR:
        multiplier = "/10000000"
    T2tt_1.setSources(dataLabels= [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1', 'efficiencyMap' ],
            dataFiles= [ 'orig/stop_obs.csv', 'orig/stop_obs_up.csv', 'orig/stop_obs_down.csv', 'orig/stop_exp.csv', 'orig/stop_exp_up.csv', 'orig/stop_exp_down.csv', ( f'orig/Acc_{SR}.csv', f'orig/Eff_{SR}.csv' ) ],
            dataFormats= ['csv']*6+["mcsv"], units= [ None ]*6 + [ multiplier ] )

    T2ttoff = dataset.addTxName('T2ttoff')
    # T2ttoff.checked ="VM"
    T2ttoff.constraint ="[[['b','W']],[['b','W']]]"
    T2ttoff.conditionDescription ="None"
    T2ttoff.condition ="None"
    T2ttoff.massConstraint = [['80 <= dm <= 169.0'], ['80 <= dm <= 169.0']] 
    T2ttoff.source = 'ATLAS'
    #+++++++ next mass plane block ++++++++++++++
    T2ttoff_1 = T2ttoff.addMassPlane(2*[[x, y]])
    #T2ttoff_1.figure = 'Fig. 01a'
    #T2ttoff_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-12/figaux_01b.png'
    T2ttoff_1.figure = ";".join ( [ f'Figaux. {fignr}' for fignr in figures[SR] ] )
    T2ttoff_1.figureUrl = ";".join ( [ f'{base}figaux_{fignr}.png' for fignr in figures[SR] ] )
    T2ttoff_1.dataUrl = 'https://www.hepdata.net/record/ins1793461?version=2'
    T2ttoff_1.setSources(dataLabels= [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1', 'efficiencyMap' ],
            dataFiles= [ 'orig/stopoff_obs.csv', 'orig/stopoff_obs_up.csv', 'orig/stopoff_obs_down.csv', 'orig/stopoff_exp.csv', 'orig/stopoff_exp_up.csv', 'orig/stopoff_exp_down.csv', ( f'orig/Acc_{SR}.csv', f'orig/Eff_{SR}.csv' ) ],
            dataFormats= ['csv']*6+["mcsv"], units= [ None ]*6 + [ multiplier ] )

    T2bbffff = dataset.addTxName('T2bbffff')
    # T2bbffff.checked ="VM"
    T2bbffff.constraint = "[[['b', 'jet','jet']],[['b', 'jet','jet']]]"
    T2bbffff.conditionDescription ="None"
    T2bbffff.condition ="None"
    T2bbffff.massConstraint = [['dm < 80'], ['dm < 80']] 
    T2bbffff.source = 'ATLAS'
    #+++++++ next mass plane block ++++++++++++++
    T2bbffff_1 = T2bbffff.addMassPlane(2*[[x, y]])
    #T2bbffff_1.figure = 'Fig. 01a'
    #T2bbffff_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-12/figaux_01b.png'
    T2bbffff_1.figure = ";".join ( [ f'Figaux. {fignr}' for fignr in figures[SR] ] )
    T2bbffff_1.figureUrl = ";".join ( [ f'{base}figaux_{fignr}.png' for fignr in figures[SR] ] )
    T2bbffff_1.dataUrl = 'https://www.hepdata.net/record/ins1793461?version=2'
    T2bbffff_1.setSources(dataLabels= [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1', 'efficiencyMap' ],
            dataFiles= [ 'orig/stop_obs.csv', 'orig/stop_obs_up.csv', 'orig/stop_obs_down.csv', 'orig/stop_exp.csv', 'orig/stop_exp_up.csv', 'orig/stop_exp_down.csv', ( f'orig/Acc_{SR}.csv', f'orig/Eff_{SR}.csv' ) ],
            dataFormats= ['csv']*6+["mcsv"], units= [ None ]*6 + [ "*2.25e-07" ] )


databaseCreator.create()
