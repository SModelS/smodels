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

# https://www.hepdata.net/record/ins1289225

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

T2tt = dataset.addTxName('T2tt')
# T2tt.checked ="VM"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure = 'Fig. 01a'
T2tt_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-12/figaux_01a.png'
T2tt_1.dataUrl = 'https://www.hepdata.net/record/ins1793461?version=2&table=stop_xSecUpperLimit_obs'
T2tt_1.setSources(dataLabels= ['obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1', 'upperLimits', 'expectedUpperLimits' ],
        dataFiles= [ 'orig/stop_obs.csv', 'orig/stop_obs_up.csv', 'orig/stop_obs_down.csv', 'orig/stop_exp.csv', 'orig/stop_exp_up.csv', 'orig/stop_exp_down.csv', 'orig/stop_xSecUpperLimit_obs.csv', 'orig/stop_xSecUpperLimit_exp.csv' ],
        dataFormats= ['csv']*8, units= [ None ] *6 + [ 'pb' ] * 2)

T2ttoff = dataset.addTxName('T2ttoff')
# T2ttoff.checked ="VM"
T2ttoff.constraint ="[[['b','W']],[['b','W']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.massConstraint = [['80 <= dm <= 169.0'], ['80 <= dm <= 169.0']]
T2ttoff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2ttoff_1 = T2ttoff.addMassPlane(2*[[x, y]])
T2ttoff_1.figure = 'Fig. 01a'
T2ttoff_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-12/figaux_01b.png'
T2ttoff_1.dataUrl = 'https://www.hepdata.net/record/ins1793461?version=2&table=stop_xSecUpperLimit_obs'
T2ttoff_1.setSources(dataLabels= ['obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1', 'upperLimits', 'expectedUpperLimits' ],
        dataFiles= [ 'orig/stopoff_obs.csv', 'orig/stopoff_obs_up.csv', 'orig/stopoff_obs_down.csv', 'orig/stopoff_exp.csv', 'orig/stopoff_exp_up.csv', 'orig/stopoff_exp_down.csv', 'orig/stop_xSecUpperLimit_obs.csv', 'orig/stop_xSecUpperLimit_exp.csv' ],
        dataFormats= ['csv']*8, units= [ None ] *6 + [ 'pb' ] * 2)

T2bbffff = dataset.addTxName('T2bbffff')
# T2bbffff.checked ="VM"
T2bbffff.constraint = "9./4.*[[['b', 'jet','jet']],[['b', 'jet','jet']]]"
T2bbffff.conditionDescription ="None"
T2bbffff.condition ="None"
T2bbffff.massConstraint = [['dm < 80'], ['dm < 80']]
T2bbffff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2bbffff_1 = T2bbffff.addMassPlane(2*[[x, y]])
T2bbffff_1.figure = 'Fig. 01a'
T2bbffff_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-12/figaux_01b.png'
T2bbffff_1.dataUrl = 'https://www.hepdata.net/record/ins1793461?version=2&table=stop_xSecUpperLimit_obs'
T2bbffff_1.setSources(dataLabels= ['obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1', 'upperLimits', 'expectedUpperLimits' ],
        dataFiles= [ 'orig/stop_obs.csv', 'orig/stop_obs_up.csv', 'orig/stop_obs_down.csv', 'orig/stop_exp.csv', 'orig/stop_exp_up.csv', 'orig/stop_exp_down.csv', 'orig/stop_xSecUpperLimit_obs.csv', 'orig/stop_xSecUpperLimit_exp.csv' ],
        dataFormats= ['csv']*8, units= [ None ] *6 + [ 'pb' ] * 2)


databaseCreator.create()
