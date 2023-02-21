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
info = MetaInfoInput('ATLAS-SUSY-2018-23')
# info.comment = ''
info.sqrts = '13.0'
info.private = False
info.lumi = '139.'
info.publication = "https://link.springer.com/article/10.1007/JHEP10(2020)005"
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-23/'
info.prettyName = 'EWK WH(2 photons)'
info.implementedBy = 'WW'

# https://www.hepdata.net/record/ins1289225

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

TChiWH = dataset.addTxName('TChiWH')
# TChiWH.checked ="VM"
TChiWH.constraint ="[[['W']],[['higgs']]]"
TChiWH.conditionDescription ="None"
TChiWH.condition ="None"
TChiWH.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
TChiWH = TChiWH.addMassPlane(2*[[x, y ]])
TChiWH.figure = 'Fig. 10'
TChiWH.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-23/fig_10.png'
TChiWH.dataUrl = 'https://doi.org/10.17182/hepdata.97041.v1/t32'
TChiWH.setSources(dataLabels= ['obsExclusion', 'expExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusionP1', 'expExclusionM1', 'upperLimits', 'expectedUpperLimits' ],
        dataFiles= [ 'orig/Figure10Observed.csv', 'orig/Figure10Expected.csv', 'orig/Figure10+1Obssigma.csv', 'orig/Figure10-1Obssigma.csv', 'orig/Figure10+1Expsigma.csv', 'orig/Figure10-1Expsigma.csv', 'orig/table7.csv', 'orig/table7.csv' ],
                 dataFormats= ['csv'] * 8, units= [ None ] *6 + [ 'pb' ] * 2,
                 indices = [ None ] * 6 + [ 3, 2 ] )

TChiHH = dataset.addTxName('TChiHH')
# TChiHH.checked ="VM"
TChiHH.constraint ="[[['higgs']],[['higgs']]]"
TChiHH.conditionDescription ="None"
TChiHH.condition ="None"
TChiHH.source = 'ATLAS'

#+++++++ next mass plane block ++++++++++++++
TChiHH1 = TChiHH.addMassPlane(2*[[x, 1. ]])
TChiHH1.figure = 'Fig. 9'
TChiHH1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-23/fig_09.png'
TChiHH1.dataUrl = 'https://www.hepdata.net/record/ins1792399'
 
TChiHH1.setSources(dataLabels= ['obsExclusion', 'upperLimits', 'expectedUpperLimits' ],
        dataFiles= [ 'orig/excl.csv', 'orig/table8.csv', 'orig/table8.csv' ],
                 dataFormats= ['csv'] * 3, units= [ None ] *1 + [ 'pb' ] * 2,
                 indices = [ None ] * 1 + [ 1, 2 ], coordinates = [ {x: 0, 'value': None} ] *1 + [ None ] * 2 )

## feed mass plane twice
TChiHH0 = TChiHH.addMassPlane(2*[[x, 0. ]])
TChiHH0.figure = 'Fig. 9'
TChiHH0.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-23/fig_09.png'
TChiHH0.dataUrl = 'https://www.hepdata.net/record/ins1792399'
 
TChiHH0.setSources(dataLabels= ['obsExclusion', 'upperLimits', 'expectedUpperLimits' ],
        dataFiles= [ 'orig/excl.csv', 'orig/table8.csv', 'orig/table8.csv' ],
                 dataFormats= ['csv'] * 3, units= [ None ] *1 + [ 'pb' ] * 2,
                 indices = [ None ] * 1 + [ 1, 2 ], coordinates = [ {x: 0, 'value': None} ] *1 + [ None ] * 2 )

databaseCreator.create()
