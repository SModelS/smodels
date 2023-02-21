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
info             = MetaInfoInput('CMS-SUS-21-002')
info.url         = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/'
info.sqrts       = 13
info.lumi        = 137
info.prettyName  = 'Hadronic EWK searches'
info.private     = False 
info.arxiv       = 'https://arxiv.org/abs/2205.09597'
info.contact     = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Accepted for publication in Phys. Lett. B'
# info.comment     = 'lumi given as 129 - 137 / fb'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#====================  T2tt  ==============================  

#+++++++txName block++++++++++++++++++++

TChiWW=dataset.addTxName('TChiWW')
TChiWW.checked=''
TChiWW.constraint="[[['W+']],[['W-']]]"
TChiWW.condition=None
TChiWW.conditionDescription = None
TChiWW.source="CMS"

TChiWZ=dataset.addTxName('TChiWZ')
TChiWZ.checked=''
TChiWZ.constraint="[[['W']],[['Z']]]"
TChiWZ.condition=None
TChiWZ.conditionDescription = None
TChiWZ.source="CMS"


TChiWH=dataset.addTxName('TChiWH')
TChiWH.checked=''
TChiWH.constraint="[[['W']],[['higgs']]]"
TChiWH.condition=None
TChiWH.conditionDescription = None
TChiWH.source="CMS"


#++++++mass plane block+++++++++

TChiWW_1 = TChiWW.addMassPlane(2*[[x,y]])
TChiWW_1.figure='Fig. 4a'
TChiWW_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure_004-a.png'
TChiWW_1.dataUrl='https://doi.org/10.17182/hepdata.127766.v1/t11'
TChiWW_1.setSources(dataLabels=['upperLimits',  'obsExclusion', 'obsExclusionP1', 'obsExclusionM1','expExclusion', 'expExclusionP1', 'expExclusionM1' ],
                    dataFiles=[ 'orig/Figure4aObserved.csv', 'orig/TChiWW_obsExclusion.csv', 'orig/TChiWW_obsExclusionP1.csv', 'orig/TChiWW_obsExclusionM1.csv',
                                 'orig/TChiWW_expExclusion.csv', 'orig/TChiWW_expExclusionP1.csv', 'orig/TChiWW_expExclusionM1.csv' ],
                    dataFormats=['csv']*7,
                    indices=[[2,5,6],None,None,None,None,None,None],
                    units=['pb']*1+[None]*6 )

TChiWZ_1 = TChiWZ.addMassPlane(2*[[x,y]])
TChiWZ_1.figure='Fig. 4b'
TChiWZ_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure_004-b.png'
TChiWZ_1.dataUrl='https://doi.org/10.17182/hepdata.127766.v1/t14'

TChiWZ_1.setSources(dataLabels=['upperLimits',  'obsExclusion', 'obsExclusionP1', 'obsExclusionM1','expExclusion', 'expExclusionP1', 'expExclusionM1' ],
                    dataFiles=[ 'orig/Figure4bObserved.csv',  'orig/TChiWZ_obsExclusion.csv', 'orig/TChiWZ_obsExclusionP1.csv', 'orig/TChiWZ_obsExclusionM1.csv',
                                 'orig/TChiWZ_expExclusion.csv', 'orig/TChiWZ_expExclusionP1.csv', 'orig/TChiWZ_expExclusionM1.csv' ],
                    dataFormats=['csv']*7, indices=[[2,5,6],None,None,None,None,None,None], units=['pb']*1+[None]*6 )

TChiWZ_1.addSource( 'expectedUpperLimits',  'orig/CMS-SUS-21-002_Figure_004-b.root', 'root', objectName = 'glimexp;1', unit = 'pb')

#objectNames = ['glimobs;1']+[None]*6

TChiWH_1 = TChiWH.addMassPlane(2*[[x,y]])
TChiWH_1.figure='Fig. 4c'
TChiWH_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure_004-c.png'
TChiWH_1.dataUrl='https://doi.org/10.17182/hepdata.127766.v1/t17'
TChiWH_1.setSources(dataLabels=['upperLimits', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1','expExclusion', 'expExclusionP1', 'expExclusionM1' ],
                    dataFiles=[ 'orig/Figure4cObserved.csv', 'orig/CMS-SUS-21-002_Figure_004-c.root', 'orig/TChiWH_obsExclusion.csv', 'orig/TChiWH_obsExclusionP1.csv', 'orig/TChiWH_obsExclusionM1.csv',
                                 'orig/TChiWH_expExclusion.csv', 'orig/TChiWH_expExclusionP1.csv', 'orig/TChiWH_expExclusionM1.csv' ],
                    dataFormats= ['csv']+['root']+['csv']*6, objectNames = [None, 'glimexp;1']+[None]*6,
                    indices=[[2,5,6]]+[None]*7,
                    units=['pb']*2+[None]*6 )


databaseCreator.create()
