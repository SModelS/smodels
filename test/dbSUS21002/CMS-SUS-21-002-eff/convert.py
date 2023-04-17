#!/usr/bin/env python3

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse
import types
import math
import uproot

argparser = argparse.ArgumentParser(description =
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath',
help = 'path to the package smodels_utils',\
type = str)
argparser.add_argument ('-smodelsPath', '--smodelsPath',
help = 'path to the package smodels_utils',\
type = str)
argparser.add_argument ('-no', '--noUpdate',
help = 'do not update the lastUpdate field.',\
action= "store_true" )
args = argparser.parse_args()

if args.noUpdate:
    os.environ["SMODELS_NOUPDATE"]="1"

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
from smodels_utils.dataPreparation.datasetCreation import aggregateDataSets, \
         createAggregationOrder
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z
from smodels_utils.dataPreparation import dataHandlerObjects

dataHandlerObjects.allowTrimming=False



#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-21-002')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/'
info.sqrts = '13.0*TeV'
info.lumi = 137
info.prettyName = 'Hadronic EWK searches'
info.publication = "Accepted for publication in Phys. Lett. B"
info.private = False
info.arxiv = 'https://arxiv.org/abs/2205.09597'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.comment = 'Root files for efficiency maps were given by the CMS collaboration, not present in HEPData. Efficiency map files have efficiency 0.0 at 7 Gev LSP Mass for all SR - removed manually from text files(take care when rerunning convert.py)'
#info.datasetOrder = " 'b_veto_SR', 'WH_SR', 'W_SR' ,'H_SR'"

aggregates = None
#aggregates = [[1,2,3,4,5,6,7,8,9], [10,11,12,13,14,15,16,17], [18,19,20,21,22,23,24,25,26],[27,28,29,30,31,32,33,34,35]]

datasets = []

b_veto_nobs = [82, 48, 24, 9, 8, 6, 9, 6, 3]
b_veto_nbg = [88.1, 48.2, 25.4, 15.5, 12.3, 7.9, 9.6, 6.2, 0.9]
b_veto_bgerr = [9.2,6.7,4.4,3.7,3.1,2.1,2.8,2.2,0.7]

WH_SR_nobs = [32, 13, 8, 8, 4, 3, 1, 1]
WH_SR_nbg =  [34.9, 16.2, 5.2, 3.1, 2.9, 1.29, 1.19, 1.02]
WH_SR_bgerr = [5.4, 3.3, 1.3, 1.3, 1.2, 0.62, 1.1, 1.1]

W_SR_nobs = [680, 312, 175, 95, 63, 23, 28, 17, 4]
W_SR_nbg = [749.4, 332.4, 176.4, 99.4, 60.7, 34.7, 26.9, 18.8, 3.5]
W_SR_bgerr = [46.1, 32.5, 13.5, 9.9, 6.9, 4.7, 3.6, 3.2, 1.7]

H_SR_nobs = [1212, 563, 282, 160, 115, 60, 65, 39, 8]
H_SR_nbg = [1244, 555, 311, 156, 93, 77.1, 55.8, 30.2, 7.5]
H_SR_bgerr = [61, 27, 20, 14, 10, 8.0, 7.1, 4.8, 2.5]


#+++++++++dataset block++++++++++++++    
#+++++++++b_veto SR++++++++++++++

order = []

for i in range(len(b_veto_nobs)):
    dsname = f'b_veto_SR{i}'
    order.append ( f"'{dsname}'" )
    dataset = DataSetInput( dsname )
    dataset.setInfo(dataType = 'efficiencyMap', dataId = dsname,
                    observedN = b_veto_nobs[i], expectedBG = b_veto_nbg[i], bgError = b_veto_bgerr[i], comment = dsname )

    #+++++++++txName block+++++++++++++++++

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

    
    #+++++++++mass plane block+++++++++++++++++
    TChiWW_1 = TChiWW.addMassPlane(2*[[x,y]])
    TChiWW_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bVeto_TChiWW.root', 'root', objectName = 'AccXEff_bVeto;1', index = i )
    TChiWW_1.addSource( 'obsExclusion', 'orig/TChiWW_obsExclusion.csv', 'csv' )  
    TChiWW_1.addSource( 'obsExclusionP1', 'orig/TChiWW_obsExclusionP1.csv', 'csv' )
    TChiWW_1.addSource( 'obsExclusionM1', 'orig/TChiWW_obsExclusionM1.csv', 'csv' ) 
    TChiWW_1.addSource( 'expExclusion', 'orig/TChiWW_expExclusion.csv', 'csv')
    TChiWW_1.addSource( 'expExclusionP1', 'orig/TChiWW_expExclusionP1.csv', 'csv')
    TChiWW_1.addSource( 'expExclusionM1', 'orig/TChiWW_expExclusionM1.csv', 'csv' )
    TChiWW_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_013.root"
    TChiWW_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_013.png"
    

    TChiWZ_1 = TChiWZ.addMassPlane(2*[[x,y]])
    TChiWZ_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bVeto_TChiWZ.root', 'root', objectName = 'AccXEff_bVeto;1', index = i )
    TChiWZ_1.addSource( 'obsExclusion', 'orig/TChiWZ_obsExclusion.csv', 'csv' )  
    TChiWZ_1.addSource( 'obsExclusionP1', 'orig/TChiWZ_obsExclusionP1.csv', 'csv' )
    TChiWZ_1.addSource( 'obsExclusionM1', 'orig/TChiWZ_obsExclusionM1.csv', 'csv' ) 
    TChiWZ_1.addSource( 'expExclusion', 'orig/TChiWZ_expExclusion.csv', 'csv')
    TChiWZ_1.addSource( 'expExclusionP1', 'orig/TChiWZ_expExclusionP1.csv', 'csv')
    TChiWZ_1.addSource( 'expExclusionM1', 'orig/TChiWZ_expExclusionM1.csv', 'csv' )
    TChiWZ_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_014-a.root"
    TChiWZ_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_014-a.png"

    
    TChiWH_1 = TChiWH.addMassPlane(2*[[x,y]])
    TChiWH_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bVeto_TChiWH.root', 'root', objectName = 'AccXEff_bVeto;1', index = i )
    TChiWH_1.addSource( 'obsExclusion', 'orig/TChiWH_obsExclusion.csv', 'csv' )
    TChiWH_1.addSource( 'obsExclusionP1', 'orig/TChiWH_obsExclusionP1.csv', 'csv' )
    TChiWH_1.addSource( 'obsExclusionM1', 'orig/TChiWH_obsExclusionM1.csv', 'csv' )
    TChiWH_1.addSource( 'expExclusion', 'orig/TChiWH_expExclusion.csv', 'csv')
    TChiWH_1.addSource( 'expExclusionP1', 'orig/TChiWH_expExclusionP1.csv', 'csv')
    TChiWH_1.addSource( 'expExclusionM1', 'orig/TChiWH_expExclusionM1.csv', 'csv' )
    TChiWH_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_015-c.root"
    TChiWH_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_015-c.png"
    

    datasets.append ( dataset )

#+++++++++++WH_SR+++++++++++++++++

for i in range(len(WH_SR_nobs)):
    dsname = f'WH_SR{i}'
    dataset = DataSetInput( dsname )
    order.append ( f"'{dsname}'" )
    dataset.setInfo(dataType = 'efficiencyMap', dataId = dsname,
                    observedN = WH_SR_nobs[i], expectedBG = WH_SR_nbg[i],bgError=WH_SR_bgerr[i],comment=dsname)

    #+++++++++txName block+++++++++++++++++

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

    #+++++++++mass plane block+++++++++++++++++
    
    TChiWW_1 = TChiWW.addMassPlane(2*[[x,y]])
    TChiWW_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bTag_TChiWW.root', 'root', objectName = 'AccXEff_WHSR;1', index = i )
    TChiWW_1.addSource( 'obsExclusion', 'orig/TChiWW_obsExclusion.csv', 'csv' )
    TChiWW_1.addSource( 'obsExclusionP1', 'orig/TChiWW_obsExclusionP1.csv', 'csv' )
    TChiWW_1.addSource( 'obsExclusionM1', 'orig/TChiWW_obsExclusionM1.csv', 'csv' )
    TChiWW_1.addSource( 'expExclusion', 'orig/TChiWW_expExclusion.csv', 'csv')
    TChiWW_1.addSource( 'expExclusionP1', 'orig/TChiWW_expExclusionP1.csv', 'csv')
    TChiWW_1.addSource( 'expExclusionM1', 'orig/TChiWW_expExclusionM1.csv', 'csv' )
    TChiWW_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_013.root"
    TChiWW_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_013.png"


    TChiWZ_1 = TChiWZ.addMassPlane(2*[[x,y]])
    TChiWZ_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bTag_TChiWZ.root', 'root', objectName = 'AccXEff_WHSR;1', index = i )
    TChiWZ_1.addSource( 'obsExclusion', 'orig/TChiWZ_obsExclusion.csv', 'csv' )  
    TChiWZ_1.addSource( 'obsExclusionP1', 'orig/TChiWZ_obsExclusionP1.csv', 'csv' )
    TChiWZ_1.addSource( 'obsExclusionM1', 'orig/TChiWZ_obsExclusionM1.csv', 'csv' ) 
    TChiWZ_1.addSource( 'expExclusion', 'orig/TChiWZ_expExclusion.csv', 'csv')
    TChiWZ_1.addSource( 'expExclusionP1', 'orig/TChiWZ_expExclusionP1.csv', 'csv')
    TChiWZ_1.addSource( 'expExclusionM1', 'orig/TChiWZ_expExclusionM1.csv', 'csv' )
    TChiWZ_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_014-b.root"
    TChiWZ_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_014-b.png"


    TChiWH_1 = TChiWH.addMassPlane(2*[[x,y]])
    TChiWH_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bTag_TChiWH.root', 'root', objectName = 'AccXEff_WHSR;1', index = i )
    TChiWH_1.addSource( 'obsExclusion', 'orig/TChiWH_obsExclusion.csv', 'csv' )  
    TChiWH_1.addSource( 'obsExclusionP1', 'orig/TChiWH_obsExclusionP1.csv', 'csv' )
    TChiWH_1.addSource( 'obsExclusionM1', 'orig/TChiWH_obsExclusionM1.csv', 'csv' ) 
    TChiWH_1.addSource( 'expExclusion', 'orig/TChiWH_expExclusion.csv', 'csv')
    TChiWH_1.addSource( 'expExclusionP1', 'orig/TChiWH_expExclusionP1.csv', 'csv')
    TChiWH_1.addSource( 'expExclusionM1', 'orig/TChiWH_expExclusionM1.csv', 'csv' )
    TChiWH_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_015-a.root"
    TChiWH_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_015-a.png"

    datasets.append ( dataset )

#+++++++++++++++W_SR+++++++++++++++++
for i in range(len(W_SR_nobs)):
    dsname = f'W_SR{i}'
    dataset = DataSetInput( dsname )
    order.append ( f"'{dsname}'" )
    dataset.setInfo(dataType = 'efficiencyMap', dataId = dsname,
                    observedN = W_SR_nobs[i], expectedBG = W_SR_nbg[i], bgError = W_SR_bgerr[i], comment = dsname )

    #+++++++++txName block+++++++++++++++++
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

    #+++++++++mass plane block+++++++++++++++++
    TChiWW_1 = TChiWW.addMassPlane(2*[[x,y]])
    TChiWW_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bTag_TChiWW.root', 'root', objectName = 'AccXEff_WSR;1', index = i )
    TChiWW_1.addSource( 'obsExclusion', 'orig/TChiWW_obsExclusion.csv', 'csv' )
    TChiWW_1.addSource( 'obsExclusionP1', 'orig/TChiWW_obsExclusionP1.csv', 'csv' )
    TChiWW_1.addSource( 'obsExclusionM1', 'orig/TChiWW_obsExclusionM1.csv', 'csv' )
    TChiWW_1.addSource( 'expExclusion', 'orig/TChiWW_expExclusion.csv', 'csv')
    TChiWW_1.addSource( 'expExclusionP1', 'orig/TChiWW_expExclusionP1.csv', 'csv')
    TChiWW_1.addSource( 'expExclusionM1', 'orig/TChiWW_expExclusionM1.csv', 'csv' )
    TChiWW_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_013.root"
    TChiWW_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_013.png"

    TChiWZ_1 = TChiWZ.addMassPlane(2*[[x,y]])
    TChiWZ_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bTag_TChiWZ.root', 'root', objectName = 'AccXEff_WSR;1', index = i )
    TChiWZ_1.addSource( 'obsExclusion', 'orig/TChiWZ_obsExclusion.csv', 'csv' )  
    TChiWZ_1.addSource( 'obsExclusionP1', 'orig/TChiWZ_obsExclusionP1.csv', 'csv' )
    TChiWZ_1.addSource( 'obsExclusionM1', 'orig/TChiWZ_obsExclusionM1.csv', 'csv' ) 
    TChiWZ_1.addSource( 'expExclusion', 'orig/TChiWZ_expExclusion.csv', 'csv')
    TChiWZ_1.addSource( 'expExclusionP1', 'orig/TChiWZ_expExclusionP1.csv', 'csv')
    TChiWZ_1.addSource( 'expExclusionM1', 'orig/TChiWZ_expExclusionM1.csv', 'csv' )
    TChiWZ_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_014-c.root"
    TChiWZ_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_014-c.png"
    
    TChiWH_1 = TChiWH.addMassPlane(2*[[x,y]])
    TChiWH_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bTag_TChiWH.root', 'root', objectName = 'AccXEff_WSR;1', index = i )
    TChiWH_1.addSource( 'obsExclusion', 'orig/TChiWH_obsExclusion.csv', 'csv' )  
    TChiWH_1.addSource( 'obsExclusionP1', 'orig/TChiWH_obsExclusionP1.csv', 'csv' )
    TChiWH_1.addSource( 'obsExclusionM1', 'orig/TChiWH_obsExclusionM1.csv', 'csv' ) 
    TChiWH_1.addSource( 'expExclusion', 'orig/TChiWH_expExclusion.csv', 'csv')
    TChiWH_1.addSource( 'expExclusionP1', 'orig/TChiWH_expExclusionP1.csv', 'csv')
    TChiWH_1.addSource( 'expExclusionM1', 'orig/TChiWH_expExclusionM1.csv', 'csv' )
    TChiWH_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_015-b.root"
    TChiWH_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_015-b.png"
    
    datasets.append ( dataset )


#++++++++++++++++H_SR++++++++++++++++++

for i in range(len(H_SR_nobs)):
    dsname = f'H_SR{i}'
    dataset = DataSetInput( dsname )
    order.append ( f"'{dsname}'" )
    dataset.setInfo(dataType = 'efficiencyMap', dataId = dsname,
                    observedN = H_SR_nobs[i], expectedBG = H_SR_nbg[i], bgError = H_SR_bgerr[i], comment = dsname )

    #+++++++++txName block+++++++++++++++++
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

    #+++++++++mass plane block+++++++++++++++++
    TChiWW_1 = TChiWW.addMassPlane(2*[[x,y]])
    TChiWW_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bTag_TChiWW.root', 'root', objectName = 'AccXEff_HSR;1', index = i )
    #TChiWW_1.efficiencyMap._unit = "0."
    TChiWW_1.addSource( 'obsExclusion', 'orig/TChiWW_obsExclusion.csv', 'csv' )
    TChiWW_1.addSource( 'obsExclusionP1', 'orig/TChiWW_obsExclusionP1.csv', 'csv' )
    TChiWW_1.addSource( 'obsExclusionM1', 'orig/TChiWW_obsExclusionM1.csv', 'csv' )
    TChiWW_1.addSource( 'expExclusion', 'orig/TChiWW_expExclusion.csv', 'csv')
    TChiWW_1.addSource( 'expExclusionP1', 'orig/TChiWW_expExclusionP1.csv', 'csv')
    TChiWW_1.addSource( 'expExclusionM1', 'orig/TChiWW_expExclusionM1.csv', 'csv' )
    TChiWW_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_013.root"
    TChiWW_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_013.png"


    TChiWZ_1 = TChiWZ.addMassPlane(2*[[x,y]])
    TChiWZ_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bTag_TChiWZ.root', 'root', objectName = 'AccXEff_HSR;1', index = i )
    TChiWZ_1.addSource( 'obsExclusion', 'orig/TChiWZ_obsExclusion.csv', 'csv' )  
    TChiWZ_1.addSource( 'obsExclusionP1', 'orig/TChiWZ_obsExclusionP1.csv', 'csv' )
    TChiWZ_1.addSource( 'obsExclusionM1', 'orig/TChiWZ_obsExclusionM1.csv', 'csv' ) 
    TChiWZ_1.addSource( 'expExclusion', 'orig/TChiWZ_expExclusion.csv', 'csv')
    TChiWZ_1.addSource( 'expExclusionP1', 'orig/TChiWZ_expExclusionP1.csv', 'csv')
    TChiWZ_1.addSource( 'expExclusionM1', 'orig/TChiWZ_expExclusionM1.csv', 'csv' )
    TChiWZ_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_014-d.root"
    TChiWZ_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_014-d.png"

    TChiWH_1 = TChiWH.addMassPlane(2*[[x,y]])
    TChiWH_1.addSource( 'efficiencyMap', f'orig/AllBinAccXEff_bTag_TChiWH.root', 'root', objectName = 'AccXEff_HSR;1', index = i )
    TChiWH_1.addSource( 'obsExclusion', 'orig/TChiWH_obsExclusion.csv', 'csv' )  
    TChiWH_1.addSource( 'obsExclusionP1', 'orig/TChiWH_obsExclusionP1.csv', 'csv' )
    TChiWH_1.addSource( 'obsExclusionM1', 'orig/TChiWH_obsExclusionM1.csv', 'csv' ) 
    TChiWH_1.addSource( 'expExclusion', 'orig/TChiWH_expExclusion.csv', 'csv')
    TChiWH_1.addSource( 'expExclusionP1', 'orig/TChiWH_expExclusionP1.csv', 'csv')
    TChiWH_1.addSource( 'expExclusionM1', 'orig/TChiWH_expExclusionM1.csv', 'csv' )
    TChiWH_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_015-c.root"
    TChiWH_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-21-002/CMS-SUS-21-002_Figure-aux_015-c.png"

    datasets.append ( dataset )

info.createCovarianceMatrix ( "orig/Covariance_matrix.root", histoname = "overall_total_covar;1", histoIsCorrelations=False, addOrder = False, max_datasets = None, datasets = datasets, aggregate = aggregates )
info.datasetOrder = ", ".join ( order )

databaseCreator.create() 
