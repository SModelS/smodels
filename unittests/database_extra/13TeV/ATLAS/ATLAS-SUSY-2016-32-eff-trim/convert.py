#!/usr/bin/env python3

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
    type = str, default = '~/smodels-utils' )
argparser.add_argument ('-smodelsPath', '--smodelsPath',
    help = 'path to the package smodels_utils',\
    type = str, default = '~/smodels' )
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

utilsPath = args.utilsPath
sys.path.append(os.path.abspath(os.path.expanduser(utilsPath)))
sys.path.append(os.path.abspath(os.path.expanduser(args.smodelsPath)))


from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z, w

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2016-32')
info.url = 'http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-32/index.html'
info.sqrts = 13
info.lumi = 31.6
info.Leff_inner = 1e-3
info.Leff_outer = 12.0
info.prettyName ='hscp search'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1902.01636'
info.contact ='Andre Lessa <lessa.a.p@gmail.com>; Jan Heisig <heisig@physik.rwth-aachen.de>; SModelS'
info.publication ='https://arxiv.org/pdf/1902.01636.pdf'
info.comment ='Search for long-lived charged particles implemented according to the ATLAS note in http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-32/hepdata_info.pdf. For the topologies with mixed MET-HSCP branches, the MET branch is irrelevant and wildcards are used. The topologies THSCPM5, THSCPM6 and THSCPM7 do not include the width information.'
info.supersedes =''
particlesFile = os.path.abspath('../../../databaseParticles.py')

#++Define list of datasets++
datasetNames = ['SR2FULL_150','SR2FULL_350','SR2FULL_575','SR2FULL_800', 'SR1FULL_175', 'SR1FULL_375', 'SR1FULL_600', 'SR1FULL_825']
#Values from arxiv:1902.01636 (Table 5).
observedNs = [0,0,0,0,227,16,1,0]
expectedBGs = [1.5,0.06,0.007,0.0017,240.,17.,2.2,0.48]
bgErrors = [0.3,0.01,0.002,0.0009,20.,2.,0.2,0.07]
#SR Upper limits from arxiv:1902.01636 (Table 5).
obsUpperLimits = ['0.09*fb','0.08*fb','0.08*fb','0.08*fb','1.26*fb','0.24*fb','0.10*fb','0.08*fb']
#Assume expected UL = observed UL (otherwise they will be overwritten)
expUpperLimits = ['0.09*fb','0.08*fb','0.08*fb','0.08*fb','1.26*fb','0.24*fb','0.10*fb','0.08*fb']



#Txname data:
txnames = {
    'THSCPM1b' : {'constraint' : "[[],[]]", 'finalState' : ['HSCP','HSCP'], 'massPlane' : [[(x,y)]]*2,
                        'coordinates' : {x: 0, y:2, 'value' : 3}},
    'THSCPM2b' : {'constraint' : "[['*'],[]]", 'finalState' : ['MET','HSCP'], 'massPlane' : [['*'],[(x,y)]],
                        'coordinates' : {x: 0, y:2, 'value' : 3}},
    'THSCPM3' : {'constraint' : "[[['*']],[['*']]]", 'finalState' : ['HSCP','HSCP'], 'massPlane' : [[x,(y,w)]]*2,
                        'coordinates' : {x: 1, y:0, w: 3, 'value' : 4}},
    'THSCPM4' : {'constraint' : "[[*],[['*']]]", 'finalState' : ['MET','HSCP'], 'massPlane' : [['*'],[x,(y,w)]],
                        'coordinates' : {x: 1, y:0, w: 3, 'value' : 4}},
    'THSCPM5' : {'constraint' :"[[['*'],['*']],[['*'],['*']]]", 'finalState' : ['HSCP','HSCP'], 'massPlane' : [[x,y,z]]*2,
                        'coordinates' : {x: 2, y:1, z: 0, 'value' : 5}},
    'THSCPM6' : {'constraint' : "[[*],[['*'],[ '*']]]", 'finalState' : ['MET','HSCP'], 'massPlane' : [['*'],[x,y,z]],
                        'coordinates' : {x: 2, y:1, z: 0, 'value' : 5}},
   'THSCPM7' : {'constraint' :  "[[['*'],['*']],[['*']]]", 'finalState' : ['HSCP','HSCP'], 'massPlane' : [[x,y,z],[x,z]],
                        'coordinates' : {x: 2, y:1, z: 0, 'value' : 5}},
    'THSCPM8' : {'constraint' : "[[['*','*']],[['*','*']]]", 'finalState' : ['HSCP','HSCP'], 'massPlane' : [[x,(y,w)]]*2,
                        'coordinates' : {x: 1, y:0, w: 3, 'value' : 4}},
   'THSCPM9' : {'constraint' : "[[*],[['*','*']]]", 'finalState' : ['MET','HSCP'], 'massPlane' : [['*'],[x,(y,w)]],
                        'coordinates' : {x: 1, y:0, w: 3, 'value' : 4}},
   'THSCPM10' : {'constraint' : "[[],[['*']]]", 'finalState' : ['HSCP','HSCP'], 'massPlane' : [[(y,w)],[x,(y,w)]],
                        'coordinates' : {x: 1, y:0, w: 3, 'value' : 4}},
   'THSCPM11' : {'constraint' : "[[],[['*','*']]]", 'finalState' : ['HSCP','HSCP'], 'massPlane' : [[(y,w)],[x,(y,w)]],
                        'coordinates' : {x: 1, y:0, w: 3, 'value' : 4}}

}


for idataset,name in enumerate(datasetNames):
#+++++++ dataset block ++++++++++++++
    dataset = DataSetInput(name)
    dataset.setInfo(dataType = 'efficiencyMap', dataId = name,
            observedN=observedNs[idataset], expectedBG=expectedBGs[idataset], bgError=bgErrors[idataset],
            upperLimit = obsUpperLimits[idataset], expectedUpperLimit = expUpperLimits[idataset])

    #+++++++ txnames ++++++++++++++++++++
    for tx in txnames:
        #Skip txnames with a single HSCP for the SR requiring 2 candidates
        if 'SR2FULL' in name and 'MET' in txnames[tx]['finalState']:
            continue
        #+++++++ next txName block ++++++++++++++
        Txname = dataset.addTxName(tx)
        Txname.checked =''
        Txname.constraint = txnames[tx]['constraint']
        Txname.condition =None
        Txname.finalState = txnames[tx]['finalState']
        Txname.massConstraints = None
        Txname.dataUrl = None
        Txname.source = 'SModelS'
        Txname.setParticlesFromFile(particlesFile)
        #+++++++ next mass plane block ++++++++++++++
        plane = Txname.addMassPlane(txnames[tx]['massPlane'])
        if not (tx in ['THSCPM1b','THSCPM2b','THSCPM3','THSCPM4','THSCPM8','THSCPM9','THSCPM10','THSCPM11']):
            dataFile = 'orig/%s_eff_mutrig_%s_noWidth.dat' %(tx,name) #Use maps without width info
        elif tx in ['THSCPM10','THSCPM11']:
            dataFile = 'orig/%s_eff_mutrig_reduced_th_%s.dat' %(tx,name)
        else:
            dataFile = 'orig/%s_eff_mutrig_%s_trimmedwidth.dat' %(tx,name)
        plane.addSource(dataLabel = 'efficiencyMap',
                    dataFile = dataFile,
                    coordinateMap = txnames[tx]['coordinates'],
                    dataFormat = 'txt')
        #++++++ exclusion (only for THSCPM1b) ++++
        if tx == 'THSCPM1b':
            plane.addSource(dataLabel='obsExclusion',dataFile='orig/Stau_ExclusionObs.csv',
                            coordinateMap = {x : 0, y: 1, 'value' : None}, dataFormat = 'csv')

databaseCreator.create()
