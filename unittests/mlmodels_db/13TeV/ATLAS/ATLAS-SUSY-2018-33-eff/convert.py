#!/usr/bin/env python3

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.
"""
import sys
import os
import argparse
import types
import numpy as np

argparser = argparse.ArgumentParser(description =
    'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath',
    help = 'path to the package smodels_utils',\
    type = str, default = '~/SModelS3/smodels-utils' )
argparser.add_argument ('-smodelsPath', '--smodelsPath',
    help = 'path to the package smodels_utils',\
    type = str, default = '~/SModelS3/smodels' )
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
from smodels_utils.dataPreparation.massPlaneObjects import x, y,z
from smodels.base.physicsUnits import GeV,TeV

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2018-33')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-33/'
info.sqrts = 13
info.lumi = 136.0
info.prettyName = 'DV with >= 1 muon'
info.private = False
info.arxiv = 'https://arxiv.org/abs/2003.11956'
info.contact ='atlas-phys-susy-conveners@cern.ch'
info.type = 'displaced'
info.publication ='Phys. Rev. D. 102 (2020) 032006'
info.publicationDOI='https://doi.org/10.1103/PhysRevD.102.032006'
info.comment ="Efficiencies extracted from patchset files in orig/likelihood folder. The bkgonly.json files, however, are hand-made by ATLAS and correspond to a simplified llhd. Moreover, there's one bkgonly.json file per SR, so no combination. Provided LLHD leads to over-exclusion; need to use sig95 values from Table 6 to match official exclusion line."
# particlesFile = os.path.abspath('../../../databaseParticles.py')
#taking out the json files as only one json per SR; also ATLAS does not combine
#info.datasetOrder = '"SR_MET", "SR_MU"'
#info.jsonFiles = '{"SRMET_bkgonly.json":["SR_MET"], "SRMU_bkgonly.json":["SR_MU"] }'


datasets = []
#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('SR_MET')
dataset.setInfo( dataType = 'efficiencyMap', dataId = 'SR_MET',
                    observedN = 0, expectedBG = 0.43, bgError=0.16*np.sqrt(2), 
                    upperLimit = '0.0228*fb', expectedUpperLimit = '0.0228*fb')

#+++++++ txnames ++++++++++++++++++++
lifetime = np.array([0.01, 0.1, 1.0])
hbar = 6.582*10**(-16)
width = hbar/lifetime

#+++++++ next txName block ++++++++++++++
TRPV1jmu = dataset.addTxName('TRPV1jmu')
#TRPV1jmu.validationTarball = "TRPV1jmu.tar.gz"
TRPV1jmu.checked =''
TRPV1jmu.constraint = "{(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > jet,mu+), (anyBSM(2) > jet,mu-)}"
TRPV1jmu.condition =None
TRPV1jmu.source = 'ATLAS'
TRPV1jmu.finalState = None
TRPV1jmu.intermediateState = None
TRPV1jmu.dataMap = {0:(1,'mass',GeV), 1:(1,'totalwidth',GeV), 2:(2,'mass',GeV), 3:(2,'totalwidth',GeV)}

#+++++++ next mass plane block ++++++++++++++
plane_1 = TRPV1jmu.addMassPlane({0 : 'x', 1 : 'y', 2 : 'x', 3 : 'y'})
plane_1.figureUrl = None # 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-33/fig_07.png'
plane_1.dataUrl = 'https://doi.org/10.17182/hepdata.91760.v2/r1'
plane_1.setSources(dataLabels= ['efficiencyMap'],
                    dataFiles=['orig/EffFromPatch_SRMET.csv'],
                    units = [ None],
                    coordinates = [{x: 0, y:1, 'value' : -1}],
                    dataFormats=['csv'])
plane_1.addSource('obsExclusion', 'orig/ObservedExclusion.csv', dataFormat = 'csv')
plane_1.addSource('expExclusion', 'orig/ExpectedExclusion.csv', dataFormat = 'csv')

datasets.append(dataset)
#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('SR_MU')
dataset.setInfo( dataType = 'efficiencyMap', dataId = 'SR_MU',
                    observedN = 1, expectedBG = 1.88, bgError=np.sqrt(0.2**2 + 0.28**2), 
                    upperLimit = '0.0272*fb', expectedUpperLimit = '0.0309*fb' )

#+++++++ next txName block ++++++++++++++
TRPV1jmu = dataset.addTxName('TRPV1jmu')
#TRPV1jmu.validationTarball = "TRPV1jmu.tar.gz"
TRPV1jmu.checked =''
TRPV1jmu.constraint = "{(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > jet,mu+), (anyBSM(2) > jet,mu-)}"
TRPV1jmu.condition =None
TRPV1jmu.source = 'ATLAS'
TRPV1jmu.finalState = None
TRPV1jmu.intermediateState = None
TRPV1jmu.dataMap = {0:(1,'mass',GeV), 1:(1,'totalwidth',GeV), 2:(2,'mass',GeV), 3:(2,'totalwidth',GeV)}

#+++++++ next mass plane block ++++++++++++++
plane_1 = TRPV1jmu.addMassPlane({0 : 'x', 1 : 'y', 2 : 'x', 3 : 'y'})
plane_1.figureUrl = None # 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-33/fig_07.png'
plane_1.dataUrl = 'https://doi.org/10.17182/hepdata.91760.v2/r1'
plane_1.setSources(dataLabels= ['efficiencyMap'],
                    dataFiles=['orig/EffFromPatch_SRMU.csv'],
                    units = [None],
                    coordinates = [{x: 0, y:1, 'value' : -1}],
                    dataFormats=['csv'])
plane_1.addSource('obsExclusion', 'orig/ObservedExclusion.csv', dataFormat = 'csv')
plane_1.addSource('expExclusion', 'orig/ExpectedExclusion.csv', dataFormat = 'csv')

datasets.append(dataset)



databaseCreator.create()
