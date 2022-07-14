#!/usr/bin/env python3

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse
import numpy as np

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
from smodels_utils.dataPreparation.dataHandlerObjects import hbar
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-EXO-19-010')
info.url = "http://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-19-010/"
info.sqrts = 13
info.lumi = 101.0
info.prettyName = "disappearing track"
info.private = False
info.type = 'displaced'
info.arxiv =  'https://arxiv.org/abs/2004.05153'
info.contact = 'CMS collaboration'
info.publication ='Phys. Lett. B 806 (2020) 135502'
info.comment = 'Combined datasets from the 2017, 20018A and 2018B runs. The chargino+chargino efficiencies were multiplied by 2, as instructed by private communications with the analysis authors.'
particlesFile = os.path.abspath('../../../databaseParticles.py')


datasets = {'SR_nlay4' : {'observedN' : 17+5+11, 'expectedBG' : 12.2+7.3+10.3,
                            'bgError' : np.sqrt(1.1**2+4.7**2+1.1**2+3.5**2+1.0**2+5.4**2)},
            'SR_nlay5' : {'observedN' : 4+0+2, 'expectedBG' : 2.1+0.6+1.0,
                            'bgError' : np.sqrt(0.4**2+0.6**2+0.6**2+0.3**2+0.7**2+0.3**2)},
            'SR_nlay6p' : {'observedN' : 6+2+1, 'expectedBG' : 6.7+1.8+5.7,
                            'bgError' : np.sqrt(1.1**2+0.7**2+0.6**2+0.2**2+1.2**2+0.6**2)}
            }

for ds in datasets:
#+++++++ dataset block ++++++++++++++
    dataset = DataSetInput(ds)
    dataset.setInfo(dataType = 'efficiencyMap', dataId = ds,
                    observedN = datasets[ds]['observedN'],
                    expectedBG = datasets[ds]['expectedBG'],
                    bgError = datasets[ds]['bgError'])
#+++++++ next txName block ++++++++++++++
    TDTM1F = dataset.addTxName('TDTM1F')
    TDTM1F.validationTarball = "TDTM1M2F.tar.gz"
    TDTM1F.setParticlesFromFile(particlesFile)
    TDTM1F.checked = ''
    TDTM1F.constraint ="{(PV > C1+(1),C1-(2)), (C1+(1) > *anySM,MET), (C1-(2) > *anySM,MET)}"
    TDTM1F.conditionDescription = None
    TDTM1F.condition = None
    TDTM1F.source = 'CMS'

#+++++++ next mass plane block ++++++++++++++
    labels = [ 'efficiencyMap' ]
    formats = [ 'csv']
    TDTM1F_1 = TDTM1F.addMassPlane(2*[[(x,y), x]])

    TDTM1F_1.setSources(dataLabels= labels,
                     dataFiles= ['./orig/TDTM1F_%s_combined.csv'%ds],
                     dataFormats= formats, coordinates = [{x : 1, y: 0, 'value' : 2}] )
    TDTM1F_1.addSource(dataLabel='obsExclusion',
                        dataFile='orig/Wino_ExclusionObs.csv',
                        coordinateMap = {x : 1, y: 0, 'value' : None}, dataFormat = 'csv')


    #Add second plane for interpolation (assume 1.0 GeV is still safe)
    TDTM1F_2 = TDTM1F.addMassPlane(2*[[(x,y), x-1.5]])
    TDTM1F_2.setSources(dataLabels= labels,
                     dataFiles= ['./orig/TDTM1F_%s_combined.csv'%ds],
                     dataFormats= formats, coordinates = [{x : 1, y: 0, 'value' : 2}] )

    #+++++++ next txName block ++++++++++++++
    TDTM2F = dataset.addTxName('TDTM2F')
    TDTM2F.validationTarball = "TDTM1M2F.tar.gz"
    TDTM2F.setParticlesFromFile(particlesFile)
    TDTM2F.checked = ''
    TDTM2F.constraint ="{(PV > C1(1),MET), (C1(1) > *anySM,MET)}"
    TDTM2F.conditionDescription = None
    TDTM2F.condition = None
    TDTM2F.source = 'CMS'

    #+++++++ next mass plane block ++++++++++++++
    TDTM2F_1 = TDTM2F.addMassPlane([[(x,y), x],[x]])

    TDTM2F_1.setSources(dataLabels= labels,
                     dataFiles= ['./orig/TDTM2F_%s_combined.csv'%ds],
                     dataFormats= formats, coordinates = [{x : 1, y: 0, 'value' : 2}] )
    TDTM2F_1.addSource(dataLabel='obsExclusion',
                        dataFile='orig/Wino_ExclusionObs.csv',
                        coordinateMap = {x : 1, y: 0, 'value' : None}, dataFormat = 'csv')

    #Add second plane for interpolation (assume 1.0 GeV is still safe)
    TDTM2F_2 = TDTM2F.addMassPlane([[(x,y), x-1.5],[x-1.5]])
    TDTM2F_2.setSources(dataLabels= labels,
                     dataFiles= ['./orig/TDTM2F_%s_combined.csv'%ds],
                     dataFormats= formats, coordinates = [{x : 1, y: 0, 'value' : 2}] )




databaseCreator.create()
