#!/usr/bin/env python

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
type = str)
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
help = 'path to the package smodels_utils',\
type = str)
args = argparser.parse_args()

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
info = MetaInfoInput('ATLAS-CONF-2013-036')
info.comment = 'Will be superseded by SUSY-2013-13 (no data available yet)'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.7'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-036/'
info.prettyName = 'ATLAS multileptons'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiChiSlepSlep = dataset.addTxName('TChiChiSlepSlep')
TChiChiSlepSlep.checked ="VM"
TChiChiSlepSlep.constraint ="[[['e+'],['e-']],[['e+'],['e-']]]+[[['e-'],['e+']],[['e-'],['e+']]]+[[['e+'],['e-']],[['e-'],['e+']]] + [[['mu+'],['mu-']],[['mu+'],['mu-']]]+[[['mu-'],['mu+']],[['mu-'],['mu+']]]+[[['mu+'],['mu-']],[['mu-'],['mu+']]] + [[['e+'],['e-']],[['mu+'],['mu-']]]+[[['e-'],['e+']],[['mu-'],['mu+']]]+[[['e+'],['e-']],[['mu-'],['mu+']]]"
TChiChiSlepSlep.conditionDescription ="[[['mu+'],['mu-']],[['e+'],['e-']]] ~ 2.*[[['mu+'],['mu-']],[['mu+'],['mu-']]], [[['mu+'],['mu-']],[['e+'],['e-']]] ~ 2.*[[['e+'],['e-']],[['e+'],['e-']]], [[['mu-'],['mu+']],[['e-'],['e+']]] ~ 2.*[[['e-'],['e+']],[['e-'],['e+']]], [[['mu-'],['mu+']],[['e-'],['e+']]] ~ 2.*[[['mu-'],['mu+']],[['mu-'],['mu+']]], [[['mu+'],['mu-']],[['e-'],['e+']]] ~ 2.*[[['e+'],['e-']],[['e-'],['e+']]], [[['mu+'],['mu-']],[['e-'],['e+']]] ~ 2.*[[['mu+'],['mu-']],[['mu-'],['mu+']]]"
TChiChiSlepSlep.condition ="Csim([[['mu+'],['mu-']],[['e+'],['e-']]],2.*[[['mu+'],['mu-']],[['mu+'],['mu-']]],2.*[[['e+'],['e-']],[['e+'],['e-']]]); Csim([[['mu-'],['mu+']],[['e-'],['e+']]],2.*[[['e-'],['e+']],[['e-'],['e+']]],2.*[[['mu-'],['mu+']],[['mu-'],['mu+']]]); Csim([[['mu+'],['mu-']],[['e-'],['e+']]],[[['e+'],['e-']],[['e-'],['e+']]],[[['mu+'],['mu-']],[['mu-'],['mu+']]])"
TChiChiSlepSlep.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
TChiChiSlepSlepD080 = TChiChiSlepSlep.addMassPlane([[x+80., x+80.-y, x],[x+75., x+80.-y, x]])
TChiChiSlepSlepD080.figure = 'Fig.(aux) 1a'
TChiChiSlepSlepD080.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-036/figaux_01a.png'
TChiChiSlepSlepD080.dataUrl = 'Not defined'
TChiChiSlepSlepD080.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/TChiChiSlepSlepD080_excl.txt', 'orig/TChiChiSlepSlepD080.txt'],
                 dataFormats= ['txt', 'txt'],objectNames= [None, 'None'],indices= [None, 'None'],units= [None, 'fb'])
#+++++++ next mass plane block ++++++++++++++
TChiChiSlepSlep050 = TChiChiSlepSlep.addMassPlane(2*[[y+x, 0.5*(y+x+x), x]])
TChiChiSlepSlep050.figure = 'Fig.(aux) 1b'
TChiChiSlepSlep050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-036/fig_07b.png'
TChiChiSlepSlep050.dataUrl = 'Not defined'
TChiChiSlepSlep050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/TChiChiSlepSlep050_excl.txt', 'orig/TChiChiSlepSlep050.txt'],
                 dataFormats= ['txt', 'txt'],objectNames= [None, 'None'],indices= [None, 'None'],units= [None, 'fb'])



databaseCreator.create()
