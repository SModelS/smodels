#!/usr/bin/env python

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
info = MetaInfoInput('ATLAS-SUSY-2013-11')
info.comment = 'TChiWW being investigated'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.publication = 'http://link.springer.com/article/10.1007/JHEP05(2014)071'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-11/'
info.arxiv = 'http://arxiv.org/abs/1403.5294'
info.contact = '?'
info.prettyName = '2 leptons (e,mu) + Etmiss'
info.supersedes = 'ATLAS-CONF-2013-049'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiWZ = dataset.addTxName('TChiWZ')
TChiWZ.constraint ="[[['W']],[['Z']]]"
TChiWZ.conditionDescription ="None"
TChiWZ.condition ="None"
TChiWZ.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChiWZ = TChiWZ.addMassPlane(2*[[x, y]])
TChiWZ.figure = 'Fig.(aux) 19b'
TChiWZ.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-11/fig_07a.png'
TChiWZ.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1286761/d64'
TChiWZ.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/expectedexclusion_TChiWZ.txt', 'orig/exclusion_TChiWZ.txt', 'orig/limit_TChiWZ.txt'],
                 dataFormats= ['txt', 'txt', 'txt'])

#+++++++ next txName block ++++++++++++++
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint ="[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription ="[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition ="Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu050 = TChipChimSlepSnu.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChipChimSlepSnu050.figure = 'Fig.(aux) 19a'
TChipChimSlepSnu050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-11/fig_05.png'
TChipChimSlepSnu050.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1286761/d63'
TChipChimSlepSnu050.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/expectedexclusion_TChipChimSlepSnu050.txt', 'orig/exclusion_TChipChimSlepSnu050.txt', 'orig/limit_TChipChimSlepSnu050.txt'],
                 dataFormats= ['txt', 'txt', 'txt'])

#+++++++ next txName block ++++++++++++++
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.conditionDescription ="None"
TChiWW.condition ="None"
TChiWW.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChiWW = TChiWW.addMassPlane(2*[[x, y]])
TChiWW.figure = 'Fig.(aux) 18'
TChiWW.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-11/figaux_18.png'
TChiWW.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1286761/d62'
TChiWW.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiWW.txt', 'orig/limit_TChiWW.txt'],
                 dataFormats= ['txt', 'txt'])

#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TSlepSlep = TSlepSlep.addMassPlane(2*[[x, y]])
TSlepSlep.figure = 'Fig.(aux) 20c'
TSlepSlep.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-11/figaux_20c.png'
TSlepSlep.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1286761/d67'
TSlepSlep.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/expectedexclusion_TSlepSlep.txt', 'orig/exclusion_TSlepSlep.txt', 'orig/limit_TSlepSlep.txt'],
                 dataFormats= ['txt', 'txt', 'txt'],units= [None, None, 'fb'])



databaseCreator.create()
