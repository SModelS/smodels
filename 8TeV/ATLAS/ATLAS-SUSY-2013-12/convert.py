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
info = MetaInfoInput('ATLAS-SUSY-2013-12')
info.comment = 'TChiWh renamed to TChiWH'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.publication = 'http://link.springer.com/article/10.1007/JHEP04(2014)169'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-12/'
info.arxiv = 'http://arxiv.org/abs/1402.7029'
info.prettyName = '3 leptons (e,mu,tau) + Etmiss'
info.supersedes = 'ATLAS-CONF-2013-035; CONF-2012-154'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiChipmStauL = dataset.addTxName('TChiChipmStauL')
TChiChipmStauL.constraint ="2.*([[['nu'],['ta']],[['ta+'],['ta-']]] + [[['ta'],['nu']] , [['ta+'],['ta-']]]+[[['nu'],['ta']],[['ta-'],['ta+']]] + [[['ta'],['nu']],[['ta-'],['ta+']]])"
TChiChipmStauL.conditionDescription ="[[['nu'],['ta']],[['ta+'],['ta-']]] ~ [[['ta'],['nu']],[['ta+'],['ta-']]],[[['nu'],['ta']],[['ta+'],['ta-']]] ~ [[['nu'],['ta']],[['ta-'],['ta+']]], [[['nu'],['ta']],[['ta+'],['ta-']]] ~ [[['ta'],['nu']],[['ta-'],['ta+']]]"
TChiChipmStauL.condition ="Csim([[['nu'],['ta']],[['ta+'],['ta-']]],[[['ta'],['nu']],[['ta+'],['ta-']]],[[['nu'],['ta']],[['ta-'],['ta+']]],[[['ta'],['nu']],[['ta-'],['ta+']]])"
TChiChipmStauL.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChiChipmStauL050 = TChiChipmStauL.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChiChipmStauL050.figure = "fig 07c"
TChiChipmStauL050.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-12/fig_07c.png"
TChiChipmStauL050.dataUrl = 'Not defined'
TChiChipmStauL050.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exc_stauL_expected.txt', 'orig/exc_stauL_obs.txt', 'orig/TChiChipmStauL.txt'],
                 dataFormats= ['txt', 'txt', 'txt'])

#+++++++ next txName block ++++++++++++++
TChiWH = dataset.addTxName('TChiWH')
TChiWH.constraint ="[[['W']],[['higgs']]]"
TChiWH.conditionDescription ="None"
TChiWH.condition ="None"
TChiWH.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChiWH = TChiWH.addMassPlane(2*[[x, y]])
TChiWH.figure = "fig 7d"
TChiWH.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-12/fig_07d.png"
TChiWH.dataUrl = 'Not defined'
TChiWH.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exc_tchiwh_expected.txt', 'orig/exc_tchiwh_obs.txt', 'orig/TChiWH.txt'],
                 dataFormats= ['txt', 'txt', 'txt'])

#+++++++ next txName block ++++++++++++++
TChiWZ = dataset.addTxName('TChiWZ')
TChiWZ.constraint ="[[['W']],[['Z']]]"
TChiWZ.conditionDescription ="None"
TChiWZ.condition ="None"
TChiWZ.source = "ATLAS"
TChiWZ.massConstraint = None
TChiWZoff = dataset.addTxName('TChiWZoff')
TChiWZoff.constraint ="71.*([[['mu+','mu-']],[['l','nu']]] + [[['e+','e-']],[['l','nu']]])"
TChiWZoff.conditionDescription ="[[['mu+','mu-']],[['l','nu']]] > [[['e+','e-']],[['l','nu']]]"
TChiWZoff.condition ="Cgtr([[['mu+','mu-']],[['l','nu']]],[[['e+','e-']],[['l','nu']]])"
TChiWZoff.massConstraint = [['dm <= 76.0'], ['dm <= 86.0']]
TChiWZoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChiWZ = TChiWZ.addMassPlane(2*[[x, y]])
TChiWZ.figure = "fig 7b"
TChiWZ.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-12/fig_07b.png"
TChiWZ.dataUrl = 'Not defined'
TChiWZ.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exc_tchiwz_expected.txt', 'orig/exc_tchiwz_obs.txt', 'orig/TChiWZ.txt'],
                 dataFormats= ['txt', 'txt', 'txt'])
TChiWZoff.addMassPlane(TChiWZ)

#+++++++ next txName block ++++++++++++++
TChiChipmSlepL = dataset.addTxName('TChiChipmSlepL')
TChiChipmSlepL.constraint ="2.*([[['L+'],['L-']],[['L'],['nu']]] + [[['L+'],['L-']],[['nu'],['L']]] + [[['L-'],['L+']],[['L'],['nu']]] + [[['L-'],['L+']],[['nu'],['L']]])"
TChiChipmSlepL.conditionDescription ="[[['L+'],['L-']],[['L'],['nu']]] ~ [[['L+'],['L-']],[['nu'],['L']]], [[['L+'],['L-']],[['nu'],['L']]] > 2.7*[[['ta+'],['ta-']],[['nu'],['L']]], [[['L+'],['L-']],[['L'],['nu']]] > 2.7*[[['ta+'],['ta-']],[['L'],['nu']]], [[['L+'],['L-']],[['nu'],['L']]] > 2.7*[[['L+'],['L-']],[['nu'],['ta']]], [[['L+'],['L-']],[['L'],['nu']]] > 2.7*[[['L+'],['L-']],[['ta'],['nu']]],[[['L+'],['L-']],[['nu'],['L']]] > 2.7*[[['e+'],['e-']],[['nu'],['L']]], [[['L+'],['L-']],[['L'],['nu']]] > 2.7*[[['e+'],['e-']],[['L'],['nu']]], [[['L+'],['L-']],[['nu'],['L']]] > 2.7*[[['L+'],['L-']],[['nu'],['e']]], [[['L+'],['L-']],[['L'],['nu']]] > 2.7*[[['L+'],['L-']],[['e'],['nu']]]"
TChiChipmSlepL.condition ="Csim([[['L+'],['L-']],[['L'],['nu']]],[[['L+'],['L-']],[['nu'],['L']]],[[['L-'],['L+']],[['nu'],['L']]],[[['L-'],['L+']],[['L'],['nu']]]); Cgtr([[['L+'],['L-']],[['nu'],['L']]],3.*[[['ta+'],['ta-']],[['nu'],['L']]]); Cgtr([[['L+'],['L-']],[['L'],['nu']]],3.*[[['ta+'],['ta-']],[['L'],['nu']]]);Cgtr([[['L+'],['L-']],[['nu'],['L']]],3.*[[['L+'],['L-']],[['nu'],['ta']]]); Cgtr([[['L+'],['L-']],[['L'],['nu']]],3.*[[['L+'],['L-']],[['ta'],['nu']]]);Cgtr([[['L+'],['L-']],[['nu'],['L']]],3.*[[['e+'],['e-']],[['nu'],['L']]]); Cgtr([[['L+'],['L-']],[['L'],['nu']]],3.*[[['e+'],['e-']],[['L'],['nu']]]); Cgtr([[['L+'],['L-']],[['nu'],['L']]],3.*[[['L+'],['L-']],[['nu'],['e']]]); Cgtr([[['L+'],['L-']],[['L'],['nu']]],3.*[[['L+'],['L-']],[['e'],['nu']]])"
TChiChipmSlepL.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepL050 = TChiChipmSlepL.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChiChipmSlepL050.figure = "fig 7a"
TChiChipmSlepL050.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-12/fig_07a.png"
TChiChipmSlepL050.dataUrl = 'Not defined'
TChiChipmSlepL050.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exc_slepL_expected.txt', 'orig/exc_slepL_obs.txt', 'orig/TChiChipmSlepL.txt'],
                 dataFormats= ['txt', 'txt', 'txt'],units= [None, None, 'fb'])



databaseCreator.create()
