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
info = MetaInfoInput('ATLAS-CONF-2013-035')
info.sqrts = '8.0'
info.private = False
info.lumi = '20.7'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-035/'
info.supersededBy = 'ATLAS-SUSY-2013-12'
info.prettyName = '3 leptons (e,mu) + Etmiss'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiWZ = dataset.addTxName('TChiWZ')
TChiWZ.checked ="VM"
TChiWZ.constraint ="[[['W']],[['Z']]]"
TChiWZ.conditionDescription ="None"
TChiWZ.condition ="None"
TChiWZ.massConstraint = None
TChiWZ.source = 'ATLAS'
TChiWZoff = dataset.addTxName('TChiWZoff')
TChiWZoff.checked ="VM"
TChiWZoff.constraint ="30.3*([[['L','nu']],[['e+','e-']]]+[[['L','nu']],[['mu+','mu-']]]+[[['L','nu']],[['ta+','ta-']]])"
TChiWZoff.conditionDescription ="[[['mu','nu']],[['e+','e-']]]>[[['e','nu']],[['e+','e-']]],[[['e','nu']],[['e+','e-']]]>[[['ta','nu']],[['e+','e-']]],[[['mu','nu']],[['mu+','mu-']]]>[[['e','nu']],[['mu+','mu-']]],[[['e','nu']],[['mu+','mu-']]]>[[['ta','nu']],[['mu+','mu-']]],[[['mu','nu']],[['ta+','ta-']]]>[[['e','nu']],[['ta+','ta-']]],[[['e','nu']],[['ta+','ta-']]]>[[['ta','nu']],[['ta+','ta-']]],[[['L','nu']],[['mu+','mu-']]]>[[['L','nu']],[['e+','e-']]],[[['L','nu']],[['e+','e-']]]>[[['L','nu']],[['ta+','ta-']]]"
TChiWZoff.condition ="Cgtr([[['mu','nu']],[['e+','e-']]],[[['e','nu']],[['e+','e-']]]);Cgtr([[['e','nu']],[['e+','e-']]],[[['ta','nu']],[['e+','e-']]]);Cgtr([[['mu','nu']],[['mu+','mu-']]],[[['e','nu']],[['mu+','mu-']]]);Cgtr([[['e','nu']],[['mu+','mu-']]],[[['ta','nu']],[['mu+','mu-']]]);Cgtr([[['mu','nu']],[['ta+','ta-']]],[[['e','nu']],[['ta+','ta-']]]);Cgtr([[['e','nu']],[['ta+','ta-']]],[[['ta','nu']],[['ta+','ta-']]]);Cgtr([[['L','nu']],[['mu+','mu-']]],[[['L','nu']],[['e+','e-']]]);Cgtr([[['L','nu']],[['e+','e-']]],[[['L','nu']],[['ta+','ta-']]])"
TChiWZoff.massConstraint = [['dm <= 76.0'], ['dm <= 86.0']]
TChiWZoff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
TChiWZ = TChiWZ.addMassPlane(2*[[x, y]])
TChiWZ.figure = 'Fig.(aux) 1b'
TChiWZ.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-035/figaux_01b.png'
TChiWZ.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-035/figaux_01b_PRELIMINARY.data'
TChiWZ.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/TChiWZ_exc.dat', 'orig/TChiWZ.txt'],
                 dataFormats= ['txt', 'txt'])
TChiWZoff.addMassPlane(TChiWZ)

#+++++++ next txName block ++++++++++++++
TChiChipmSlepL = dataset.addTxName('TChiChipmSlepL')
TChiChipmSlepL.checked ="VM"
TChiChipmSlepL.constraint ="2.*([[['L'],['L']],[['L'],['nu']]] + [[['L'],['L']],[['nu'],['L']]])"
TChiChipmSlepL.conditionDescription ="[[['L'],['L']],[['L'],['nu']]] ~ [[['L'],['L']],[['nu'],['L']]], [[['L'],['L']],[['nu'],['L']]] > 2.7*[[['ta'],['ta']],[['nu'],['L']]], [[['L'],['L']],[['L'],['nu']]] > 2.7*[[['ta'],['ta']],[['L'],['nu']]], [[['L'],['L']],[['nu'],['L']]] > 2.7*[[['L'],['L']],[['nu'],['ta']]], [[['L'],['L']],[['L'],['nu']]] > 2.7*[[['L'],['L']],[['ta'],['nu']]],[[['L'],['L']],[['nu'],['L']]] > 2.7*[[['e'],['e']],[['nu'],['L']]], [[['L'],['L']],[['L'],['nu']]] > 2.7*[[['e'],['e']],[['L'],['nu']]], [[['L'],['L']],[['nu'],['L']]] > 2.7*[[['L'],['L']],[['nu'],['e']]], [[['L'],['L']],[['L'],['nu']]] > 2.7*[[['L'],['L']],[['e'],['nu']]]"
TChiChipmSlepL.condition ="Csim([[['L'],['L']],[['L'],['nu']]],[[['L'],['L']],[['nu'],['L']]]); Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['ta'],['ta']],[['nu'],['L']]]); Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['ta'],['ta']],[['L'],['nu']]]);Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['L'],['L']],[['nu'],['ta']]]); Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['L'],['L']],[['ta'],['nu']]]);Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['e'],['e']],[['nu'],['L']]]); Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['e'],['e']],[['L'],['nu']]]); Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['L'],['L']],[['nu'],['e']]]); Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['L'],['L']],[['e'],['nu']]])"
TChiChipmSlepL.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepL050 = TChiChipmSlepL.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChiChipmSlepL050.figure = 'Fig.(aux) 1a'
TChiChipmSlepL050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-035/figaux_01a.png'
TChiChipmSlepL050.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-035/figaux_01a_PRELIMINARY.data'
TChiChipmSlepL050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiChipmSlepL.txt', 'orig/TChiChipmSlepL.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])




databaseCreator.create()
