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
info = MetaInfoInput('ATLAS-CONF-2013-062')
info.comment = ' T5WWLSP060 and T6WWLSP060 originally have xvalue on y-axes, changed by us to M2. Will be supers. by ATLAS-SUSY-2013-20 (no data yet). Other topologies available but not with 3 mass planes.'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-062/ http://cds.cern.ch/record/1557779'
info.prettyName = '1 lepton + jets + ETmiss'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.checked ="VM"
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription ="None"
T5WW.condition ="None"
T5WW.source = "ATLAS"
T5WW.massConstraint = None
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.3 * [[['jet','jet'],['L','nu']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription = "[[['jet','jet'],['L','nu']],[['jet','jet'],['jet','jet']]] > 2.7* [[['jet','jet'],['ta','nu']],[['jet','jet'],['jet','jet']]], [[['jet','jet'],['L','nu']],[['jet','jet'],['jet','jet']]] > 2.7* [[['jet','jet'],['e','nu']],[['jet','jet'],['jet','jet']]]" 
T5WWoff.condition = "Cgtr([[['jet','jet'],['L','nu']],[['jet','jet'],['jet','jet']]],3.*[[['jet','jet'],['ta','nu']],[['jet','jet'],['jet','jet']]]);Cgtr([[['jet','jet'],['L','nu']],[['jet','jet'],['jet','jet']]],3.*[[['jet','jet'],['e','nu']],[['jet','jet'],['jet','jet']]])"
T5WWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T5WWoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T5WW050 = T5WW.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
T5WW050.figure = 'Fig. 13a'
T5WW050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-062/fig_13a.png'
T5WW050.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-062/fig_13a_PRELIMINARY.data'
T5WW050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T5WW.txt', 'orig/Fig13a.txt'],
                 dataFormats= ['txt', 'txt'])
T5WWoff.addMassPlane(T5WW050)



databaseCreator.create()
