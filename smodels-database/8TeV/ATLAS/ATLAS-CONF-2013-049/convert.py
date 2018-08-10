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
info = MetaInfoInput('ATLAS-CONF-2013-049')
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-049/ https://cds.cern.ch/record/1547565'
info.supersededBy = 'ATLAS-SUSY-2013-11'
info.prettyName = '2 leptons (e,mu) + Etmiss'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.checked ="VM"
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TSlepSlep = TSlepSlep.addMassPlane(2*[[x, y]])
TSlepSlep.figure = 'Fig.(aux) 16a'
TSlepSlep.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-049/fig_08c.png'
TSlepSlep.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-049/figaux_16a_PRELIMINARY.data'
TSlepSlep.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/TSlepSlep_exc.dat', 'orig/TSlepSlep.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])



databaseCreator.create()
