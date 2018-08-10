#!/usr/bin/env python


#TO RUN THIS FILE: python templateNew.py -smodelsPath <path-to-the-smodels-folder>  -utilsPath <path-to-the-smodels-utils-folder>
#This file should be run inside the experimental result forlder


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
type = types.StringType)
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
help = 'path to the package smodels_utils',\
type = types.StringType)
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
info = MetaInfoInput('<Experimental Result ID (e.g. ATLAS-CONF-2013-007)')
info.comment = '<Optional comments about the experimental result>'
info.sqrts = '<center-of-mass energy in TeV (e.g. 8)'
info.private = <True/False> #Flag for private (internal, non-public) results
info.lumi = '<luminosity in fb^-1> (e.g. 20.1)'
info.url = '<URL address to the experimental result (e.g. https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-007/)>'
info.prettyName = '<Some pretty name for the result (e.g. ATLAS SS+b)>'
info.implementedBy = '<Implementation author>'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('<name of dataset folder (e.g. data)>')
dataset.setInfo(dataType = '<type of data (upperLimit or efficiencyMap)', dataId = '<datasetId, e.g. None or SR3b>')
#For efficiencyMap results it is also possible to already set the statistics at this level:
#..., observedN = 1, expectedBG = 2.2 , bgError = 0.8, upperLimit = '1.927E-01*fb', expectedUpperLimit = '2.438E-01*fb')

#+++++++ next txName block ++++++++++++++
Tx = dataset.addTxName('<txname label (e.g. T6ttWW)')
Tx.checked ='<Name of whoever checked this txname (or None)>'
Tx.constraint ="<string describing the txname (e.g. [[['t+'],['W-']],[['t-'],['W+']]])>"
Tx.conditionDescription ="<outdated? (set to None)>"
Tx.condition ="<string describing the conditions (e.g. cGtr([[['mu+','mu-']],[['l','nu']]],[[['e+','e-']],[['l','nu']]]))"
Tx.massConstraint = "list with mass conditions (e.g. None or [['dm >= 169.0', 'dm <= 76.0'], ['dm >= 169.0', 'dm <= 76.0']])"
Tx.source = '<source of the data (e.g. CMS, ATLAS, fastlim,...)'
#+++++++ next txName block ++++++++++++++
#Repeat above block structure for additional txnames (if they exist)

#Define mass planes:
#+++++++ next mass plane block ++++++++++++++
Tx_massplane = Tx.addMassPlane("nested list describing the plane axes (e.g. 2*[[x, y, 60.0]] means m_mother =x, m_intermediate = y, m_LSP = 60 GeV in both branches)")
Tx_massplane.figure = '<figure number from official resul> (e.g. Fig.15a)'
Tx_massplane.figureUrl = '<URL address for the figure (e.g. https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-007/fig_15a.png)>'
Tx_massplane.dataUrl = '<URL address for the data (e.g. HEPDATA link or https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-007/)>'
#Associate data files with their respective objects:
#A few different ways of doing it:
Tx_massplane.addSource(dataLabel = "type of data (e.g. 'obsExclusion', 'upperLimits', 'efficiencyMap',...)",
                        dataFile = "address to the data file (e.g. orig/exclusion_T6ttWWLSP060.txt)",
                        dataFormat = "data format  (e.g. 'root', 'txt', 'svg')",
                        objectName = "(Optional) For root files, specify the name of the object contained in the root file here (e.g. 'interpret'", 
                        index = "(Optional) For root files, specify the index for objects in listOfPrimitives of ROOT.TCanvas",
                        unit =  "(Optional) Specify the data unit (e.g. 'fb', None,...). The default it none.",
                        coordinateMap = "(Optional) Specify how the axes variables defined for the plane map to the columns in the data file (e.g. {x : 0, y : 1, 'ul' : 2})",
                        scale = "(Optional) Float number to rescale the data (e.g. 0.001 will divide all the upperlimit/efficiencies by 1000")

#You can also set several sources at once:
Tx_massplane.setSources(dataLabels= ['<list of data labels'],
                 dataFiles= ['<list of data files>'],
                 dataFormats= ['<list of data formats>'],units= ['<list of units>'],...)



#+++++++ next mass plane block ++++++++++++++
#Repeat above block structure for additional mass planes (if they exist)

#If exactly the same plane applies to distinct txnames (for instance the on-shell and off-shell ones), one can do:
#(to avoid repeating the same information twice)
TxB.addMassPlane(Tx_massplane)



#Repeat everything for other datasets/txnames...


#Finally, create the database entry:
databaseCreator.create()
