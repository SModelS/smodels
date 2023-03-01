#!/usr/bin/env python3

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
from smodels_utils.dataPreparation import dataHandlerObjects
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z

dataHandlerObjects.trimmingFactor = 3


#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-20-001')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/index.html'
info.sqrts = 13
info.lumi = 137
info.prettyName = '2 OSSF leptons'
info.private = False
info.arxiv = 'https://arxiv.org/abs/2012.08600'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'JHEP 04 (2021) 123'
info.comment = ''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++txName block +++++++++++++++++
T5ZZ=dataset.addTxName('T5ZZ')
T5ZZ.checked=''
T5ZZ.constraint="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.condition=None
T5ZZ.conditionDescription = None
T5ZZ.source="CMS"
T5ZZ.massConstraint=None

#++++++ mass plane block+++++++++

T5ZZ_1 = T5ZZ.addMassPlane(2*[[x,y,1.0]])
T5ZZ_1.figure='Fig. 10'
T5ZZ_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_010.png'
T5ZZ_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_010.root'
T5ZZ_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_010.root'
T5ZZ_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_010.root'

labels = ['expectedUpperLimits','upperLimits','expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1']
formats= [ 'root' ]*8
objNames = ['ul_histo_exp;1','ul_histo;1','ex_exp_smoothed_graph;1','ex_exp_m1s_smoothed_graph;1','ex_exp_p1s_smoothed_graph;1','ex_obs_smoothed_graph;1','ex_obs_m1s_smoothed_graph;1','ex_obs_p1s_smoothed_graph;1']
units = [ "pb" ] + [ "pb" ]+[ None ] * 6 


T5ZZ_1.setSources(dataLabels=labels,
                    dataFiles=['orig/CMS-SUS-20-001_Figure_010.root']*8,
                    dataFormats=formats,objectNames=objNames, units=units )
                    

T5ZZ_2 = T5ZZ.addMassPlane(2*[[x,y,0.]])
T5ZZ_2.figure='Fig. 10'
T5ZZ_2.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_010.png'
T5ZZ_2.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_010.root'
T5ZZ_2.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_010.root'
T5ZZ_2.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_010.root'


T5ZZ_2.setSources(dataLabels=labels,
                    dataFiles=['orig/CMS-SUS-20-001_Figure_010.root']*8,
                    dataFormats=formats,objectNames=objNames, units=units )
T5ZZ_2.axes = None # ""


#+++++txName block +++++++++++++++++
TChiWZ=dataset.addTxName('TChiWZ')
TChiWZ.checked=''
TChiWZ.constraint="[[['W']],[['Z']]]"
TChiWZ.condition=None
TChiWZ.conditionDescription = None
TChiWZ.source="CMS"
TChiWZ.massConstraint=None

#++++++ mass plane block+++++++++

TChiWZ_1 = TChiWZ.addMassPlane(2*[[x,y]])
TChiWZ_1.figure='Fig. 11'
TChiWZ_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_011.png'
TChiWZ_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_011.root'
TChiWZ_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_011.root'
TChiWZ_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_011.root'



TChiWZ_1.setSources(dataLabels=labels,
                    dataFiles=['orig/CMS-SUS-20-001_Figure_011.root']*8,
                    dataFormats=formats,objectNames=objNames, units=units )
                    
                    
                    
lsp_masses = [0., 1.]
#+++++txName block +++++++++++++++++
TChiZZ=dataset.addTxName('TChiZZ')
TChiZZ.validationTarball = "TChipmZZ.tar.gz"
TChiZZ.checked=''
TChiZZ.constraint="[[['Z']],[['Z']]]"
TChiZZ.condition=None
TChiZZ.conditionDescription = None
TChiZZ.source="CMS"
#TChiZZ.massConstraint=None




#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	plane 						= TChiZZ.addMassPlane(2*[[x,lsp]])	
	plane.figure 				= 'Fig. 12a'
	plane.figureUrl 			= 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_012-a.png'
	plane.dataUrl 				= 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_012-a.root'
	plane.setSources(dataLabels=['expectedUpperLimits','upperLimits','obsExclusion'],
                    dataFiles=['orig/CMS-SUS-20-001_Figure_012-a.root','orig/CMS-SUS-20-001_Figure_012-a.root','orig/TChiZZ_excl.csv'],
                    dataFormats=['root','root','csv'],objectNames=['expected_limit;1','observed_limit;1',None],
                    units=['pb','pb',None])

 





#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
# TSlepSlep.validationTarball = "TSlepSlepLeft.tar.gz"
TSlepSlep.validationTarball = "TSmuSmu.tar.gz" ## given per particle (but left+right-handed)
TSlepSlep.checked =""
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane(2*[[x, y]])
TSlepSlep_1.figure = "Fig. 14"
TSlepSlep_1.figureUrl = "https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_014.png"
TSlepSlep_1.dataUrl = 'https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-001/CMS-SUS-20-001_Figure_014.root'

TSlepSlep_1.setSources(dataLabels= ['obsExclusion', 'expExclusion','expExclusionP1','expExclusionM1','obsExclusionP1','obsExclusionM1', 'upperLimits','expectedUpperLimits'],
                 dataFiles= ['orig/CMS-SUS-20-001_Figure_014.root', 'orig/CMS-SUS-20-001_Figure_014.root','orig/CMS-SUS-20-001_Figure_014.root', 'orig/CMS-SUS-20-001_Figure_014.root','orig/CMS-SUS-20-001_Figure_014.root', 'orig/CMS-SUS-20-001_Figure_014.root','orig/CMS-SUS-20-001_Figure_014.root', 'orig/CMS-SUS-20-001_Figure_014.root'],
                 dataFormats= ['root']*8, 
                 objectNames= ['ex_obs_smoothed_graph;1', 'ex_exp_smoothed_graph;1','ex_exp_p1s_smoothed_graph;1', 'ex_exp_m1s_smoothed_graph;1','ex_obs_p1s_smoothed_graph;1', 'ex_obs_m1s_smoothed_graph;1','ul_histo;1',  'ul_histo_exp;1'],
                 units= [None]*6+['pb']*2) #,
#                 indices= [8,3,4,5,9,10,1,2])



databaseCreator.create()

