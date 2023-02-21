#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse
from json import loads

argparser = argparse.ArgumentParser(description =
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath',
help = 'path to the package smodels_utils',\
type = str )
argparser.add_argument ('-smodelsPath', '--smodelsPath',
help = 'path to the package smodels_utils',\
type = str )
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


SRlist = ['SRDF_0a_cuts','SRDF_0b_cuts','SRDF_0c_cuts','SRDF_0d_cuts','SRDF_0e_cuts','SRDF_0f_cuts','SRDF_0g_cuts','SRDF_0h_cuts','SRDF_0i_cuts','SRDF_1a_cuts','SRDF_1b_cuts','SRDF_1c_cuts','SRDF_1d_cuts','SRDF_1e_cuts','SRDF_1f_cuts','SRDF_1g_cuts','SRDF_1h_cuts','SRDF_1i_cuts','SRSF_0a_cuts','SRSF_0b_cuts','SRSF_0c_cuts','SRSF_0d_cuts','SRSF_0e_cuts','SRSF_0f_cuts','SRSF_0g_cuts','SRSF_0h_cuts','SRSF_0i_cuts','SRSF_1a_cuts','SRSF_1b_cuts','SRSF_1c_cuts','SRSF_1d_cuts','SRSF_1e_cuts','SRSF_1f_cuts','SRSF_1g_cuts','SRSF_1h_cuts','SRSF_1i_cuts']

#+++++++ global info block ++++++++++++++
info 				= MetaInfoInput('ATLAS-SUSY-2018-32')
info.url 			= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-32/'
info.sqrts 			= 13
info.lumi 			= 139
info.prettyName 	= '2 OS leptons + Etmiss'
info.private 		= False
info.arxiv 			= 'arXiv:1908.08215'
info.contact 		= 'atlas-phys-susy-conveners@cern.ch'
info.publication 	= 'https://epjc.epj.org/articles/epjc/abs/2020/02/10052_2019_Article_7594/10052_2019_Article_7594.html'
info.implementedBy  = "Philipp Neuhuber, Timothée Pascal"
info.jsonFiles      = "{'bkgonly_simplified.json' :"+str(SRlist)+"}"
info.includeCRs     = True

TChiWW = {
'name' 		 : 'TChiWW',
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1750597-v1-Exclusion_contour_(exp)_1.csv',
			   'obsExcl'		: 'orig/HEPData-ins1750597-v1-Exclusion_contour_(obs)_1.csv',
			   'effMap'			: 'orig/EffFromPatch_'},

'constraint' : "[[['W']],[['W']]]",
'massPlane'  : 2*[[x, y]]}


srb = ['SR-DF-0J','SR-DF-1J','SR-SF-0J','SR-SF-1J']
mt2 = ['100,105','105,110','110,120','120,140','140,160','160,180','180,220','220,260','260,inf']

TX 	  = TChiWW
obsN  = [[14.,14.,19.,16.,11.,8.,9.,0.,4.],			[12.,12.,14.,15.,7.,4.,5.,3.,3.],					[14.,15.,24.,37.,20.,12.,12.,5.,8.],				[12.,13.,30.,21.,15.,11.,8.,5.,5.]]
expN  = [[14.,11.3,20.,21.7,11.,6.3,6.5,3.2,2.9], [14.7,9.9,14.3,14.8,6.5,4.4,5.6,2.4,2.7],	[15.8,13.9,26.9,33.3,17.4,10.4,13.4,6.7,6.8],	  [17.,13.7,17.,23.3,17.,9.4,12.3,6.5,8.]]
bgErr = [[4.,2.9,4.,3.,2.,0.9,1.2,0.6,0.5],	 [2.9,2.,2.4,2.2,1.,0.9,1.,0.7,0.8],			[2.3,1.9,3.2,3.4,2.4,1.2,1.6,1.1,1.],	     [4.,2.5,4.,2.8,2.3,1.9,1.6,1.4,2.7]]

for i in range(len(srb)):
	for j in range(len(mt2)):
		sr = srb[:][i] + '_mt2=[' + mt2[:][j] + ')'
		print(sr)

		char = [['e','f'],['g','h'],['a','b'],['c','d']]

		if j < 2: num = '0' + str(j+8)
		else: num = str(j+8)

		fig = 'Fig' + num + char[i][0] + '\n' + 'Fig' + num + char[i][1]
		figUrl1 = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-32/figaux_' + num + char[i][0] + '.png'
		figUrl2 = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-32/figaux_' + num + char[i][1] + '.png'

		figUrl = figUrl1 + '\n' + figUrl2

		#+++++++ dataset block ++++++++++++++
		dataset = DataSetInput(SRlist[9*i+j])
		dataset.setInfo(dataType = 'efficiencyMap', dataId = SRlist[9*i+j], observedN = obsN[i][j], expectedBG = expN[i][j], bgError = bgErr[i][j])
		dataset.comment = sr
		#+++++++ next txName block ++++++++++++++
		newTx 							= dataset.addTxName(TX['name'])
		newTx.checked 					= 'False'
		newTx.constraint 				= TX['constraint']
		newTx.conditionDescription 		= None
		newTx.condition 				= None
		newTx.source 					= 'ATLAS'
		#+++++++ next mass plane block ++++++++++++++
		newPlane 						= newTx.addMassPlane(TX['massPlane'])
		newPlane.figure 				= fig
		newPlane.figureUrl 				= figUrl
		newPlane.dataUrl 				= 'https://www.hepdata.net/record/ins1750597?version=1&table=Acceptance ' + srb[i] + '-[' + mt2[j] + ') for C1C1WW grid'
		newPlane.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
							dataFiles 	= [TX['sources']['expExcl'], TX['sources']['obsExcl'], TX['sources']['effMap'] + sr + '.txt'],
                            dataFormats	= ['csv', 'csv', 'txt'])

databaseCreator.create()
