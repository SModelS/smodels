#!/usr/bin/env python3

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""

import sys
from smodels_utils.dataPreparation.argParser import getParserArgs
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z

getParserArgs()

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2018-41')
# info.comment = ''
info.sqrts = '13.0'
info.private = False
info.lumi = '139.'
info.publication = "https://doi.org/10.1103/PhysRevD.104.112010"
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-41/'
info.arxiv = 'https://arxiv.org/abs/2108.07586'
info.prettyName = 'Boosted hadronic EWK searches'
info.implementedBy = 'WW'

# DatayieldsandbackgroundbreakdowninSR.csv
bgrounds = """
SR-4Q-WW,1.9,0.4,-0.4
SR-4Q-WZ,3.4,0.7,-0.7
SR-4Q-ZZ,1.9,0.5,-0.5
SR-4Q-VV,3.9,0.8,-0.8
SR-2B2Q-WZ,1.6,0.4,-0.4
SR-2B2Q-Wh,1.9,0.7,-0.7
SR-2B2Q-ZZ,1.7,0.5,-0.5
SR-2B2Q-Zh,1.6,0.5,-0.5
SR-2B2Q-VZ,2.2,0.6,-0.6
SR-2B2Q-Vh,2.5,0.8,-0.8
"""
data="""
SR-4Q-WW,2.0
SR-4Q-WZ,3.0
SR-4Q-ZZ,1.0
SR-4Q-VV,3.0
SR-2B2Q-WZ,2.0
SR-2B2Q-Wh,0.0
SR-2B2Q-ZZ,2.0
SR-2B2Q-Zh,1.0
SR-2B2Q-VZ,2.0
SR-2B2Q-Vh,1.0
"""

yields = {}
for line in data.split("\n"):
    line = line.strip()
    if len(line)==0:
        continue
    tokens = line.split(",")
    yields[ tokens[0] ] = { "obsN": int( float ( tokens[1] ) ) }

for line in bgrounds.split("\n"):
    line = line.strip()
    if len(line)==0:
        continue
    tokens = line.split(",")
    yields[ tokens[0] ] [ "bgExp" ] = float( tokens[1] )
    err = max ( float ( tokens[2] ), abs ( float ( tokens[3] ) ) )
    yields[ tokens[0] ] [ "expErr" ] = err

# tab_06
from smodels.tools.physicsUnits import fb
yields["SR-4Q-WW"]["ul"]=0.032*fb
yields["SR-4Q-WZ"]["ul"]=0.036*fb
yields["SR-4Q-ZZ"]["ul"]=0.025*fb
yields["SR-4Q-VV"]["ul"]=0.034*fb
yields["SR-2B2Q-WZ"]["ul"]=0.033*fb
yields["SR-2B2Q-Wh"]["ul"]=0.022*fb
yields["SR-2B2Q-ZZ"]["ul"]=0.033*fb
yields["SR-2B2Q-Zh"]["ul"]=0.026*fb
yields["SR-2B2Q-VZ"]["ul"]=0.032*fb
yields["SR-2B2Q-Vh"]["ul"]=0.026*fb
yields["SR-4Q-WW"]["eul"]=0.032*fb*4.2/4.5
yields["SR-4Q-WZ"]["eul"]=0.036*fb*5.1/5.0
yields["SR-4Q-ZZ"]["eul"]=0.025*fb*4.1/3.6
yields["SR-4Q-VV"]["eul"]=0.034*fb*5.3/4.7
yields["SR-2B2Q-WZ"]["eul"]=0.033*fb*4.0/4.7
yields["SR-2B2Q-Wh"]["eul"]=0.022*fb*3.9/3.1
yields["SR-2B2Q-ZZ"]["eul"]=0.033*fb*4.1/4.5
yields["SR-2B2Q-Zh"]["eul"]=0.026*fb*3.9/3.6
yields["SR-2B2Q-VZ"]["eul"]=0.032*fb*4.4/4.4
yields["SR-2B2Q-Vh"]["eul"]=0.026*fb*4.4/3.6

for k,v in yields.items():
    yields[k]["eul"] = round ( float(yields[k]["eul"]/fb),4 )*fb

acceptancetables = { "SR-4Q-WW": "t124", "SR-2B2Q-Wh": "t127",
                     "SR-4Q-WZ": "t125", "SR-4Q-ZZ": "t128", 
                     "SR-4Q-VV": "t128", "SR-2B2Q-WZ": "t126",
                     "SR-2B2Q-ZZ": "t126", "SR-2B2Q-Zh": "t130", 
                     "SR-2B2Q-VZ": "t126", "SR-2B2Q-Vh": "t130" }
acceptances = { "SR-4Q-WW": "AcceptanceofC1C1-WWsignalsbySR-4Q-VV.csv",
    "SR-4Q-WZ": "AcceptanceofC1N2-WZsignalsbySR-4Q-VV.csv",
    "SR-4Q-ZZ": "AcceptanceofN2N3-ZZsignalsbySR-4Q-VV.csv",
    "SR-4Q-VV": "AcceptanceofN2N3-ZZsignalsbySR-4Q-VV.csv",
    "SR-2B2Q-Wh": "AcceptanceofC1N2-WhsignalsbySR-2B2Q-Vh.csv",
    "SR-2B2Q-WZ": "AcceptanceofC1N2-WZsignalsbySR-2B2Q-VZ.csv",
    "SR-2B2Q-ZZ": "AcceptanceofC1N2-WZsignalsbySR-2B2Q-VZ.csv",
    "SR-2B2Q-Zh": "AcceptanceofN2N3-ZhsignalsbySR-2B2Q-Vh.csv",
    "SR-2B2Q-VZ": "AcceptanceofC1N2-WZsignalsbySR-2B2Q-VZ.csv",
    "SR-2B2Q-Vh": "AcceptanceofN2N3-hhsignalsbySR-2B2Q-Vh.csv",
}

for SR,data in yields.items():
    #+++++++ dataset block ++++++++++++++
    dataset = DataSetInput( SR )
    dataset.setInfo(dataType = 'efficiencyMap', dataId = SR, 
       observedN = data["obsN"], expectedBG=data["bgExp"], 
       bgError = data["expErr"],
       upperLimit = data["ul"], expectedUpperLimit = data["eul"]
)
        
    acc = "orig/"+ acceptances[SR]
    eff = acc.replace("Acceptance","Efficiency")
    acct = acceptancetables[SR]
    units = [ None ] * 4 + [ "/10000" ]

    if "WZ" in SR:
        TChiWZ = dataset.addTxName('TChiWZ')
        # TChiWZ.checked ="VM"
        TChiWZ.constraint ="[[['W']],[['Z']]]"
        TChiWZ.conditionDescription ="None"
        TChiWZ.condition ="None"
        TChiWZ.source = 'ATLAS'
        #+++++++ next mass plane block ++++++++++++++
        TChiWZ = TChiWZ.addMassPlane(2*[[x, y ]])
        TChiWZ.figure = 'Figaux. 07b'
        TChiWZ.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-41/figaux_07b.png'
        TChiWZ.dataUrl = f'https://doi.org/10.17182/hepdata.104458.v1/{acct}'
        # TChiWZ.exclusionDataUrl = 'https://doi.org/10.17182/hepdata.104458.v1/t34'
        TChiWZ.setSources(dataLabels= ['obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 
                'expExclusion', 'efficiencyMap' ],
                dataFiles= [ 'orig/Obslimiton(W~,B~)simplifiedmodel(C1N2-WZ).csv', 'orig/Obslimit(+1sig)on(W~,B~)simplifiedmodel(C1N2-WZ).csv', 'orig/Obslimit(-1sig)on(W~,B~)simplifiedmodel(C1N2-WZ).csv', 'orig/Explimiton(W~,B~)simplifiedmodel(C1N2-WZ).csv', (acc,eff) ],
                 dataFormats= ['csv'] * 5, units= units, indices = [ None ] * 5 )

    if "Wh" in SR:
        TChiWH = dataset.addTxName('TChiWH')
        # TChiWH.checked ="VM"
        TChiWH.constraint ="[[['W']],[['higgs']]]"
        TChiWH.conditionDescription ="None"
        TChiWH.condition ="None"
        TChiWH.source = 'ATLAS'
        #+++++++ next mass plane block ++++++++++++++
        TChiWH = TChiWH.addMassPlane(2*[[x, y ]])
        TChiWH.figure = 'Figaux. 07c'
        TChiWH.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-41/figaux_07c.png'
        TChiWH.dataUrl = 'https://doi.org/10.17182/hepdata.104458.v1/{acct}'
        # TChiWH.exclusionDataUrl = 'https://doi.org/10.17182/hepdata.104458.v1/t40'
        TChiWH.setSources(dataLabels= ['obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 
                'expExclusion', 'efficiencyMap' ],
                dataFiles= [ 'orig/Obslimiton(W~,B~)simplifiedmodel(C1N2-Wh).csv', 'orig/Obslimit(+1sig)on(W~,B~)simplifiedmodel(C1N2-Wh).csv', 'orig/Obslimit(-1sig)on(W~,B~)simplifiedmodel(C1N2-Wh).csv', 'orig/Explimiton(W~,B~)simplifiedmodel(C1N2-Wh).csv', ( eff, acc ) ],
                 dataFormats= ['csv'] * 5, units=units, indices = [ None ] * 5 )

    if "WW" in SR:
        TChiWW = dataset.addTxName('TChiWW')
        # TChiWW.checked ="VM"
        TChiWW.constraint ="[[['W']],[['W']]]"
        TChiWW.conditionDescription ="None"
        TChiWW.condition ="None"
        TChiWW.source = 'ATLAS'
        #+++++++ next mass plane block ++++++++++++++
        TChiWW = TChiWW.addMassPlane(2*[[x, y ]])
        TChiWW.figure = 'Figaux. 07a'
        TChiWW.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-41/figaux_07a.png'
        TChiWW.dataUrl = 'https://doi.org/10.17182/hepdata.104458.v1/{acct}'
        # TChiWW.exclusionDataUrl = 'https://doi.org/10.17182/hepdata.104458.v1/t40'
        TChiWW.setSources(dataLabels= ['obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 
                'expExclusion', 'efficiencyMap' ],
                dataFiles= [ 'orig/Obslimiton(W~,B~)simplifiedmodel(C1C1-WW).csv', 'orig/Obslimit(+1sig)on(W~,B~)simplifiedmodel(C1C1-WW).csv', 'orig/Obslimit(-1sig)on(W~,B~)simplifiedmodel(C1C1-WW).csv', 'orig/Explimiton(W~,B~)simplifiedmodel(C1C1-WW).csv', (eff,acc) ],
                 dataFormats= ['csv'] * 5, units= units, indices = [ None ] * 5 )

databaseCreator.create()
