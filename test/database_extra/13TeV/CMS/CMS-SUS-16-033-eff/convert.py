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
argparser.add_argument ('-no', '--noUpdate',
help = 'do not update the lastUpdate field.',\
action= "store_true" )
argparser.add_argument ('-r', '--resetValidation',
help = 'reset the validation flag',\
action= "store_true" )
args = argparser.parse_args()

if args.noUpdate:
    os.environ["SMODELS_NOUPDATE"]="1"

if args.resetValidation:                                                                  os.environ["SMODELS_RESETVALIDATION"]="1"

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
from smodels_utils.dataPreparation import dataHandlerObjects
dataHandlerObjects.allowTrimming = True
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput,getSignalRegionsEMBaked,getStatsEMBaked
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z


#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-16-033')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-033/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = '0L + jets + Etmiss (using MHT)'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1704.07781'
info.contact = 'smodels-users@lists.oeaw.ac.at'
info.publication = 'Phys. Rev D 96 (2017) 032003, http://dx.doi.org/10.1103/PhysRevD.96.032003'
info.comment = 'Recast with MadAnalysis5, PAD code http://doi.org/10.7484/INSPIREHEP.DATA.77YH.NBR3 (aggregate regions)'

topos = [ "T1", "T2", "T1tttt", "T1ttttoff", "T1bbbb", "T2tt", "T2ttoff", "T2bb", "TGQ" ]# 
topos.append ( "T3GQ" )
topos.append ( "T5GQ" )

constraints = { "T2": "[[[jet]],[[jet]]]", "T1": "[[[jet,jet]],[[jet,jet]]]", \
                "T2bb": "[[[b]],[[b]]]", "T1bbbb": "[[[b,b]],[[b,b]]]", \
                "T2tt": "[[[t]],[[t]]]", "T1tttt": "[[[t,t]],[[t,t]]]", \
                "T2ttoff": "[[[W,b]],[[W,b]]]",
                "T1ttttoff": "[[[b,b,W,W]],[[b,b,W,W]]]",
                "T5GQ": "[[[jet],[jet,jet]],[[jet,jet]]]",
                "TGQ": "[[[jet]],[[jet,jet]]]", "T3GQ": "[[[jet]],[[jet],[jet]]]" }

figures = { "T1": '12-c', "T1bbbb": '12-b', "T1tttt": '12-a', "T3GQ": "None", \
            "T1ttttoff": '12-a', "T5GQ": "None", "T2ttoff": "13-a",
            "T2bb": '13-b', "T2tt": '13-a', "T2": '13-c', "TGQ": "None" }

baseUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-033/'

nofigures = [ "TGQ", "T3GQ", "T5GQ" ]

def getFigure ( topo ):
    if topo in nofigures:
        return ""
    return "Fig. %s" % figures[topo]

def getFigureUrl ( topo ):
    if topo in nofigures:
        return ""
    return baseUrl+'/CMS-SUS-16-033_Figure_0%s.png' % figures[topo]

def getDataUrl ( topo ):
    if topo in nofigures:
        return ""
    return baseUrl+'/CMS-SUS-16-033_Figure_0%s.root' % figures[topo]

def getLocalFile ( topo ):
    if topo in nofigures:
        return ""
    return 'orig/CMS-SUS-16-033_Figure_0%s.root' % figures[topo]

def getAxis ( topo ):
    if topo in [ "TGQ" ]:
        # return [[x,0.],[y,0.]]
        return [[x,z],[y,z]]
    if topo in [ "T3GQ" ]:
        # return [[x,0.],[y,0.]]
        return [[y,z],[x,y,z]]
    if topo in [ "T5GQ" ]:
        # return [[x,0.],[y,0.]]
        return [[x,y,z],[y,z]]
    return 2*[[x,y]]

stats = getStatsEMBaked()
#print ( "[convert] stats", stats )

datasets={}

dsnames=getSignalRegionsEMBaked("orig/%s.embaked" % "T1" )
#print ( "[convert] dsnames", dsnames )

for dsname in dsnames:
    #+++++++ dataset block ++++++++++++++
    if not dsname in stats:
        print ( "cannot find stats for %s" % dsname )
        continue
    dataset = DataSetInput( dsname )
    dst = stats[dsname]
    dataset.setInfo(dataType = 'efficiencyMap', dataId = dsname,
            observedN = dst["nobs"], expectedBG = dst["nb"] , bgError = dst["deltanb"] )
    datasets[dsname]=dataset
    #+++++txName block +++++++++++++++++

for topo in topos:
    for dsname in dsnames:
        Tx=datasets[dsname].addTxName( topo )
        Tx.checked=''
        Tx.constraint=constraints[topo]
        Tx.condition=None
        Tx.conditionDescription = None
        if topo == "T1ttttoff":
                Tx.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
        Tx.source="SModelS"

        #++++++next mass plane block+++++++++

        Tx_1 = Tx.addMassPlane(getAxis(topo))
        Tx_1.figure=getFigure(topo)
        Tx_1.figureUrl=getFigureUrl(topo)
        if topo in [ "TGQ" ]:
#                Tx_1.axes = "[[x, 0.0], [y, 0.0]]"
                Tx_1.axes = "[[x, 0.0], [y, 0.0]]; [[x, 695.0], [y, 695.0]]; [[x, 995.0], [y, 995.0]]"
        if topo in [ "T3GQ" ]:
#                Tx_1.axes = "[[y, 0.0], [x, y, 0.0]]"
                Tx_1.axes = "[[y, 0.0], [x, y, 0.0]]; [[y, 695.0], [x, y, 695.0]]; [[y, 995.0], [x, y, 995.0]]"
        if topo in [ "T5GQ" ]:
                Tx_1.axes = "[[x, y, 0.0], [y, 0.0]]; [[x, y, 695.0], [y, 695.0]]; [[x, y, 995.0], [y, 995.0]]"
        Tx_1.dataUrl=None # getDataUrl(topo)
        Tx_1.histoDataUrl=getDataUrl(topo)
        Tx_1.exclusionDataUrl=getDataUrl(topo)
        names = ['ExpLim;1','ExpLimSdn;1','ExpLimSup;1','ObsLim;1','ObsLimSdn;1','ObsLimSup;1',dsname ]
        if topo in [ "T2", "T2tt", "T2ttoff" ]:
                names = ['ExpLim2;1','ExpLimSdn2;1','ExpLimSup2;1','ObsLim2;1','ObsLimSdn2;1','ObsLimSup2;1',dsname ]

        #n=len(names)
            #n=7
        #formats=['root']*n
        # names = [None]*6 + [ dsname ]
        formats=['root']*6+['embaked']
        labels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','efficiencyMap']
        units = [None]*7 # (n-1)+['pb']
        exclroot = getLocalFile ( topo )
        if exclroot != "":
            Tx_1.setSources(dataLabels=labels,
                dataFiles=[exclroot]*6+['orig/%s.embaked' % topo ],
                dataFormats=formats,objectNames=names,units=units )
        else:
            Tx_1.setSources(dataLabels=labels[-1:],
                dataFiles=['orig/%s.embaked' % topo ],
                dataFormats=formats[-1:],objectNames=names[-1:],units=units[-1:] )

databaseCreator.create()
