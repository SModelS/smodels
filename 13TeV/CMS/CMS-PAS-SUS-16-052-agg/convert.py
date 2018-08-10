#!/usr/bin/env python3

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
from smodels_utils.dataPreparation.datasetCreation import DatasetsFromLatex
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z

DataSetInput.ntoys = 50000

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-PAS-SUS-16-052-agg')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/index.html'
info.sqrts = '13.0*TeV'
info.lumi = 35.9
info.prettyName = 'soft lepton, <= 2 jets'
info.private = False
# info.arxiv = 'https://arxiv.org/abs/1704.07781'
info.comment = 'https://cds.cern.ch/record/2273394, http://inspirehep.net/record/1609006'
info.implementedBy = 'WW'
info.contact = 'SModelS'
# info.combinableWith = "CMS-PAS-SUS-16-050" ## FIXME can be combined with a lot, since it vetoes > 2 jets!

max_datasets = None

aggregate = None

## all 44 signal regions. agreement 97%. around 3.3s per point.

## corr > .6: 23 agg regions, agreement 86%. around 2s per point on commodore.
# aggregate = [[1], [2, 3, 4, 18, 19, 21, 22], [5], [6, 7, 8], [9], [10, 11, 12], [13], [14, 15, 16], [17, 20], [23], [24, 25, 26], [27], [28, 29, 30], [31], [32, 33, 34], [35], [36, 37, 38], [39], [40], [41], [42], [43], [44]]

## corr > .7: 28 agg regions. agreement 94%. around X.Xs per point on commodore.
## aggregate = [[1], [2, 3, 4, 19, 21, 22], [5], [6, 7, 8], [9], [10, 11, 12], [13], [14], [15], [16], [17], [18], [20], [23], [24, 25, 26], [27], [28], [29, 30], [31], [32, 33, 34], [35], [36, 37, 38], [39], [40], [41], [42], [43], [44]]

## corr > .7, removing 0 points (18,28): 26 agg regions. agreement 90%. around 2.2s per point on commodore.
# aggregate = [[1], [2, 3, 4, 19, 21, 22], [5], [6, 7, 8], [9], [10, 11, 12], [13], [14], [15], [16], [17], [20], [23], [24, 25, 26], [27], [29, 30], [31], [32, 33, 34], [35], [36, 37, 38], [39], [40], [41], [42], [43], [44]]

## corr > .7, removing low points (13,16,17,18,23,27,28,31): 20 agg regions. agreement 92%. around 2.5s per point on higgs.
# aggregate = [[1], [2, 3, 4, 19, 21, 22], [5], [6, 7, 8], [9], [10, 11, 12], [14], [15], [20], [24, 25, 26], [29, 30], [32, 33, 34], [35], [36, 37, 38], [39], [40], [41], [42], [43], [44]]

## corr > .7, removing low points (13,16,17,18,23,27,28,29,30,31,42,43): 17 agg regions. agreement 91%. around 1.5s per point on commodore.
# aggregate = [[1], [2, 3, 4, 19, 21, 22], [5], [6, 7, 8], [9], [10, 11, 12], [14], [15], [20], [24, 25, 26], [32, 33, 34], [35], [36, 37, 38], [39], [40], [41], [44]]

## corr > .7, removing low points (13,16,17,18,23,24,25,26,27,28,29,30,31,42,43): 16 agg regions. agreement 91%. around 1.4s per point on commodore.
#aggregate = [[1], [2, 3, 4, 19, 21, 22], [5], [6, 7, 8], [9], [10, 11, 12], [14], [15], [20],  [32, 33, 34], [35], [36, 37, 38], [39], [40], [41], [44]]

## corr > .7, moving low points (13,16,17,18,23,24,25,26,27,28,29,30,31,42,43) to separate region: 17 agg regions. agreement 94%. around 1.6s per point on commodore.
aggregate = [[1], [2, 3, 4, 19, 21, 22], [5], [6, 7, 8], [9], [10, 11, 12], [14], [15], [20],  [32, 33, 34], [35], [36, 37, 38], [39], [40], [41], [44], [13,16,17,18,23,24,25,26,27,28,29,30,31,42,43]]

info.createCovarianceMatrix ( "orig/CMS-PAS-SUS-16-052_Figure-aux_001.root", \
                    "c1/CovarianceSRs", addOrder=False, max_datasets = max_datasets,
                    aggregate = aggregate )

## fixme for now i use the root files of CMS16033, as CMS16050 currently
## only publishes the numbers in a table, not in a digitized format.
#creator = DatasetCreator ( "orig/CMS-SUS-16-033_Figure_009.root:DataObs", \
#            "orig/PostFitHistograms.root:PostFitTotal", readDatasetNames=False )
creator = DatasetsFromLatex ( "orig/tables.tex", max_datasets = max_datasets,
            c_obs=6, c_bg=5, ds_name = "#0", aggregate = aggregate )

dsnames = [] ## collect the dataset names
for agg in creator.aggregate:
    tmp  = []
    for a in agg:
        s = creator.origDataSets[a-1].dataId #.lower()
        # s = s.replace("x","X").replace("y","Y")
        tmp.append ( s )
    dsnames.append ( tmp )

#  T2bbWWoff.txt  T6bbWWoff.txt
for ctr,dataset in enumerate(creator):
    #+++++++ next txName block ++++++++++++++
    T2bbWWoff = dataset.addTxName('T2bbWWoff')
    T2bbWWoff.checked =''
    T2bbWWoff.constraint = "[[['b','l','nu']],[['b','jet','jet']]]"
    # T2bbWWoff.constraint = "[[['b','mu','nu']],[['b','jet','jet']]]+[[['b','e','nu']],[['b','jet','jet']]]"
    T2bbWWoff.conditionDescription = None
    T2bbWWoff.condition =None
    T2bbWWoff.source ='SModelS'                                                                #+++++++ next mass plane block ++++++++++++++
    T2bbWWoff_1 = T2bbWWoff.addMassPlane([[x,x-y]]*2)
    #----exclusion source----
    #dataId = dataset.dataId.lower()
    #dataId = dataId.replace("x","X").replace("y","Y")
    dataId = dsnames[ctr]
    T2bbWWoff_1.addSource( 'efficiencyMap', 'orig/AccpEffMap_T2tt.root', 'root',
                        objectName = dataId, scale=3.47 )
    T2bbWWoff_1.addSource( 'obsExclusion', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                        'canvas', objectName = "cCONT_", index=6 )
    T2bbWWoff_1.addSource( 'obsExclusionP1', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                        'canvas', objectName = "cCONT_", index=7 )
    T2bbWWoff_1.addSource( 'obsExclusionM1', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                        'canvas', objectName = "cCONT_", index=8 )
    T2bbWWoff_1.addSource( 'expExclusion', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                        'canvas', objectName = "cCONT_", index=3 )
    T2bbWWoff_1.addSource( 'expExclusionP1', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                        'canvas', objectName = "cCONT_", index=4 )
    T2bbWWoff_1.addSource( 'expExclusionM1', 'orig/CMS-PAS-SUS-16-052_Figure_004.root',\
                        'canvas', objectName = "cCONT_", index=5 )
    T2bbWWoff_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/AccpEffMap_T2tt.root"
    T2bbWWoff_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/CMS-PAS-SUS-16-052_Figure-aux_007.png"

    #+++++++ next txName block ++++++++++++++
    T6bbWWoff = dataset.addTxName('T6bbWWoff')
    T6bbWWoff.checked =''
    T6bbWWoff.constraint = "[[['b'],['l','nu']],[['b'],['jet','jet']]]"
    # T6bbWWoff.constraint = "[[['b'],['e','nu']],[['b'],['jet','jet']]]+[[['b'],['mu','nu']],[['b'],['jet','jet']]]"
    T6bbWWoff.conditionDescription = None
    T6bbWWoff.condition =None
    T6bbWWoff.source ='SModelS'                                                                #+++++++ next mass plane block ++++++++++++++
    T6bbWWoff_1 = T6bbWWoff.addMassPlane([[x,x-.5*y,x-y]]*2)
    #----exclusion source----
    T6bbWWoff_1.addSource( 'efficiencyMap', 'orig/AccpEffMap_T2bW.root', 'root',
                        objectName = dataId, scale=3.47 )
    T6bbWWoff_1.addSource( 'obsExclusion', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                        'canvas', objectName = "cCONT_", index=6 )
    T6bbWWoff_1.addSource( 'obsExclusionP1', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                        'canvas', objectName = "cCONT_", index=7 )
    T6bbWWoff_1.addSource( 'obsExclusionM1', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                        'canvas', objectName = "cCONT_", index=8 )
    T6bbWWoff_1.addSource( 'expExclusion', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                        'canvas', objectName = "cCONT_", index=3 )
    T6bbWWoff_1.addSource( 'expExclusionP1', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                        'canvas', objectName = "cCONT_", index=4 )
    T6bbWWoff_1.addSource( 'expExclusionM1', 'orig/CMS-PAS-SUS-16-052_Figure_005.root',\
                        'canvas', objectName = "cCONT_", index=5 )
    T6bbWWoff_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/AccpEffMap_T2bW.root"
    T6bbWWoff_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/CMS-PAS-SUS-16-052_Figure-aux_008.png"

databaseCreator.create()
