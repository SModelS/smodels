smodelsOutputDefault = {
'OutputStatus' : {
    'sigmacut' : 0.03,
    'minmassgap' : 5.0,
    'maxcond' : 0.2,
    'ncpus' : -6,
    'model' : 'share.models.mssm',
    'checkinput' : True,
    'doinvisible' : True,
    'docompress' : True,
    'testcoverage' : True,
    'computestatistics' : True,
    'combineanas' : 'CMS-SUS-16-050-agg,CMS-SUS-13-012',
    'combinesrs' : True,
    'file status' : 1,
    'decomposition status' : 1,
    'warnings' : 'Input file ok',
    'input file' : './testFiles/slha/gluino_squarks.slha',
    'database version' : 'unittest312',
    'smodels version' : '3.1.2'
},
'ExptRes' : [
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 6.660644,
        'upper limit (fb)' : 17.1845,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T6WW'],
        'Mass (GeV)' : [[992.87, 269.0, 129.0], [991.72, 269.0, 129.0]],
        'AnalysisID' : 'ATLAS-SUSY-2013-02',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 0.387596,
        'r_expected' : None,
        'Width (GeV)' : [
            [14.81, 0.0013847, 'stable'],
            [14.778, 0.0013847, 'stable']
        ]
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 1.956671,
        'upper limit (fb)' : 6.09651,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T2'],
        'Mass (GeV)' : [[991.43, 129.0], [991.31, 129.0]],
        'AnalysisID' : 'ATLAS-SUSY-2013-02',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 0.3209494,
        'r_expected' : None,
        'Width (GeV)' : [[6.6463, 'stable'], [6.8054, 'stable']]
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 9.06887,
        'upper limit (fb)' : 32.3863,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T5WW'],
        'Mass (GeV)' : [[865.0, 269.0, 129.0], [865.0, 269.0, 129.0]],
        'AnalysisID' : 'ATLAS-SUSY-2013-02',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 0.2800218,
        'r_expected' : None,
        'Width (GeV)' : [
            [0.0456539663, 0.00138466665, 'stable'],
            [0.0456539663, 0.00138466665, 'stable']
        ]
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.1561866,
        'upper limit (fb)' : 1.205,
        'expected upper limit (fb)' : 0.7514,
        'TxNames' : ['T2', 'T1', 'T1tttt'],
        'Mass (GeV)' : None,
        'AnalysisID' : 'CMS-SUS-13-012',
        'DataSetID' : '3NJet6_1000HT1250_600MHTinf',
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 19.5,
        'dataType' : 'efficiencyMap',
        'r' : 0.1296154,
        'r_expected' : 0.2078608,
        'Width (GeV)' : None,
        'nll' : 3.906258,
        'nll_min' : 3.573349,
        'nll_SM' : 4.330763
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 7.529518,
        'upper limit (fb)' : 68.0639,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T1'],
        'Mass (GeV)' : [[865.0, 129.0], [865.0, 129.0]],
        'AnalysisID' : 'CMS-PAS-SUS-15-002',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 2.2,
        'dataType' : 'upperLimit',
        'r' : 0.1106242,
        'r_expected' : None,
        'Width (GeV)' : [
            [0.0456539663, 'stable'],
            [0.0456539663, 'stable']
        ]
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 2115.13,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'Mass (GeV)' : [[269.0, 129.0], [268.9, 129.0]],
        'AnalysisID' : 'CMS-SUS-16-039',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.09792791,
        'r_expected' : None,
        'Width (GeV)' : [
            [0.00138466665, 'stable'],
            [0.00112364506, 'stable']
        ]
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.004140203,
        'upper limit (fb)' : 0.04531092,
        'expected upper limit (fb)' : 0.07675983,
        'TxNames' : ['T6bbHH'],
        'Mass (GeV)' : [[959.4, 268.9, 129.0], [959.4, 268.9, 129.0]],
        'AnalysisID' : 'ATLAS-SUSY-2018-31',
        'DataSetID' : '(combined)',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'combined',
        'r' : 0.0913732,
        'r_expected' : 0.05393711,
        'Width (GeV)' : [
            [13.4903869, 0.00112364506, 'stable'],
            [13.4903869, 0.00112364506, 'stable']
        ],
        'nll' : 48.54424,
        'nll_min' : 48.38229,
        'nll_SM' : 48.38229
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 467.97,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [[269.0, 129.0], [268.9, 129.0]],
        'AnalysisID' : 'CMS-SUS-16-039',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.08526429,
        'r_expected' : None,
        'Width (GeV)' : [
            [0.00138466665, 'stable'],
            [0.00112364506, 'stable']
        ]
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.6619329,
        'upper limit (fb)' : 10.8194,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T1'],
        'Mass (GeV)' : [[865.0, 129.0], [865.0, 129.0]],
        'AnalysisID' : 'ATLAS-SUSY-2013-02',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 0.06118018,
        'r_expected' : None,
        'Width (GeV)' : [
            [0.0456539663, 'stable'],
            [0.0456539663, 'stable']
        ]
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 18.47093,
        'upper limit (fb)' : 481.119,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [[269.0, 129.0], [268.9, 129.0]],
        'AnalysisID' : 'ATLAS-SUSY-2013-12',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 0.03839161,
        'r_expected' : None,
        'Width (GeV)' : [
            [0.00138466665, 'stable'],
            [0.00112364506, 'stable']
        ]
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.0007318257,
        'upper limit (fb)' : 0.4,
        'expected upper limit (fb)' : 0.35,
        'TxNames' : ['T1qqtt', 'T1tttt'],
        'Mass (GeV)' : None,
        'AnalysisID' : 'ATLAS-CONF-2013-037',
        'DataSetID' : 'SRtN3',
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.7,
        'dataType' : 'efficiencyMap',
        'r' : 0.001829564,
        'r_expected' : 0.002090931,
        'Width (GeV)' : None,
        'nll' : 3.015969,
        'nll_min' : 2.989467,
        'nll_SM' : 3.019043
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.03496325,
        'upper limit (fb)' : 42.0614,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T1tttt'],
        'Mass (GeV)' : [[865.0, 129.0], [865.0, 129.0]],
        'AnalysisID' : 'CMS-PAS-SUS-15-002',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 2.2,
        'dataType' : 'upperLimit',
        'r' : 0.000831243,
        'r_expected' : None,
        'Width (GeV)' : [
            [0.0456539663, 'stable'],
            [0.0456539663, 'stable']
        ]
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.003947952,
        'upper limit (fb)' : 32.47149,
        'expected upper limit (fb)' : 0.7914109,
        'TxNames' : ['T1tttt'],
        'Mass (GeV)' : [[865.0, 129.0], [865.0, 129.0]],
        'AnalysisID' : 'CMS-SUS-16-050-agg',
        'DataSetID' : '(combined)',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'combined',
        'r' : 0.0001215821,
        'r_expected' : 0.004988498,
        'Width (GeV)' : [
            [0.0456539663, 'stable'],
            [0.0456539663, 'stable']
        ],
        'nll' : 147.9699,
        'nll_min' : 147.9563,
        'nll_SM' : 148.0373
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.003073679,
        'upper limit (fb)' : 43.074,
        'expected upper limit (fb)' : 55.9236,
        'TxNames' : ['T1tttt'],
        'Mass (GeV)' : [[865.0, 129.0], [865.0, 129.0]],
        'AnalysisID' : 'CMS-PAS-SUS-12-026',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 9.2,
        'dataType' : 'upperLimit',
        'r' : 7.135811e-05,
        'r_expected' : 5.496211e-05,
        'Width (GeV)' : [
            [0.0456539663, 'stable'],
            [0.0456539663, 'stable']
        ]
    }
],
'CombinedRes' : [
    {
        'AnalysisID' : 'CMS-SUS-13-012,CMS-SUS-16-050-agg',
        'r' : 0.1183885,
        'r_expected' : 0.2081383,
        'nll' : 151.8762,
        'nll_min' : 151.572,
        'nll_SM' : 152.3681,
        'Txnames' : ['T2', 'T1', 'T1tttt']
    }
],
'Total xsec for missing topologies (fb)' : 3067.167,
'missing topologies' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 148.1634,
        'element' : '[[[jet,W]],[[jet,jet,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 144.5845,
        'element' : '[[[W]],[[W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 144.0016,
        'element' : '[[[jet,jet,W]],[[jet,jet,jet,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 104.7137,
        'element' : '[[[jet,jet,W]],[[t,b,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 87.27849,
        'element' : '[[[jet]],[[jet,jet,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 86.53583,
        'element' : '[[[jet,jet,higgs]],[[jet,jet,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 75.19827,
        'element' : '[[[jet,W]],[[t,b,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 73.08604,
        'element' : '[[[jet,jet,W]],[[jet,t,b,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 73.08604,
        'element' : '[[[t,b,W]],[[jet,jet,jet,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 64.66376,
        'element' : '[[[jet,W]],[[jet,jet,jet,W]]] (MET,MET)'
    }
],
'Total xsec for missing topologies with displaced decays (fb)' : 0.0,
'missing topologies with displaced decays' : [],
'Total xsec for missing topologies with prompt decays (fb)' : 3067.167,
'missing topologies with prompt decays' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 148.1634,
        'element' : '[[[jet,W]],[[jet,jet,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 144.5845,
        'element' : '[[[W]],[[W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 144.0016,
        'element' : '[[[jet,jet,W]],[[jet,jet,jet,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 104.7137,
        'element' : '[[[jet,jet,W]],[[t,b,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 87.27849,
        'element' : '[[[jet]],[[jet,jet,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 86.53583,
        'element' : '[[[jet,jet,higgs]],[[jet,jet,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 75.19827,
        'element' : '[[[jet,W]],[[t,b,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 73.08604,
        'element' : '[[[jet,jet,W]],[[jet,t,b,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 73.08604,
        'element' : '[[[t,b,W]],[[jet,jet,jet,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 64.66376,
        'element' : '[[[jet,W]],[[jet,jet,jet,W]]] (MET,MET)'
    }
],
'Total xsec for topologies outside the grid (fb)' : 28.62023,
'topologies outside the grid' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 23.32642,
        'element' : '[[[jet]],[[jet,jet]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 3.761819,
        'element' : '[[[b,W]],[[b,W]]] (MET,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 1.531994,
        'element' : '[[[jet]],[[t,t]]] (MET,MET)'
    }
]
}
