smodelsOutputDefault = {
'OutputStatus' : {
    'sigmacut' : 0.03,
    'minmassgap' : 5.0,
    'minmassgapisr' : 5.0,
    'maxcond' : 0.2,
    'ncpus' : -6,
    'model' : 'share.models.mssm',
    'checkinput' : True,
    'doinvisible' : True,
    'docompress' : True,
    'computestatistics' : False,
    'testcoverage' : True,
    'combineanas' : 'CMS-SUS-13-012,ATLAS-SUSY-2018-31,ATLAS-CONF-2013-037',
    'file status' : 1,
    'decomposition status' : 1,
    'warnings' : 'Input file ok',
    'input file' : './testFiles/slha/simplyGluino.slha',
    'database version' : 'unittest211',
    'smodels version' : '2.2.0'
},
'ExptRes' : [
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 4309.035,
        'upper limit (fb)' : 167.6413,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T1'],
        'Mass (GeV)' : [[675.0, 200.0], [675.0, 200.0]],
        'AnalysisID' : 'CMS-PAS-SUS-15-002',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 2.2,
        'dataType' : 'upperLimit',
        'r' : 25.70389,
        'r_expected' : None,
        'Width (GeV)' : [[1.0, 'stable'], [1.0, 'stable']]
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 572.1689,
        'upper limit (fb)' : 38.13784,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T1'],
        'Mass (GeV)' : [[675.0, 200.0], [675.0, 200.0]],
        'AnalysisID' : 'ATLAS-SUSY-2013-02',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 15.00266,
        'r_expected' : None,
        'Width (GeV)' : [[1.0, 'stable'], [1.0, 'stable']]
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 1.716507,
        'upper limit (fb)' : 0.3857116,
        'expected upper limit (fb)' : 0.2637784,
        'TxNames' : ['T1'],
        'Mass (GeV)' : [[675.0, 200.0], [675.0, 200.0]],
        'AnalysisID' : 'CMS-SUS-13-012',
        'DataSetID' : '6NJet8_1000HT1250_450MHTinf',
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 19.5,
        'dataType' : 'efficiencyMap',
        'r' : 4.450234,
        'r_expected' : 6.507382,
        'Width (GeV)' : [[1.0, 'stable'], [1.0, 'stable']],
        'nll' : 25.091652322862235,
        'nll_min' : 2.6471251067851,
        'nll_SM' : 3.0208989134775686
    }
],
'Total xsec for missing topologies (fb)' : 0.0,
'missing topologies' : [],
'Total xsec for missing topologies with displaced decays (fb)' : 0.0,
'missing topologies with displaced decays' : [],
'Total xsec for missing topologies with prompt decays (fb)' : 0.0,
'missing topologies with prompt decays' : [],
'Total xsec for topologies outside the grid (fb)' : 0.0,
'topologies outside the grid' : [],
'CombinedRes' : [
    {
        'AnalysisID' : 'CMS-SUS-13-012',
        'r' : 4.450234,
        'r_expected' : 6.507382,
        'nll' : 25.091652322862235,
        'nll_min' : 2.6471251067851,
        'nll_SM' : 3.0208989134775686
    }
]
}

