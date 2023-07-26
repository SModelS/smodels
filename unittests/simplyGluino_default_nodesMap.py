smodelsOutputDefault = {
'OutputStatus' : {
    'sigmacut' : 10.0,
    'minmassgap' : 5.0,
    'maxcond' : 0.2,
    'ncpus' : 1,
    'model' : 'share.models.mssm',
    'promptwidth' : 1e-11,
    'stablewidth' : 1e-25,
    'eraseprompt' : 'spin,eCharge,colordim',
    'checkinput' : True,
    'doinvisible' : True,
    'docompress' : True,
    'computestatistics' : True,
    'testcoverage' : True,
    'combinesrs' : False,
    'combineanas' : 'ATLAS-CONF-2013-037,CMS-SUS-13-012',
    'reportallsrs' : False,
    'experimentalfeatures' : False,
    'file status' : 1,
    'decomposition status' : 1,
    'warnings' : 'Input file ok',
    'input file' : './testFiles/slha/simplyGluino.slha',
    'database version' : 'unittest300-beta',
    'smodels version' : '3.0.0-beta'
},
'SMS Decomposition' : [
    {
        'ID' : 1,
        'SMS' : '(PV > gluino(1),gluino(2)), (gluino(1) > N1,q,q), (gluino(2) > N1,q,q)',
        'Masses (GeV)' : [
            ('gluino', 675.0),
            ('gluino', 675.0),
            ('N1', 200.0),
            ('N1', 200.0)
        ],
        'PIDs' : [
            ('gluino', 1000021),
            ('gluino', 1000021),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 572.0, 'xsec 13.0 TeV': 4310.0}
    }
],
'ExptRes' : [
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 4309.035,
        'upper limit (fb)' : 167.641,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T1'],
        'Mass (GeV)' : [('gluino', 675.0), ('N1', 200.0)],
        'AnalysisID' : 'CMS-PAS-SUS-15-002',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 2.2,
        'dataType' : 'upperLimit',
        'r' : 25.70394,
        'r_expected' : None,
        'Width (GeV)' : [('gluino', float('inf')), ('N1', 0.0)],
        'TxNames weights (fb)' : {'T1': 4309.03465},
        'Nodes Map' : {
            0 : 'PV',
            1 : 'gluino',
            2 : 'gluino',
            5 : 'N1',
            3 : 'q',
            4 : 'q',
            8 : 'N1',
            6 : 'q',
            7 : 'q'
        }
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 572.1689,
        'upper limit (fb)' : 38.1378,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T1'],
        'Mass (GeV)' : [('gluino', 675.0), ('N1', 200.0)],
        'AnalysisID' : 'ATLAS-SUSY-2013-02',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 15.00267,
        'r_expected' : None,
        'Width (GeV)' : [('gluino', float('inf')), ('N1', 0.0)],
        'TxNames weights (fb)' : {'T1': 572.168935},
        'Nodes Map' : {
            0 : 'PV',
            1 : 'gluino',
            2 : 'gluino',
            5 : 'N1',
            3 : 'q',
            4 : 'q',
            8 : 'N1',
            6 : 'q',
            7 : 'q'
        }
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 1.716507,
        'upper limit (fb)' : 0.3857116,
        'expected upper limit (fb)' : 0.2637784,
        'TxNames' : ['T1'],
        'Mass (GeV)' : [('gluino', 675.0), ('N1', 200.0)],
        'AnalysisID' : 'CMS-SUS-13-012',
        'DataSetID' : '6NJet8_1000HT1250_450MHTinf',
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 19.5,
        'dataType' : 'efficiencyMap',
        'r' : 4.450234,
        'r_expected' : 6.507382,
        'Width (GeV)' : [('gluino', float('inf')), ('N1', 0.0)],
        'TxNames weights (fb)' : {'T1': 1.716506805},
        'Nodes Map' : {
            0 : 'PV',
            1 : 'gluino',
            2 : 'gluino',
            5 : 'N1',
            3 : 'q',
            4 : 'q',
            8 : 'N1',
            6 : 'q',
            7 : 'q'
        },
        'likelihood' : 1.267167e-11,
        'l_max' : 0.07085462,
        'l_SM' : 0.04875737
    }
],
'CombinedRes' : [
    {
        'AnalysisID' : 'CMS-SUS-13-012',
        'r' : 4.450234,
        'r_expected' : 6.507382,
        'likelihood' : 1.267167e-11,
        'l_max' : 0.07085462,
        'l_SM' : 0.04875737
    }
],
'Total xsec for missing topologies (fb)' : 0.0,
'missing topologies' : [],
'Total xsec for missing topologies with displaced decays (fb)' : 0.0,
'missing topologies with displaced decays' : [],
'Total xsec for missing topologies with prompt decays (fb)' : 0.0,
'missing topologies with prompt decays' : [],
'Total xsec for topologies outside the grid (fb)' : 0.0,
'topologies outside the grid' : []
}
