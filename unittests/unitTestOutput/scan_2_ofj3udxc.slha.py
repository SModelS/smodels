smodelsOutput = {
'OutputStatus' : {
    'sigmacut' : 0.0,
    'minmassgap' : 5.0,
    'maxcond' : 0.2,
    'ncpus' : 30,
    'promptwidth' : 1e-11,
    'stablewidth' : 1e-25,
    'eraseprompt' : 'eCharge,colordim',
    'checkinput' : True,
    'doinvisible' : True,
    'docompress' : True,
    'computestatistics' : True,
    'testcoverage' : True,
    'combinesrs' : True,
    'combineanas' : 'ATLAS-SUSY-2018-22-multibin, CMS-EXO-20-004',
    'reportallsrs' : False,
    'experimentalfeatures' : False,
    'file status' : 1,
    'decomposition status' : 1,
    'warnings' : 'Input file ok',
    'input file' : './testFiles/rfails/scan_2_ofj3udxc.slha',
    'database version' : '3.0.0',
    'smodels version' : '3.0.0'
},
'ExptRes' : [
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.3683661,
        'upper limit (fb)' : 0.7769,
        'expected upper limit (fb)' : 0.7927,
        'TxNames' : ['TRV1'],
        'Mass (GeV)' : [('zp', 2000.0), ('chi', 65.0)],
        'AnalysisID' : 'ATLAS-SUSY-2018-22',
        'DataSetID' : 'SR2j_2200',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'efficiencyMap',
        'r' : 0.4741486,
        'r_expected' : 0.4646979,
        'Width (GeV)' : [('zp', 112.3956), ('chi', 'stable')],
        'nll' : 9.709546,
        'nll_min' : 9.200812,
        'nll_SM' : 9.200812
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 2.004366,
        'upper limit (fb)' : None,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TRV1'],
        'Mass (GeV)' : [('zp', 2000.0), ('chi', 65.0)],
        'AnalysisID' : 'ATLAS-SUSY-2018-22-multibin',
        'DataSetID' : '(combined)',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'combined',
        'r' : None,
        'r_expected' : None,
        'Width (GeV)' : [('zp', 112.3956), ('chi', 'stable')]
    }
],
'Total xsec for missing topologies (fb)' : 3664.733,
'missing topologies' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 2326.642,
        'SMS' : 'PV > (W,W)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 938.6953,
        'SMS' : 'PV > (Z,Z)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 234.168,
        'SMS' : 'PV > (jet,jet)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 58.54201,
        'SMS' : 'PV > (b,b)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 58.5226,
        'SMS' : 'PV > (t,t)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 47.56281,
        'SMS' : 'PV > (MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.6001723,
        'SMS' : 'PV > (ta,ta)'
    }
],
'Total xsec for missing topologies with displaced decays (fb)' : 0.0,
'missing topologies with displaced decays' : [],
'Total xsec for missing topologies with prompt decays (fb)' : 3664.733,
'missing topologies with prompt decays' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 2326.642,
        'SMS' : 'PV > (W,W)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 938.6953,
        'SMS' : 'PV > (Z,Z)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 234.168,
        'SMS' : 'PV > (jet,jet)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 58.54201,
        'SMS' : 'PV > (b,b)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 58.5226,
        'SMS' : 'PV > (t,t)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 47.56281,
        'SMS' : 'PV > (MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.6001723,
        'SMS' : 'PV > (ta,ta)'
    }
],
'Total xsec for topologies outside the grid (fb)' : 0.0,
'topologies outside the grid' : []
}