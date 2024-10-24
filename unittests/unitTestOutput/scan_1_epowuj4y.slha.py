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
    'input file' : './testFiles/rfails/scan_1_epowuj4y.slha',
    'database version' : '3.0.0',
    'smodels version' : '3.0.0'
},
'ExptRes' : [
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.0006565354,
        'upper limit (fb)' : 0.7769,
        'expected upper limit (fb)' : 0.7927,
        'TxNames' : ['TRV1'],
        'Mass (GeV)' : [('zp', 2791.7), ('chi', 1378.4)],
        'AnalysisID' : 'ATLAS-SUSY-2018-22',
        'DataSetID' : 'SR2j_2200',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'efficiencyMap',
        'r' : 0.0008450707,
        'r_expected' : 0.0008282268,
        'Width (GeV)' : [('zp', 30.28039), ('chi', 'stable')],
        'nll' : 9.201049,
        'nll_min' : 9.200812,
        'nll_SM' : 9.200812
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.003112708,
        'upper limit (fb)' : None,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TRV1'],
        'Mass (GeV)' : [('zp', 2791.7), ('chi', 1378.4)],
        'AnalysisID' : 'ATLAS-SUSY-2018-22-multibin',
        'DataSetID' : '(combined)',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'combined',
        'r' : None,
        'r_expected' : None,
        'Width (GeV)' : [('zp', 30.28039), ('chi', 'stable')]
    }
],
'Total xsec for missing topologies (fb)' : 34.56013,
'missing topologies' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 22.64352,
        'SMS' : 'PV > (jet,jet)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 5.684396,
        'SMS' : 'PV > (t,t)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 5.660879,
        'SMS' : 'PV > (b,b)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.2838615,
        'SMS' : 'PV > (W,W)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.1463862,
        'SMS' : 'PV > (higgs,higgs)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.141084,
        'SMS' : 'PV > (Z,Z)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 9.383267e-07,
        'SMS' : 'PV > (ta,ta)'
    }
],
'Total xsec for missing topologies with displaced decays (fb)' : 0.0,
'missing topologies with displaced decays' : [],
'Total xsec for missing topologies with prompt decays (fb)' : 34.56013,
'missing topologies with prompt decays' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 22.64352,
        'SMS' : 'PV > (jet,jet)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 5.684396,
        'SMS' : 'PV > (t,t)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 5.660879,
        'SMS' : 'PV > (b,b)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.2838615,
        'SMS' : 'PV > (W,W)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.1463862,
        'SMS' : 'PV > (higgs,higgs)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.141084,
        'SMS' : 'PV > (Z,Z)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 9.383267e-07,
        'SMS' : 'PV > (ta,ta)'
    }
],
'Total xsec for topologies outside the grid (fb)' : 0.0,
'topologies outside the grid' : []
}