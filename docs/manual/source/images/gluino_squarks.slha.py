smodelsOutput = {
'OutputStatus' : {
    'sigmacut' : 10.0,
    'minmassgap' : 5.0,
    'minmassgapisr' : 1.0,
    'maxcond' : 0.2,
    'ncpus' : 1,
    'model' : 'share.models.mssm',
    'promptwidth' : 1e-11,
    'stablewidth' : 1e-25,
    'ignorepromptqnumbers' : 'spin,eCharge,colordim',
    'checkinput' : True,
    'doinvisible' : True,
    'docompress' : True,
    'testcoverage' : True,
    'computestatistics' : True,
    'combinesrs' : False,
    'combineanas' : 'ATLAS-SUSY-2018-10,CMS-SUS-13-012',
    'pyhfbackend' : 'pytorch',
    'reportallsrs' : False,
    'file status' : 1,
    'decomposition status' : 1,
    'warnings' : 'Input file ok',
    'input file' : 'inputFiles/slha/gluino_squarks.slha',
    'database version' : '3.1.0',
    'smodels version' : '3.1.1'
},
'SMS Decomposition' : [
    {
        'ID' : 1,
        'SMS' : '(PV > N2(1),C1-(2)), (N2(1) > N1,higgs), (C1-(2) > N1~,W-)',
        'Masses (GeV)' : [
            ('N2', 268.9),
            ('C1-', 269.0),
            ('N1', 129.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('N2', 1000023),
            ('C1-', -1000024),
            ('N1', 1000022),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 320.0, 'xsec 13.0 TeV': 712.0}
    },
    {
        'ID' : 2,
        'SMS' : '(PV > N2(1),C1+(2)), (N2(1) > N1,higgs), (C1+(2) > N1,W+)',
        'Masses (GeV)' : [
            ('N2', 268.9),
            ('C1+', 269.0),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('N2', 1000023),
            ('C1+', 1000024),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 68.4, 'xsec 13.0 TeV': 136.0}
    },
    {
        'ID' : 3,
        'SMS' : '(PV > N2(1),C1-(2)), (N2(1) > N1,Z), (C1-(2) > N1~,W-)',
        'Masses (GeV)' : [
            ('N2', 268.9),
            ('C1-', 269.0),
            ('N1', 129.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('N2', 1000023),
            ('C1-', -1000024),
            ('N1', 1000022),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 5.29, 'xsec 13.0 TeV': 13.7}
    },
    {
        'ID' : 4,
        'SMS' : '(PV > N2(1),C1+(2)), (N2(1) > N1,Z), (C1+(2) > N1,W+)',
        'Masses (GeV)' : [
            ('N2', 268.9),
            ('C1+', 269.0),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('N2', 1000023),
            ('C1+', 1000024),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 13.2, 'xsec 13.0 TeV': 26.2}
    },
    {
        'ID' : 5,
        'SMS' : '(PV > C1-(1),C1+(2)), (C1-(1) > N1~,W-), (C1+(2) > N1,W+)',
        'Masses (GeV)' : [
            ('C1-', 269.0),
            ('C1+', 269.0),
            ('N1~', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('C1-', -1000024),
            ('C1+', 1000024),
            ('N1~', -1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 59.6, 'xsec 13.0 TeV': 145.0}
    },
    {
        'ID' : 6,
        'SMS' : '(PV > C1+(1),sd_L(2)), (C1+(1) > N1,W+), (sd_L(2) > q,C1-(6)), (C1-(6) > N1~,W-)',
        'Masses (GeV)' : [
            ('C1+', 269.0),
            ('sd_L', 994.3),
            ('N1', 129.0),
            ('C1-', 269.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('C1+', 1000024),
            ('sd_L', 1000001),
            ('N1', 1000022),
            ('C1-', -1000024),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 4.69, 'xsec 13.0 TeV': 35.1}
    },
    {
        'ID' : 7,
        'SMS' : '(PV > su_R(1),su_L(2)), (su_R(1) > N1,q), (su_L(2) > q,C1+(6)), (C1+(6) > N1,W+)',
        'Masses (GeV)' : [
            ('su_R', 991.3),
            ('su_L', 991.5),
            ('N1', 129.0),
            ('C1+', 269.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_R', 2000002),
            ('su_L', 1000002),
            ('N1', 1000022),
            ('C1+', 1000024),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.51, 'xsec 13.0 TeV': 11.8}
    },
    {
        'ID' : 8,
        'SMS' : '(PV > su_R(1),gluino(2)), (su_R(1) > N1,q), (gluino(2) > q,q,C1-(7)), (C1-(7) > N1~,W-)',
        'Masses (GeV)' : [
            ('su_R', 991.3),
            ('gluino', 865.0),
            ('N1', 129.0),
            ('C1-', 269.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('su_R', 2000002),
            ('gluino', 1000021),
            ('N1', 1000022),
            ('C1-', -1000024),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 21.6, 'xsec 13.0 TeV': 194.0}
    },
    {
        'ID' : 9,
        'SMS' : '(PV > su_R(1),gluino(2)), (su_R(1) > N1,q), (gluino(2) > q,q,C1+(7)), (C1+(7) > N1,W+)',
        'Masses (GeV)' : [
            ('su_R', 991.3),
            ('gluino', 865.0),
            ('N1', 129.0),
            ('C1+', 269.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_R', 2000002),
            ('gluino', 1000021),
            ('N1', 1000022),
            ('C1+', 1000024),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.96, 'xsec 13.0 TeV': 17.5}
    },
    {
        'ID' : 10,
        'SMS' : '(PV > su_R(1),gluino(2)), (su_R(1) > N1,q), (gluino(2) > q,c,C1-(7)), (C1-(7) > N1~,W-)',
        'Masses (GeV)' : [
            ('su_R', 991.3),
            ('gluino', 865.0),
            ('N1', 129.0),
            ('C1-', 269.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('su_R', 2000002),
            ('gluino', 1000021),
            ('N1', 1000022),
            ('C1-', -1000024),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.96, 'xsec 13.0 TeV': 17.5}
    },
    {
        'ID' : 11,
        'SMS' : '(PV > su_R(1),gluino(2)), (su_R(1) > N1,q), (gluino(2) > q,c,C1+(7)), (C1+(7) > N1,W+)',
        'Masses (GeV)' : [
            ('su_R', 991.3),
            ('gluino', 865.0),
            ('N1', 129.0),
            ('C1+', 269.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_R', 2000002),
            ('gluino', 1000021),
            ('N1', 1000022),
            ('C1+', 1000024),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.96, 'xsec 13.0 TeV': 17.5}
    },
    {
        'ID' : 12,
        'SMS' : '(PV > su_R(1),gluino(2)), (su_R(1) > N1,q), (gluino(2) > t-,b,C1+(7)), (C1+(7) > N1,W+)',
        'Masses (GeV)' : [
            ('su_R', 991.3),
            ('gluino', 865.0),
            ('N1', 129.0),
            ('C1+', 269.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_R', 2000002),
            ('gluino', 1000021),
            ('N1', 1000022),
            ('C1+', 1000024),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.98, 'xsec 13.0 TeV': 17.8}
    },
    {
        'ID' : 13,
        'SMS' : '(PV > su_R(1),gluino(2)), (su_R(1) > N1,q), (gluino(2) > b,t+,C1-(7)), (C1-(7) > N1~,W-)',
        'Masses (GeV)' : [
            ('su_R', 991.3),
            ('gluino', 865.0),
            ('N1', 129.0),
            ('C1-', 269.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('su_R', 2000002),
            ('gluino', 1000021),
            ('N1', 1000022),
            ('C1-', -1000024),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.98, 'xsec 13.0 TeV': 17.8}
    },
    {
        'ID' : 14,
        'SMS' : '(PV > gluino(1),su_L(2)), (gluino(1) > N1,q,q), (su_L(2) > q,C1+(7)), (C1+(7) > N1,W+)',
        'Masses (GeV)' : [
            ('gluino', 865.0),
            ('su_L', 991.5),
            ('N1', 129.0),
            ('C1+', 269.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('gluino', 1000021),
            ('su_L', 1000002),
            ('N1', 1000022),
            ('C1+', 1000024),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 3.69, 'xsec 13.0 TeV': 30.1}
    },
    {
        'ID' : 15,
        'SMS' : '(PV > gluino(1),su_L(2)), (gluino(1) > N1,c,c), (su_L(2) > q,C1+(7)), (C1+(7) > N1,W+)',
        'Masses (GeV)' : [
            ('gluino', 865.0),
            ('su_L', 991.5),
            ('N1', 129.0),
            ('C1+', 269.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('gluino', 1000021),
            ('su_L', 1000002),
            ('N1', 1000022),
            ('C1+', 1000024),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.23, 'xsec 13.0 TeV': 10.0}
    },
    {
        'ID' : 16,
        'SMS' : '(PV > su_L(1),su_L(2)), (su_L(1) > q,C1+(4)), (su_L(2) > q,C1+(6)), (C1+(4) > N1,W+), (C1+(6) > N1,W+)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('su_L', 991.5),
            ('C1+', 269.0),
            ('C1+', 269.0),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('su_L', 1000002),
            ('C1+', 1000024),
            ('C1+', 1000024),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 7.9, 'xsec 13.0 TeV': 56.3}
    },
    {
        'ID' : 17,
        'SMS' : '(PV > su_L(1),sd_L(2)), (su_L(1) > q,C1+(4)), (sd_L(2) > q,C1-(6)), (C1+(4) > N1,W+), (C1-(6) > N1~,W-)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('sd_L', 994.3),
            ('C1+', 269.0),
            ('C1-', 269.0),
            ('N1', 129.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('sd_L', 1000001),
            ('C1+', 1000024),
            ('C1-', -1000024),
            ('N1', 1000022),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 2.4, 'xsec 13.0 TeV': 20.5}
    },
    {
        'ID' : 18,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,N2(4)), (gluino(2) > q,q,C1-(7)), (N2(4) > N1,higgs), (C1-(7) > N1~,W-)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('N2', 268.9),
            ('C1-', 269.0),
            ('N1', 129.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('N2', 1000023),
            ('C1-', -1000024),
            ('N1', 1000022),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 57.6, 'xsec 13.0 TeV': 469.0}
    },
    {
        'ID' : 19,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,N2(4)), (gluino(2) > q,q,C1+(7)), (N2(4) > N1,higgs), (C1+(7) > N1,W+)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('N2', 268.9),
            ('C1+', 269.0),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('N2', 1000023),
            ('C1+', 1000024),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.23, 'xsec 13.0 TeV': 10.1}
    },
    {
        'ID' : 20,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,N2(4)), (gluino(2) > q,c,C1-(7)), (N2(4) > N1,higgs), (C1-(7) > N1~,W-)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('N2', 268.9),
            ('C1-', 269.0),
            ('N1', 129.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('N2', 1000023),
            ('C1-', -1000024),
            ('N1', 1000022),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.23, 'xsec 13.0 TeV': 10.1}
    },
    {
        'ID' : 21,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,N2(4)), (gluino(2) > q,c,C1+(7)), (N2(4) > N1,higgs), (C1+(7) > N1,W+)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('N2', 268.9),
            ('C1+', 269.0),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('N2', 1000023),
            ('C1+', 1000024),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.23, 'xsec 13.0 TeV': 10.1}
    },
    {
        'ID' : 22,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,N2(4)), (gluino(2) > t-,b,C1+(7)), (N2(4) > N1,higgs), (C1+(7) > N1,W+)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('N2', 268.9),
            ('C1+', 269.0),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('N2', 1000023),
            ('C1+', 1000024),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.25, 'xsec 13.0 TeV': 10.2}
    },
    {
        'ID' : 23,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,N2(4)), (gluino(2) > b,t+,C1-(7)), (N2(4) > N1,higgs), (C1-(7) > N1~,W-)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('N2', 268.9),
            ('C1-', 269.0),
            ('N1', 129.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('N2', 1000023),
            ('C1-', -1000024),
            ('N1', 1000022),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.25, 'xsec 13.0 TeV': 10.2}
    },
    {
        'ID' : 24,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,C1+(4)), (gluino(2) > q,q,N2(7)), (C1+(4) > N1,W+), (N2(7) > N1,higgs)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('C1+', 269.0),
            ('N2', 268.9),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('C1+', 1000024),
            ('N2', 1000023),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.26, 'xsec 13.0 TeV': 10.3}
    },
    {
        'ID' : 25,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,C1+(4)), (gluino(2) > q,q,C1-(7)), (C1+(4) > N1,W+), (C1-(7) > N1~,W-)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('C1+', 269.0),
            ('C1-', 269.0),
            ('N1', 129.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('C1+', 1000024),
            ('C1-', -1000024),
            ('N1', 1000022),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 2.96, 'xsec 13.0 TeV': 24.1}
    },
    {
        'ID' : 26,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,C1+(4)), (gluino(2) > q,q,C1+(7)), (C1+(4) > N1,W+), (C1+(7) > N1,W+)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('C1+', 269.0),
            ('C1+', 269.0),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('C1+', 1000024),
            ('C1+', 1000024),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 2.96, 'xsec 13.0 TeV': 24.1}
    },
    {
        'ID' : 27,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,C1+(4)), (gluino(2) > q,c,C1-(7)), (C1+(4) > N1,W+), (C1-(7) > N1~,W-)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('C1+', 269.0),
            ('C1-', 269.0),
            ('N1', 129.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('C1+', 1000024),
            ('C1-', -1000024),
            ('N1', 1000022),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 2.96, 'xsec 13.0 TeV': 24.1}
    },
    {
        'ID' : 28,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,C1+(4)), (gluino(2) > q,c,C1+(7)), (C1+(4) > N1,W+), (C1+(7) > N1,W+)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('C1+', 269.0),
            ('C1+', 269.0),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('C1+', 1000024),
            ('C1+', 1000024),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 2.96, 'xsec 13.0 TeV': 24.1}
    },
    {
        'ID' : 29,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,C1+(4)), (gluino(2) > c,c,N2(7)), (C1+(4) > N1,W+), (N2(7) > N1,higgs)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('C1+', 269.0),
            ('N2', 268.9),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('C1+', 1000024),
            ('N2', 1000023),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.26, 'xsec 13.0 TeV': 10.3}
    },
    {
        'ID' : 30,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,C1+(4)), (gluino(2) > t-,b,C1+(7)), (C1+(4) > N1,W+), (C1+(7) > N1,W+)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('C1+', 269.0),
            ('C1+', 269.0),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('C1+', 1000024),
            ('C1+', 1000024),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 3.0, 'xsec 13.0 TeV': 24.5}
    },
    {
        'ID' : 31,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,C1+(4)), (gluino(2) > b,b,N2(7)), (C1+(4) > N1,W+), (N2(7) > N1,higgs)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('C1+', 269.0),
            ('N2', 268.9),
            ('N1', 129.0),
            ('N1', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('C1+', 1000024),
            ('N2', 1000023),
            ('N1', 1000022),
            ('N1', 1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 1.6, 'xsec 13.0 TeV': 13.0}
    },
    {
        'ID' : 32,
        'SMS' : '(PV > su_L(1),gluino(2)), (su_L(1) > q,C1+(4)), (gluino(2) > b,t+,C1-(7)), (C1+(4) > N1,W+), (C1-(7) > N1~,W-)',
        'Masses (GeV)' : [
            ('su_L', 991.5),
            ('gluino', 865.0),
            ('C1+', 269.0),
            ('C1-', 269.0),
            ('N1', 129.0),
            ('N1~', 129.0)
        ],
        'PIDs' : [
            ('su_L', 1000002),
            ('gluino', 1000021),
            ('C1+', 1000024),
            ('C1-', -1000024),
            ('N1', 1000022),
            ('N1~', -1000022)
        ],
        'Weights (fb)' : {'xsec 8.0 TeV': 3.0, 'xsec 13.0 TeV': 24.5}
    }
],
'ExptRes' : [
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.1797601,
        'upper limit (fb)' : 0.0582,
        'expected upper limit (fb)' : 0.0582,
        'TxNames' : ['T6WW'],
        'FinalStates' : ['PV > (W,jet,MET),(W,jet,MET)'],
        'Mass (GeV)' : [
            ('sd_L/su_L', 993.11),
            ('su_L', 991.5),
            ('C1+/C1-', 269.0),
            ('C1+', 269.0),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2018-10',
        'DataSetID' : '4j0blowx_3',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'efficiencyMap',
        'r' : 3.088661,
        'r_expected' : 3.088661,
        'Width (GeV)' : [
            ('sd_L/su_L', 14.817),
            ('su_L', 14.771),
            ('C1+/C1-', 0.0013847),
            ('C1+', 0.0013847),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'T6WW': 0.17976008654165604},
        'nll' : 15.49185,
        'nll_min' : 3.183696,
        'nll_SM' : 3.183696
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 2.184581,
        'upper limit (fb)' : 1.27,
        'expected upper limit (fb)' : 1.07,
        'TxNames' : ['T6WW'],
        'FinalStates' : ['PV > (W,jet,MET),(W,jet,MET)'],
        'Mass (GeV)' : [
            ('sd_L/su_L', 993.11),
            ('su_L', 991.5),
            ('C1+/C1-', 269.0),
            ('C1+', 269.0),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2016-07',
        'DataSetID' : '5j_Meff_1600',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 36.1,
        'dataType' : 'efficiencyMap',
        'r' : 1.720142,
        'r_expected' : 2.041664,
        'Width (GeV)' : [
            ('sd_L/su_L', 14.817),
            ('su_L', 14.771),
            ('C1+/C1-', 0.0013847),
            ('C1+', 0.0013847),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'T6WW': 2.184580974615298},
        'nll' : 14.25994,
        'nll_min' : 6.930189,
        'nll_SM' : 7.004632
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 35.78982,
        'upper limit (fb)' : 53.7511,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T6WW'],
        'FinalStates' : ['PV > (W,jet,MET),(W,jet,MET)'],
        'Mass (GeV)' : [
            ('sd_L/su_L', 993.11),
            ('su_L', 991.5),
            ('C1+/C1-', 269.0),
            ('C1+', 269.0),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2016-07',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 36.1,
        'dataType' : 'upperLimit',
        'r' : 0.6658436,
        'r_expected' : None,
        'Width (GeV)' : [
            ('sd_L/su_L', 14.817),
            ('su_L', 14.771),
            ('C1+/C1-', 0.0013847),
            ('C1+', 0.0013847),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'T6WW': 35.78982375338465}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 5.500066,
        'upper limit (fb)' : 11.5698,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T6WW'],
        'FinalStates' : ['PV > (W,jet,MET),(W,jet,MET)'],
        'Mass (GeV)' : [
            ('sd_L/su_L', 992.72),
            ('su_L', 991.5),
            ('C1+/C1-', 269.0),
            ('C1+', 269.0),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2013-20',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.0,
        'dataType' : 'upperLimit',
        'r' : 0.4753812,
        'r_expected' : None,
        'Width (GeV)' : [
            ('sd_L/su_L', 14.806),
            ('su_L', 14.771),
            ('C1+/C1-', 0.0013847),
            ('C1+', 0.0013847),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'T6WW': 5.500065906194978}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 5.500066,
        'upper limit (fb)' : 17.1895,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T6WW'],
        'FinalStates' : ['PV > (W,jet,MET),(W,jet,MET)'],
        'Mass (GeV)' : [
            ('sd_L/su_L', 992.72),
            ('su_L', 991.5),
            ('C1+/C1-', 269.0),
            ('C1+', 269.0),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2013-02',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 0.3199666,
        'r_expected' : None,
        'Width (GeV)' : [
            ('sd_L/su_L', 14.806),
            ('su_L', 14.771),
            ('C1+/C1-', 0.0013847),
            ('C1+', 0.0013847),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'T6WW': 5.500065906194978}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 144.5845,
        'upper limit (fb)' : 554.995,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWW'],
        'FinalStates' : ['PV > (W,MET),(W,MET)'],
        'Mass (GeV)' : [
            ('C1-', 269.0),
            ('C1+', 269.0),
            ('N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2018-32',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.2605149,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1-', 0.00138466665),
            ('C1+', 0.00138466665),
            ('N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWW': 144.584475}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.0146227,
        'upper limit (fb)' : 0.0589,
        'expected upper limit (fb)' : 0.0925,
        'TxNames' : ['TChiWW'],
        'FinalStates' : ['PV > (W,MET),(W,MET)'],
        'Mass (GeV)' : [
            ('C1-', 269.0),
            ('C1+', 269.0),
            ('N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2018-32',
        'DataSetID' : 'SRDF_0d_cuts',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'efficiencyMap',
        'r' : 0.2482631,
        'r_expected' : 0.1580832,
        'Width (GeV)' : [
            ('C1-', 0.00138466665),
            ('C1+', 0.00138466665),
            ('N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWW': 0.0146226954636},
        'nll' : 5.377263,
        'nll_min' : 4.917842,
        'nll_SM' : 4.917842
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 980.798,
        'expected upper limit (fb)' : 811.908,
        'TxNames' : ['TChiWH'],
        'FinalStates' : ['PV > (higgs,MET),(W,MET)'],
        'Mass (GeV)' : [('N2', 268.9), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2018-23',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.2111855,
        'r_expected' : 0.2551154,
        'Width (GeV)' : [
            ('N2', 0.00112364506),
            ('C1+', 0.00138466665),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 194.632,
        'expected upper limit (fb)' : 243.047,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2019-09',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.2050081,
        'r_expected' : 0.1641704,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.01354272,
        'upper limit (fb)' : 0.0672,
        'expected upper limit (fb)' : 0.1072,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2019-09',
        'DataSetID' : 'SRWZ_5',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'efficiencyMap',
        'r' : 0.2015286,
        'r_expected' : 0.1263314,
        'Width (GeV)' : [
            ('C1+/C1-', 0.0013847),
            ('N2', 0.0011236),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.013542722171439512},
        'nll' : 7.453417,
        'nll_min' : 7.047479,
        'nll_SM' : 7.047479
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.2640166,
        'upper limit (fb)' : 1.43,
        'expected upper limit (fb)' : 0.934,
        'TxNames' : ['T6WW', 'TChiWW', 'TChiWZ'],
        'FinalStates' : [
            'PV > (W,jet,MET),(W,jet,MET)',
            'PV > (W,MET),(W,MET)',
            'PV > (W,MET),(Z,MET)'
        ],
        'Mass (GeV)' : None,
        'AnalysisID' : 'CMS-SUS-13-012',
        'DataSetID' : '3NJet6_1500HTinf_300MHTinf',
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 19.5,
        'dataType' : 'efficiencyMap',
        'r' : 0.184627,
        'r_expected' : 0.282673,
        'Width (GeV)' : None,
        'TxNames weights (fb)' : {
            'T6WW' : 0.2624601049703651,
            'TChiWW' : 0.001204890709847,
            'TChiWZ' : 0.0003515868180097819
        },
        'nll' : 5.550729,
        'nll_min' : 5.429653,
        'nll_SM' : 6.052394
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.009064725,
        'upper limit (fb)' : 0.04949,
        'expected upper limit (fb)' : 0.05109,
        'TxNames' : ['TChiWH'],
        'FinalStates' : ['PV > (higgs,MET),(W,MET)'],
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2019-08',
        'DataSetID' : 'SR_LM_High_MCT',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'efficiencyMap',
        'r' : 0.1831628,
        'r_expected' : 0.1774266,
        'Width (GeV)' : [
            ('N2', 0.0011236),
            ('C1+/C1-', 0.0013847),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 0.009064724873000252},
        'nll' : 3.349255,
        'nll_min' : 3.23298,
        'nll_SM' : 3.23298
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 1148.34,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'FinalStates' : ['PV > (higgs,MET),(W,MET)'],
        'Mass (GeV)' : [('N2', 268.9), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2019-08',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.1803736,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 0.00112364506),
            ('C1+', 0.00138466665),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 1572.24,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'FinalStates' : ['PV > (higgs,MET),(W,MET)'],
        'Mass (GeV)' : [('N2', 268.9), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-17-004',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.1317421,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 0.00112364506),
            ('C1+', 0.00138466665),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 1626.51,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'FinalStates' : ['PV > (higgs,MET),(W,MET)'],
        'Mass (GeV)' : [('N2', 268.9), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-18-007',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 77.5,
        'dataType' : 'upperLimit',
        'r' : 0.1273464,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 0.00112364506),
            ('C1+', 0.00138466665),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 367.751,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-17-004',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.1085004,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 2115.13,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'FinalStates' : ['PV > (higgs,MET),(W,MET)'],
        'Mass (GeV)' : [('N2', 268.9), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-16-039',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.09792791,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 0.00112364506),
            ('C1+', 0.00138466665),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.01495825,
        'upper limit (fb)' : 0.1572,
        'expected upper limit (fb)' : 0.1934,
        'TxNames' : ['TChiWW', 'TChiWZ', 'TChiWH'],
        'FinalStates' : [
            'PV > (W-,MET),(W+,MET)',
            'PV > (W,MET),(Z,MET)',
            'PV > (higgs,MET),(W,MET)'
        ],
        'Mass (GeV)' : None,
        'AnalysisID' : 'CMS-SUS-21-002',
        'DataSetID' : 'b_veto_SR0',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 137.0,
        'dataType' : 'efficiencyMap',
        'r' : 0.09515423,
        'r_expected' : 0.07734357,
        'Width (GeV)' : None,
        'TxNames weights (fb)' : {
            'TChiWW' : 0.010798855394827498,
            'TChiWZ' : 0.0021157892362806815,
            'TChiWH' : 0.002043601101045336
        },
        'nll' : 6.462508,
        'nll_min' : 6.376614,
        'nll_SM' : 6.376614
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 467.97,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-16-039',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.08526429,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 2437.37,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'FinalStates' : ['PV > (higgs,MET),(W,MET)'],
        'Mass (GeV)' : [('N2', 268.9), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-16-045',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.08498105,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 0.00112364506),
            ('C1+', 0.00138466665),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 95.88424,
        'upper limit (fb)' : 1136.02,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'FinalStates' : ['PV > (higgs,MET),(W,MET)'],
        'Mass (GeV)' : [('N2', 268.9), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-13-006',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 19.5,
        'dataType' : 'upperLimit',
        'r' : 0.08440365,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 0.00112364506),
            ('C1+', 0.00138466665),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 95.88423927189842}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.04476938,
        'upper limit (fb)' : 0.531,
        'expected upper limit (fb)' : 0.438,
        'TxNames' : ['TChiWW'],
        'FinalStates' : ['PV > (W-,MET),(W+,MET)'],
        'Mass (GeV)' : [
            ('C1-', 269.0),
            ('C1+', 269.0),
            ('N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2013-11',
        'DataSetID' : 'WWc-DF',
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'efficiencyMap',
        'r' : 0.08431146,
        'r_expected' : 0.1022132,
        'Width (GeV)' : [
            ('C1-', 0.00138466665),
            ('C1+', 0.00138466665),
            ('N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWW': 0.044769383356433996},
        'nll' : 3.871663,
        'nll_min' : 3.832856,
        'nll_SM' : 3.966765
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 2487.7,
        'expected upper limit (fb)' : 3375.71,
        'TxNames' : ['TChiWH'],
        'FinalStates' : ['PV > (higgs,MET),(W,MET)'],
        'Mass (GeV)' : [('N2', 268.9), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-21-002',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 137.0,
        'dataType' : 'upperLimit',
        'r' : 0.08326176,
        'r_expected' : 0.06135902,
        'Width (GeV)' : [
            ('N2', 0.00112364506),
            ('C1+', 0.00138466665),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.01256271,
        'upper limit (fb)' : 0.15899,
        'expected upper limit (fb)' : 0.14317,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2018-06',
        'DataSetID' : 'SR_low',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'efficiencyMap',
        'r' : 0.07901573,
        'r_expected' : 0.08774681,
        'Width (GeV)' : [
            ('C1+/C1-', 0.0013847),
            ('N2', 0.0011236),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.012562710565159365},
        'nll' : 5.485862,
        'nll_min' : 5.414862,
        'nll_SM' : 5.584257
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 15.24553,
        'upper limit (fb)' : 197.645,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T6WW'],
        'FinalStates' : ['PV > (W,jet,MET),(W,jet,MET)'],
        'Mass (GeV)' : [('su_L', 991.5), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2018-10',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.07713592,
        'r_expected' : None,
        'Width (GeV)' : [
            ('su_L', 14.7712529),
            ('C1+', 0.00138466665),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'T6WW': 15.245528607859374}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.007638752,
        'upper limit (fb)' : 0.1,
        'expected upper limit (fb)' : 0.124,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2016-24',
        'DataSetID' : 'WZ-0Jb',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 36.1,
        'dataType' : 'efficiencyMap',
        'r' : 0.07638752,
        'r_expected' : 0.06160284,
        'Width (GeV)' : [
            ('C1+/C1-', 0.0013847),
            ('N2', 0.0011236),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.007638751758050727},
        'nll' : 2.057516,
        'nll_min' : 1.884691,
        'nll_SM' : 1.884691
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.006656745,
        'upper limit (fb)' : 0.09,
        'expected upper limit (fb)' : 0.125,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2017-03',
        'DataSetID' : 'SR2l_Int',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 36.1,
        'dataType' : 'efficiencyMap',
        'r' : 0.07396383,
        'r_expected' : 0.05325396,
        'Width (GeV)' : [
            ('C1+/C1-', 0.0013847),
            ('N2', 0.0011236),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.006656745095367584},
        'nll' : 2.344793,
        'nll_min' : 2.219183,
        'nll_SM' : 2.219183
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 584.765,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2016-24',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 36.1,
        'dataType' : 'upperLimit',
        'r' : 0.06823447,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.01717407,
        'upper limit (fb)' : 0.256,
        'expected upper limit (fb)' : 0.291,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2013-12',
        'DataSetID' : 'SR0tau_a_Bin16',
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'efficiencyMap',
        'r' : 0.06708621,
        'r_expected' : 0.05901742,
        'Width (GeV)' : [
            ('C1+/C1-', 0.0013847),
            ('N2', 0.0011236),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.01717407016330772},
        'nll' : 3.232349,
        'nll_min' : 3.143051,
        'nll_SM' : 3.143051
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 601.75,
        'expected upper limit (fb)' : 391.04,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2018-06',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.06630848,
        'r_expected' : 0.1020385,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.006744727,
        'upper limit (fb)' : 0.107,
        'expected upper limit (fb)' : 0.0923,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2018-05-ewk',
        'DataSetID' : 'SRInt_1_cuts',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'efficiencyMap',
        'r' : 0.06303483,
        'r_expected' : 0.07307397,
        'Width (GeV)' : [
            ('C1+/C1-', 0.0013847),
            ('N2', 0.0011236),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.006744727082739748},
        'nll' : 4.684092,
        'nll_min' : 4.683139,
        'nll_SM' : 4.703294
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 18.47093,
        'upper limit (fb)' : 330.645,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2013-12',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 0.05586334,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 18.47093282810156}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.05091224,
        'upper limit (fb)' : 0.954862,
        'expected upper limit (fb)' : 0.931994,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-16-039-agg',
        'DataSetID' : 'AR1',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'efficiencyMap',
        'r' : 0.05331895,
        'r_expected' : 0.05462722,
        'Width (GeV)' : [
            ('C1+/C1-', 0.0013847),
            ('N2', 0.0011236),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.050912243359359},
        'nll' : 7.542752,
        'nll_min' : 7.542241,
        'nll_SM' : 7.542987
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 144.5845,
        'upper limit (fb)' : 3152.54,
        'expected upper limit (fb)' : 2945.62,
        'TxNames' : ['TChiWW'],
        'FinalStates' : ['PV > (W-,MET),(W+,MET)'],
        'Mass (GeV)' : [
            ('C1-', 269.0),
            ('C1+', 269.0),
            ('N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-21-002',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 137.0,
        'dataType' : 'upperLimit',
        'r' : 0.04586285,
        'r_expected' : 0.04908456,
        'Width (GeV)' : [
            ('C1-', 0.00138466665),
            ('C1+', 0.00138466665),
            ('N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWW': 144.584475}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 881.725,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2018-05',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.04525348,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 18.47093,
        'upper limit (fb)' : 465.626,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-13-006',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 19.5,
        'dataType' : 'upperLimit',
        'r' : 0.03966903,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 18.47093282810156}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 1141.88,
        'expected upper limit (fb)' : 867.141,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-20-001',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 137.0,
        'dataType' : 'upperLimit',
        'r' : 0.03494336,
        'r_expected' : 0.04601458,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 15.24553,
        'upper limit (fb)' : 520.288,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T6WW'],
        'FinalStates' : ['PV > (W,jet,MET),(W,jet,MET)'],
        'Mass (GeV)' : [('su_L', 991.5), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2018-22',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.0293021,
        'r_expected' : None,
        'Width (GeV)' : [
            ('su_L', 14.7712529),
            ('C1+', 0.00138466665),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'T6WW': 15.245528607859374}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 1571.83,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-16-034',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.02538514,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 18.47093,
        'upper limit (fb)' : 924.951,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2013-11',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 0.01996963,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 18.47093282810156}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 12712.3,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'FinalStates' : ['PV > (higgs,MET),(W,MET)'],
        'Mass (GeV)' : [('N2', 268.9), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-16-043',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.01629369,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 0.00112364506),
            ('C1+', 0.00138466665),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 2713.7,
        'expected upper limit (fb)' : 2526.8,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'CMS-SUS-21-002',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 137.0,
        'dataType' : 'upperLimit',
        'r' : 0.01470359,
        'r_expected' : 0.01579117,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 4177.63,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'FinalStates' : ['PV > (W,MET),(Z,MET)'],
        'Mass (GeV)' : [('C1+', 269.0), ('N2', 268.9), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2017-03',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 36.1,
        'dataType' : 'upperLimit',
        'r' : 0.00955114,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+', 0.00138466665),
            ('N2', 0.00112364506),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    }
],
'CombinedRes' : [
    {
        'AnalysisID' : 'ATLAS-SUSY-2018-10,CMS-SUS-13-012',
        'r' : 3.110038,
        'r_expected' : 3.169283,
        'nll' : 21.04258,
        'nll_min' : 9.236051,
        'nll_SM' : 9.23609,
        'Txnames' : ['T6WW', 'TChiWW', 'TChiWZ']
    }
],
'Total xsec for missing topologies (fb)' : 388.8044,
'missing topologies' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 96.48054,
        'SMS' : 'PV > (W,jet,MET),(W,jet,jet,MET)',
        'SMS IDs' : [25, 26, 27, 28]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 70.17688,
        'SMS' : 'PV > (jet,MET),(W,jet,jet,MET)',
        'SMS IDs' : [8, 9, 10, 11]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 48.96736,
        'SMS' : 'PV > (W,jet,MET),(W,b,t,MET)',
        'SMS IDs' : [30, 32]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 40.20174,
        'SMS' : 'PV > (higgs,jet,MET),(W,jet,jet,MET)',
        'SMS IDs' : [18, 19, 20, 21]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 35.6173,
        'SMS' : 'PV > (jet,MET),(W,b,t,MET)',
        'SMS IDs' : [12, 13]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.57984,
        'SMS' : 'PV > (W,jet,MET),(higgs,jet,jet,MET)',
        'SMS IDs' : [24, 29]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.40384,
        'SMS' : 'PV > (higgs,jet,MET),(W,b,t,MET)',
        'SMS IDs' : [22, 23]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.05211,
        'SMS' : 'PV > (jet,jet,MET),(W,jet,MET)',
        'SMS IDs' : [14, 15]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 13.04106,
        'SMS' : 'PV > (W,jet,MET),(b,b,higgs,MET)',
        'SMS IDs' : [31]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 11.82386,
        'SMS' : 'PV > (jet,MET),(W,jet,MET)',
        'SMS IDs' : [7]
    }
],
'Total xsec for missing topologies with displaced decays (fb)' : 0.0,
'missing topologies with displaced decays' : [],
'Total xsec for missing topologies with prompt decays (fb)' : 388.8044,
'missing topologies with prompt decays' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 96.48054,
        'SMS' : 'PV > (W,jet,MET),(W,jet,jet,MET)',
        'SMS IDs' : [25, 26, 27, 28]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 70.17688,
        'SMS' : 'PV > (jet,MET),(W,jet,jet,MET)',
        'SMS IDs' : [8, 9, 10, 11]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 48.96736,
        'SMS' : 'PV > (W,jet,MET),(W,b,t,MET)',
        'SMS IDs' : [30, 32]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 40.20174,
        'SMS' : 'PV > (higgs,jet,MET),(W,jet,jet,MET)',
        'SMS IDs' : [18, 19, 20, 21]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 35.6173,
        'SMS' : 'PV > (jet,MET),(W,b,t,MET)',
        'SMS IDs' : [12, 13]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.57984,
        'SMS' : 'PV > (W,jet,MET),(higgs,jet,jet,MET)',
        'SMS IDs' : [24, 29]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.40384,
        'SMS' : 'PV > (higgs,jet,MET),(W,b,t,MET)',
        'SMS IDs' : [22, 23]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.05211,
        'SMS' : 'PV > (jet,jet,MET),(W,jet,MET)',
        'SMS IDs' : [14, 15]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 13.04106,
        'SMS' : 'PV > (W,jet,MET),(b,b,higgs,MET)',
        'SMS IDs' : [31]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 11.82386,
        'SMS' : 'PV > (jet,MET),(W,jet,MET)',
        'SMS IDs' : [7]
    }
],
'Total xsec for topologies outside the grid (fb)' : 0.0,
'topologies outside the grid' : []
}
