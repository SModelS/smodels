smodelsOutput = {
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
    'combineanas' : 'ATLAS-SUSY-2018-10,CMS-SUS-13-012',
    'reportallsrs' : False,
    'experimentalfeatures' : False,
    'file status' : 1,
    'decomposition status' : 1,
    'warnings' : 'Input file ok',
    'input file' : 'inputFiles/slha/gluino_squarks.slha',
    'database version' : '3.0.0-beta',
    'smodels version' : '3.0.0-beta'
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
        'Weights (fb)' : {'xsec 8.0 TeV': 174.0, 'xsec 13.0 TeV': 392.0}
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
        'Weights (fb)' : {'xsec 8.0 TeV': 3.18, 'xsec 13.0 TeV': 23.3}
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
        'Weights (fb)' : {'xsec 8.0 TeV': 11.8, 'xsec 13.0 TeV': 106.0}
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
        'Weights (fb)' : {'xsec 8.0 TeV': 2.46, 'xsec 13.0 TeV': 20.1}
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
        'Weights (fb)' : {'xsec 8.0 TeV': 5.5, 'xsec 13.0 TeV': 35.8}
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
        'Weights (fb)' : {'xsec 8.0 TeV': 29.4, 'xsec 13.0 TeV': 240.0}
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
            ('sd_L/su_L', 'prompt'),
            ('su_L', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('C1+', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'T6WW': 0.17976008654165604},
        'nll' : 15.491848,
        'nll_min' : 3.183696,
        'nll_SM' : 3.183696
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 2.184581,
        'upper limit (fb)' : 1.27,
        'expected upper limit (fb)' : 1.07,
        'TxNames' : ['T6WW'],
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
            ('sd_L/su_L', 'prompt'),
            ('su_L', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('C1+', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'T6WW': 2.184580974615298},
        'likelihood' : 6.411921e-07,
        'l_max' : 0.000977816,
        'l_SM' : 0.0009076683
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 35.78982,
        'upper limit (fb)' : 53.7511,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T6WW'],
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
            ('sd_L/su_L', 'prompt'),
            ('su_L', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('C1+', 'prompt'),
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
            ('sd_L/su_L', 'prompt'),
            ('su_L', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('C1+', 'prompt'),
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
            ('sd_L/su_L', 'prompt'),
            ('su_L', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('C1+', 'prompt'),
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
            ('C1-', 'prompt'),
            ('C1+', 'prompt'),
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
            ('C1-', 'prompt'),
            ('C1+', 'prompt'),
            ('N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWW': 0.0146226954636},
        'likelihood' : 0.004620449,
        'l_max' : 0.0073149,
        'l_SM' : 0.0073149
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 980.798,
        'expected upper limit (fb)' : 811.908,
        'TxNames' : ['TChiWH'],
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2018-23',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.2111855,
        'r_expected' : 0.2551154,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 194.632,
        'expected upper limit (fb)' : 243.047,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2019-09',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.2050081,
        'r_expected' : 0.1641704,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.01354272,
        'upper limit (fb)' : 0.0713,
        'expected upper limit (fb)' : 0.146,
        'TxNames' : ['TChiWZ'],
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
        'r' : 0.18994,
        'r_expected' : 0.09275837,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.013542722171439512},
        'likelihood' : 0.0005794582,
        'l_max' : 0.0008695981,
        'l_SM' : 0.0008695981
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.2640824,
        'upper limit (fb)' : 1.43,
        'expected upper limit (fb)' : 0.934,
        'TxNames' : ['T6WW', 'TChiWW', 'TChiWZ'],
        'Mass (GeV)' : None,
        'AnalysisID' : 'CMS-SUS-13-012',
        'DataSetID' : '3NJet6_1500HTinf_300MHTinf',
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 19.5,
        'dataType' : 'efficiencyMap',
        'r' : 0.184673,
        'r_expected' : 0.2827435,
        'Width (GeV)' : None,
        'TxNames weights (fb)' : {
            'T6WW' : 0.2624601049703651,
            'TChiWW' : 0.001204890709847,
            'TChiWZ' : 0.0004174098342360046
        },
        'likelihood' : 0.003884919,
        'l_max' : 0.004384619,
        'l_SM' : 0.002352225
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 1148.34,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2019-08',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.1803736,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 1148.34,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2019-08-grp',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.1803736,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.009064725,
        'upper limit (fb)' : 0.0531,
        'expected upper limit (fb)' : 0.0531,
        'TxNames' : ['TChiWH'],
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
        'r' : 0.1707104,
        'r_expected' : 0.1707104,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 0.009064724873000252},
        'likelihood' : 0.03511049,
        'l_max' : 0.03943978,
        'l_SM' : 0.03943978
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 1572.24,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'CMS-PAS-SUS-17-004',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.1317421,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 1572.24,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-17-004',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.1317421,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 1626.39,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-18-007',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 77.5,
        'dataType' : 'upperLimit',
        'r' : 0.1273558,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 367.751,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'CMS-PAS-SUS-17-004',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.1085004,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 367.751,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-17-004',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.1085004,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
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
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-16-039',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.09792791,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.01495825,
        'upper limit (fb)' : 0.1572,
        'expected upper limit (fb)' : 0.1934,
        'TxNames' : ['TChiWH', 'TChiWW', 'TChiWZ'],
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
        'likelihood' : 0.001560876,
        'l_max' : 0.001700873,
        'l_SM' : 0.001700873
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.01448299,
        'upper limit (fb)' : 0.15899,
        'expected upper limit (fb)' : 0.14317,
        'TxNames' : ['TChiWZ'],
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
        'r' : 0.09109373,
        'r_expected' : 0.1011594,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.014482992254171961},
        'likelihood' : 0.004191917,
        'l_max' : 0.004449953,
        'l_SM' : 0.003756542
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 467.967,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-16-039',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.08526483,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 2440.76,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-16-045',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.08486302,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 95.88424,
        'upper limit (fb)' : 1136.02,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWH'],
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-13-006',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 19.5,
        'dataType' : 'upperLimit',
        'r' : 0.08440365,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 95.88423927189842}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 207.1303,
        'upper limit (fb)' : 2487.7,
        'expected upper limit (fb)' : 3375.71,
        'TxNames' : ['TChiWH'],
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-21-002',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 137.0,
        'dataType' : 'upperLimit',
        'r' : 0.08326176,
        'r_expected' : 0.06135902,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.04389609,
        'upper limit (fb)' : 0.531,
        'expected upper limit (fb)' : 0.438,
        'TxNames' : ['TChiWW'],
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
        'r' : 0.08266684,
        'r_expected' : 0.1002194,
        'Width (GeV)' : [
            ('C1-', 'prompt'),
            ('C1+', 'prompt'),
            ('N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWW': 0.043896094537829604},
        'likelihood' : 0.02079683,
        'l_max' : 0.02164771,
        'l_SM' : 0.01893459
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 15.24553,
        'upper limit (fb)' : 197.645,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T6WW'],
        'Mass (GeV)' : [('su_L', 991.5), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2018-10',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.07713592,
        'r_expected' : None,
        'Width (GeV)' : [
            ('su_L', 'prompt'),
            ('C1+', 'prompt'),
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
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.007638751758050727},
        'likelihood' : 0.127771,
        'l_max' : 0.151876,
        'l_SM' : 0.151876
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 584.765,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2016-24',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 36.1,
        'dataType' : 'upperLimit',
        'r' : 0.06823447,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
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
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.01717407016330772},
        'likelihood' : 0.03946467,
        'l_max' : 0.04315095,
        'l_SM' : 0.04315095
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 18.47093,
        'upper limit (fb)' : 281.789,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'CMS-PAS-SUS-12-022',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 9.2,
        'dataType' : 'upperLimit',
        'r' : 0.06554881,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 18.47093282810156}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.006744727,
        'upper limit (fb)' : 0.107,
        'expected upper limit (fb)' : 0.0923,
        'TxNames' : ['TChiWZ'],
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
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.006744727082739748},
        'likelihood' : 0.009241119,
        'l_max' : 0.009249933,
        'l_SM' : 0.009065371
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 18.47093,
        'upper limit (fb)' : 330.645,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2013-12',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 0.05586334,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
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
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.050912243359359},
        'likelihood' : 0.0005299372,
        'l_max' : 0.0005302082,
        'l_SM' : 0.0005298125
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 18.47093,
        'upper limit (fb)' : 354.511,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-CONF-2013-035',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.7,
        'dataType' : 'upperLimit',
        'r' : 0.05210257,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 18.47093282810156}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 825.36,
        'expected upper limit (fb)' : 519.9,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2018-06',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.04834391,
        'r_expected' : 0.0767477,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 144.5845,
        'upper limit (fb)' : 3152.54,
        'expected upper limit (fb)' : 2934.9,
        'TxNames' : ['TChiWW'],
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
        'r_expected' : 0.04926385,
        'Width (GeV)' : [
            ('C1-', 'prompt'),
            ('C1+', 'prompt'),
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
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2018-05',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.04525348,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 18.47093,
        'upper limit (fb)' : 464.427,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-13-006',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 19.5,
        'dataType' : 'upperLimit',
        'r' : 0.03977144,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 18.47093282810156}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 1141.22,
        'expected upper limit (fb)' : 866.225,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-20-001',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 137.0,
        'dataType' : 'upperLimit',
        'r' : 0.03496357,
        'r_expected' : 0.04606324,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.05030016,
        'upper limit (fb)' : 1.56,
        'expected upper limit (fb)' : 1.53,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-16-039',
        'DataSetID' : 'SR1',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'efficiencyMap',
        'r' : 0.03224369,
        'r_expected' : 0.03287592,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.05030016005491955},
        'likelihood' : 0.0005299514,
        'l_max' : 0.0005302082,
        'l_SM' : 0.0005298125
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 15.24553,
        'upper limit (fb)' : 520.288,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T6WW'],
        'Mass (GeV)' : [('su_L', 991.5), ('C1+', 269.0), ('N1', 129.0)],
        'AnalysisID' : 'ATLAS-SUSY-2018-22',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 139.0,
        'dataType' : 'upperLimit',
        'r' : 0.0293021,
        'r_expected' : None,
        'Width (GeV)' : [
            ('su_L', 'prompt'),
            ('C1+', 'prompt'),
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
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-16-034',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.02538514,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.01298155,
        'upper limit (fb)' : 0.53,
        'expected upper limit (fb)' : 0.26,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2017-03',
        'DataSetID' : 'SR3l_Low',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 36.1,
        'dataType' : 'efficiencyMap',
        'r' : 0.0244935,
        'r_expected' : 0.04992905,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 0.01298155270640454},
        'likelihood' : 0.00160016,
        'l_max' : 0.01772008,
        'l_SM' : 0.001214341
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 18.47093,
        'upper limit (fb)' : 924.951,
        'expected upper limit (fb)' : None,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2013-11',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 0.01996963,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
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
        'Mass (GeV)' : [
            ('N2', 268.9),
            ('C1+/C1-', 269.0),
            ('N1', 129.0),
            ('N1/N1~', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-16-043',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 35.9,
        'dataType' : 'upperLimit',
        'r' : 0.01629369,
        'r_expected' : None,
        'Width (GeV)' : [
            ('N2', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('N1', 'stable'),
            ('N1/N1~', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWH': 207.13027034006157}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 39.90113,
        'upper limit (fb)' : 2713.7,
        'expected upper limit (fb)' : 2511.59,
        'TxNames' : ['TChiWZ'],
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'CMS-SUS-21-002',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 137.0,
        'dataType' : 'upperLimit',
        'r' : 0.01470359,
        'r_expected' : 0.0158868,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
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
        'Mass (GeV)' : [
            ('C1+/C1-', 269.0),
            ('N2', 268.9),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-SUSY-2017-03',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 36.1,
        'dataType' : 'upperLimit',
        'r' : 0.00955114,
        'r_expected' : None,
        'Width (GeV)' : [
            ('C1+/C1-', 'prompt'),
            ('N2', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'TChiWZ': 39.9011280599384}
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 5.500066,
        'upper limit (fb)' : 40471.1,
        'expected upper limit (fb)' : None,
        'TxNames' : ['T6WW'],
        'Mass (GeV)' : [
            ('sd_L/su_L', 992.72),
            ('su_L', 991.5),
            ('C1+/C1-', 269.0),
            ('C1+', 269.0),
            ('N1/N1~', 129.0),
            ('N1', 129.0)
        ],
        'AnalysisID' : 'ATLAS-CONF-2013-089',
        'DataSetID' : None,
        'AnalysisSqrts (TeV)' : 8.0,
        'lumi (fb-1)' : 20.3,
        'dataType' : 'upperLimit',
        'r' : 0.0001359011,
        'r_expected' : None,
        'Width (GeV)' : [
            ('sd_L/su_L', 'prompt'),
            ('su_L', 'prompt'),
            ('C1+/C1-', 'prompt'),
            ('C1+', 'prompt'),
            ('N1/N1~', 'stable'),
            ('N1', 'stable')
        ],
        'TxNames weights (fb)' : {'T6WW': 5.500065906194978}
    }
],
'CombinedRes' : [
    {
        'AnalysisID' : 'ATLAS-SUSY-2018-10,CMS-SUS-13-012',
        'r' : 3.129969,
        'r_expected' : 3.169335,
        'nll' : 21.042502,
        'nll_min' : 9.236050,
        'nll_SM' : 9.236090
    }
],
'Total xsec for missing topologies (fb)' : 388.8044,
'missing topologies' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 96.48054,
        'SMS' : 'PV > (W,jet,MET), (W,jet,jet,MET)',
        'SMS IDs' : [25, 26, 27, 28]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 70.17688,
        'SMS' : 'PV > (jet,MET), (W,jet,jet,MET)',
        'SMS IDs' : [8, 9, 10, 11]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 48.96736,
        'SMS' : 'PV > (W,jet,MET), (W,b,t,MET)',
        'SMS IDs' : [30, 32]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 40.20174,
        'SMS' : 'PV > (higgs,jet,MET), (W,jet,jet,MET)',
        'SMS IDs' : [18, 19, 20, 21]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 35.6173,
        'SMS' : 'PV > (jet,MET), (W,b,t,MET)',
        'SMS IDs' : [12, 13]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.57984,
        'SMS' : 'PV > (W,jet,MET), (higgs,jet,jet,MET)',
        'SMS IDs' : [24, 29]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.40384,
        'SMS' : 'PV > (higgs,jet,MET), (W,b,t,MET)',
        'SMS IDs' : [22, 23]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.05211,
        'SMS' : 'PV > (jet,jet,MET), (W,jet,MET)',
        'SMS IDs' : [14, 15]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 13.04106,
        'SMS' : 'PV > (W,jet,MET), (b,b,higgs,MET)',
        'SMS IDs' : [31]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 11.82386,
        'SMS' : 'PV > (jet,MET), (W,jet,MET)',
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
        'SMS' : 'PV > (W,jet,MET), (W,jet,jet,MET)',
        'SMS IDs' : [25, 26, 27, 28]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 70.17688,
        'SMS' : 'PV > (jet,MET), (W,jet,jet,MET)',
        'SMS IDs' : [8, 9, 10, 11]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 48.96736,
        'SMS' : 'PV > (W,jet,MET), (W,b,t,MET)',
        'SMS IDs' : [30, 32]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 40.20174,
        'SMS' : 'PV > (higgs,jet,MET), (W,jet,jet,MET)',
        'SMS IDs' : [18, 19, 20, 21]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 35.6173,
        'SMS' : 'PV > (jet,MET), (W,b,t,MET)',
        'SMS IDs' : [12, 13]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.57984,
        'SMS' : 'PV > (W,jet,MET), (higgs,jet,jet,MET)',
        'SMS IDs' : [24, 29]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.40384,
        'SMS' : 'PV > (higgs,jet,MET), (W,b,t,MET)',
        'SMS IDs' : [22, 23]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 20.05211,
        'SMS' : 'PV > (jet,jet,MET), (W,jet,MET)',
        'SMS IDs' : [14, 15]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 13.04106,
        'SMS' : 'PV > (W,jet,MET), (b,b,higgs,MET)',
        'SMS IDs' : [31]
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 11.82386,
        'SMS' : 'PV > (jet,MET), (W,jet,MET)',
        'SMS IDs' : [7]
    }
],
'Total xsec for topologies outside the grid (fb)' : 0.0,
'topologies outside the grid' : []
}
