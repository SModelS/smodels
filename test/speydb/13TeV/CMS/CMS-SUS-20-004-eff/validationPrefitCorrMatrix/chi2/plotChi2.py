#!/usr/bin/env python3

from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.theory import decomposer
from smodels.tools.theoryPredictionsCombiner import TheoryPredictionsCombiner
from smodels.theory.model import Model
from smodels.share.models.SMparticles import SMList
from smodels.share.models.mssm import BSMList
from smodels.experiment.databaseObj import Database
from smodels.tools.physicsUnits import fb, GeV
from smodels_utils.plotting import mpkitty as plt
import numpy as np

def fetch():
    dbpath = "../../../../../"
    database = Database( dbpath )
    dTypes = ["efficiencyMap"]
    anaids = [ 'CMS-SUS-20-004', 'CMS-SUS-20-004-slv1' ]
    dsids = [ 'all' ]
    results = database.getExpResults(analysisIDs=anaids,
            datasetIDs=dsids, dataTypes=dTypes, useNonValidated = True )
    print ( "results", results )
    ret = {}
    for r in results:
        ret[r.globalInfo.id]=r
    return ret

def getTheoryPrediction( res, slhafile ):
    model = Model(BSMparticles=BSMList, SMparticles=SMList)
    model.updateParticles(inputFile=slhafile)
    sigmacut = 0.005*fb
    mingap = 5.*GeV
    smstopos = decomposer.decompose(model, sigmacut, doCompress=True,
           doInvisible=True, minmassgap=mingap )
    ts = theoryPredictionsFor(res, smstopos,
        combinedResults=True, useBestDataset=False )
    return ts[0]

def createRange ( xrange ):
    return np.arange ( xrange["min"], xrange["max"], xrange["delta"] )

def computeChi2s( tp, xrange : dict  ):
    chi2s = {}
    for i in createRange ( xrange ):
        llhd = tp.likelihood ( i )
        chi2 = -2.*np.log ( llhd )
        chi2s[i]= chi2
    lmax = -2. * np.log ( tp.lmax() )
    chi2s["lmax"] = lmax
    return chi2s


def plot ( chi2v2, chi2v1, setup ):
    plt.clf()
    slhafile = setup["slhafile"]
    xrange = setup["xrange"]
    x = createRange ( xrange )
    slv1 = setup["SLv1"]
    slv2 = setup["SLv2"]
    min1 = chi2v1.pop ( "lmax" )
    min2 = chi2v2.pop ( "lmax" )
    #minChi2v1 = min ( chi2v1.values() )
    #print ( "min1",min1,minChi2v1 )
    # minChi2v2 = min ( chi2v2.values() )
    #print ( "min2",min2,minChi2v2 )
    valuesv1 = np.array ( list ( chi2v1.values() ) ) - min1 # minChi2v1
    valuesv2 = np.array ( list ( chi2v2.values() ) ) - min2 # minChi2v2
    if "fullx" in setup:
        fullx = setup["fullx"]
        fully = setup["fully"]
        plt.plot ( fullx, fully, c="black", label="Bill, full" )
    plt.plot ( x, slv1, c="green", label="Bill, SLv1" )
    plt.plot ( x, slv2, c="red", label="Bill, SLv2" )

    plt.plot ( chi2v1.keys(), valuesv1, c="darkgreen", linestyle="dashed", label="SModelS SLv1" )
    plt.plot ( chi2v2.keys(), valuesv2, c="darkred", linestyle="dashed", label="SModelS SLv2" )
    plt.legend()
    ax = plt.gca()
    ax.set_ylim ( [0,10.] )
    plt.title ( rf"$\Delta\chi^2$, {slhafile}" )
    plt.ylabel ( r"$\Delta\chi^2$")
    plt.xlabel ( r"signal strength $\mu$")
    plt.kittyPlot ( f"chi2_{slhafile}.png" )

def getSetup( i=1 ):
    bills = {}
    nan = float("nan")
    fullx = list ( np.linspace(0,1,16+1) )
    bills[1] = { "slhafile": "TChiHH_300_1_300_1.slha",
              "xrange": { "min": -.3, "max": 2.5, "delta": .1 }, "fullx": fullx }
    bills[1]["fully"] = [0.0009691, 0.0411193, 0.1329959, 0.2649894, 0.452793, 0.6973918, 0.9997511, 1.3608103, 1.7814781, 2.2626238, 2.8050761, 3.4096148, 4.0769662, 4.8079994, 5.6025722, 6.4611683, 7.3840528]
    bills[1]["SLv1"] = [10.703636741203411, 7.765331170861913, 5.429922872544864, 3.5606053885559277, 2.229789961884478, 1.260335686905023, 0.5948967019269276, 0.19015382682113113, 0.013105705092357312, 0.03818822415058776, 0.2452012565470909, 0.6178577895101967, 1.1427724452224197, 1.8087528690991803, 2.6062939522198576, 3.527213285077295, 4.564382902877043, 5.7115283794849745, 6.963076307657644, 8.314036267358688, 9.759908402400413, 11.296609810749032, 12.920415346775968, 14.627909443326814, 16.4159463718986, 18.281617259999905, 20.222222333821747, 22.235247645327092]
    bills[1]["SLv2"] =  [nan, 4.739163622381085, 3.170939369951469, 1.957615092913187, 1.0712456270721873, 0.47268278821661625, 0.1270341933289103, 0.001616844675908169, 0.06915466396418424, 0.3079452647667438, 0.7009525554437062, 1.2347611187945802, 1.8986918542215392, 2.6841252144388363, 3.5840022457906855, 4.592462125350124, 5.704578884902361, 6.916170479623162, 8.223659270284799, 9.623968888293149, 11.114446275596748, 12.692798566434703, 14.357036401540569, 16.10541593011604, 17.936374556951506, 19.8484600254057, 21.84025932253914, 23.910337915497934]

    bills[2] = { "slhafile": "TChiHH_750_1_750_1.slha",
              "xrange": { "min": -.8, "max": 4., "delta": .1 } }
    bills[2]["SLv1"] = [3.6740493314549667, 3.3003495620199885, 3.0071787346453505, 2.680146956659314, 2.3980122381012734, 2.1274027457113505, 1.9607380050316294, 1.2529243968812693, 0.9618908025073551, 0.7828665146406308, 0.6237094311207727, 0.4839111398625562, 0.362958594317746, 0.26033536542308866, 0.17552288615570433, 0.10800225460496904, 0.05725549891508308, 0.02276737431495235, 0.004027023592556134, 0.0005294824431985035, 0.011777258538586466, 0.03728190069270454, 0.07656511063618154, 0.12916008074961383, 0.19461238127018987, 0.272480944802993, 0.36233830977815273, 0.4641402295536352, 0.5785746690100382, 0.7050367437235252, 0.8428885220737641, 0.99155924820937, 1.1505341124862696, 1.319346135998643, 1.4975692439113573, 1.6848126809073563, 1.8807163209257567, 2.0849470525998015, 2.2971953515371126, 2.517172970760896, 2.7446102663934937, 2.9792545535691204, 3.220868451685476, 3.469228153609407, 3.724122616582008, 3.9853519198365177, 4.252726995151761, 4.526068109937086]
    bills[2]["SLv2"] = [nan, nan, 0.6141081753405615, 0.5334588852087734, 0.49818052165827, 0.5061825307853667, 0.18165910415987696, 0.06620007136390882, 0.0046273374704526304, 5.8842375750600695e-05, 0.006910564078708603, 0.025060304607450234, 0.05435953944387961, 0.09463345607477436, 0.14568208145433914, 0.20728178340877434, 0.27918750208579013, 0.36113550693724505, 0.45284661902678636, 0.5540297939235472, 0.6643856002441737, 0.7836101524909793, 0.9113983685183769, 1.0474473551419692, 1.1914591401718155, 1.3431432266722538, 1.502218517389423, 1.668763052116475, 1.8434262922949358, 2.025742563971704, 2.2152410886799885, 2.4115141290772044, 2.6142064066774253, 2.8230069412734906, 3.037642209447455, 3.2578707209925994, 3.483478560536298, 3.7142758884058367, 3.950093626505293, 4.190781252356345, 4.436204416664026, 4.6862432192878885, 4.940790655822667, 5.199751327326908, 5.4630401770399715, 5.730581602676892, 6.002308227983548, 6.278160479433808]
    bills[3] = { "slhafile": "TChiHH_450_1_450_1.slha",
              "xrange": { "min": -.1, "max": 3., "delta": .1 } }
    bills[3]["SLv1"] = [12.455313913615527, 8.466680957231063, 5.858570682170921, 3.8658519724903897, 2.3815319129254817, 1.3173325341952307, 0.6024819222031113, 0.18098047578507703, 0.008599491688471517, 0.05022558186800552, 0.27772814306936766, 0.668324639680975, 1.203352307802163, 1.8673506389454246, 2.6473745385791005, 3.5324765717631976, 4.513313163379422, 5.5818420191634175, 6.731087276894499, 7.954955010852274, 9.248087277628969, 10.605745071653587, 12.023713918004717, 13.49822689998274, 15.025901607050486, 16.603687750003218, 18.228937974537615, 19.899932716004884, 21.614398810517628, 23.370098736677875, 25.164988422430497]
    bills[3]["SLv2"] = [7.621163625503328, 5.959306284217405, 4.551005246900019, 3.3372145922760126, 2.312882332428984, 1.473881933358257, 0.8192274353651783, 0.3523070549940144, 0.07890020076092696, 0.0008607059305347775, 0.11118433987104481, 0.3957993189250999, 0.8381518069398339, 1.4222066054934714, 2.1335949052340766, 2.959821949862942, 3.8901911253993546, 4.915244280571017, 6.027176080030415, 7.2189774463665515, 8.484586640987999, 9.818681843207571, 11.216564595038335, 12.674064603774383, 14.187461720266384, 15.753534398536232, 17.370048308701826, 19.034288610860642, 20.74363975407246, 22.495730511223172, 24.288403676200858]
    return bills[i]

def main():
    res = fetch()
    nrs = [ 1,2 , 3 ]
    # nrs = [ 1]
    for i in nrs:
        setup = getSetup( i)
        slhafile = setup["slhafile"]
        print ( "\n[plotChi2] now processing", slhafile )
        tpV2 = getTheoryPrediction ( res["CMS-SUS-20-004"], slhafile )
        tpV1 = getTheoryPrediction ( res["CMS-SUS-20-004-slv1"], slhafile )
        xrange = setup["xrange"]
        chi2v2 = computeChi2s( tpV2, xrange )
        chi2v1 = computeChi2s( tpV1, xrange )
        plot ( chi2v2, chi2v1, setup )

if __name__ == "__main__":
    main()
