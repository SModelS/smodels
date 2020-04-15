#!/usr/bin/env python3

import glob, copy, subprocess, gzip, os, sys
from smodels.tools.physicsUnits import GeV, TeV, pb
from smodels.particlesLoader import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model

def cloneExclLine():
    import ROOT
    target = ROOT.gDirectory
    f = ROOT.TFile("orig.root","read")
    lines = [ "TGQ/obsExclusion_[[x, 0.0], [y, 0.0]];1", "TGQ/obsExclusion_[[x, 695.0], [y, 695.0]];1", "TGQ/obsExclusion_[[x, 995.0], [y, 995.0]];1" ]
    lines += [ "TGQ/expExclusion_[[x, 0.0], [y, 0.0]];1", "TGQ/expExclusion_[[x, 695.0], [y, 695.0]];1", "TGQ/expExclusion_[[x, 995.0], [y, 995.0]];1" ]
    objs = []
    for eline in lines:
        orig = f.Get( eline )
        tgq = orig.Clone()
        objs.append ( tgq )
        print ( "tgq", tgq, tgq.GetN() )
        #eline = eline.replace("obj","exp" )
        #orig = f.Get( eline )
        #tgq = orig.Clone()
        #objs.append ( tgq )
        #print ( "tgq", tgq, tgq.GetN() )
    g = ROOT.TFile("sms.root","update")
    g.mkdir ( "TGQ12" )
    g.cd( "TGQ12" )
    for o in objs:
        o.Write()
    g.Write()
    g.Close()
    f.Close()
    target.cd()
    #cmd = "cp sms.root orig.root"
    #subprocess.getoutput ( cmd )
    cmd = "cp new.root sms.root"
    subprocess.getoutput ( cmd )

def checkNewExclLine():
    import ROOT
    print ( "Now check excl line" )
    ROOT.gDirectory.cd()
    f = ROOT.TFile("new.root")
    tgq = f.Get("TGQ12")
    print ( "check")
    tgq.ls()
    f.Close()

def getEfficiencies( fromPickle=True, debug=False, writePickle=True ):
    if fromPickle and os.path.exists ( "efficiencies.pcl" ):
        import pickle
        f=open("efficiencies.pcl","rb")
        efficiencies = pickle.load(f)
        f.close()
        return efficiencies
    from smodels.experiment.databaseObj import Database
    db = Database ( "../../../../smodels-database/" )
    txnames = [ "T1", "T2", "TGQ", "T3GQ", "T5GQ" ]
    efficiencies = {}
    expRes = db.getExpResults( analysisIDs = "ATLAS-SUSY-2016-07",
                               txnames = txnames,
                               dataTypes = [ "efficiencyMap" ],
                               useNonValidated=True )[0]
    gluinos = range ( 400, 5010, 50 )
    squarks = range ( 400, 5010, 50 )
    neutralinos = [ 0, 695, 995 ]
    if debug:
        gluinos = [ 1000 ]
        squarks = [  900 ]
        neutralinos = [ 695 ]
    datasets = expRes.datasets
    for dataset in datasets:
        dataId = dataset.dataInfo.dataId
        efficiencies [ dataId ] = {}
        print ( "dataset", dataset )
        for mgluino in gluinos:
            for msquark in squarks:
                for mN in neutralinos:
                    efficiencies[dataId][(mgluino,msquark,mN)]={}
                    for txname in dataset.txnameList:
                        txn = txname.txName
                        mass = [[mgluino*GeV,mN*GeV],[msquark*GeV,mN*GeV]]
                        if txn == "T1":
                            mass = [[mgluino*GeV,mN*GeV],[mgluino*GeV,mN*GeV]]
                        if txn == "T2":
                            mass = [[msquark*GeV,mN*GeV],[msquark*GeV,mN*GeV]]
                        if txn == "T3GQ":
                            mass = [[msquark*GeV,mN*GeV],[mgluino*GeV,msquark*GeV,mN*GeV]]
                        if txn == "T5GQ":
                            mass = [[msquark*GeV,mgluino*GeV,mN*GeV],[mgluino*GeV,mN*GeV]]
                        eff = txname.getEfficiencyFor ( mass )
                        efficiencies[dataId][(mgluino,msquark,mN)][txn]=eff
    if writePickle:
        import pickle
        f=open("efficiencies.pcl","wb")
        pickle.dump(efficiencies,f)
        f.close()
    return efficiencies
                
def retrieveMassesFromFileName ( slhafile ):
    """ get the particle masses from the slha file name """
    tokens = slhafile.split("_" )
    gluino = int(tokens[3])
    squark = int(tokens[1])
    neutralino = int(tokens[2])
    return gluino,squark,neutralino

def getWeights( fromPickle=True, debug=True, writePickle=True ):
    if fromPickle and os.path.exists ( "weights.pcl" ):
        import pickle
        f=open("weights.pcl","rb")
        weights = pickle.load(f)
        f.close()
        return weights
    txnames = [ "T1", "T2", "TGQ", "T3GQ", "T5GQ" ]
    weights = {}
    from smodels.theory.decomposer import decompose
    slhadir = "%s/git/smodels-utils/slha/" % os.environ["HOME"]
    files = glob.glob("%s/TGQ12*slha" % slhadir )
    if len(files)==0:
        print ( "[createTGQ12] cannot find TGQ12*slha. Maybe untar tarball? Will try!" )
        cmd = "cd %s; tar xzvf TGQ12.tar.gz" % slhadir
        o = subprocess.getoutput ( cmd )
        print ( "[createTGQ12] %s" % o )
        files = glob.glob("%s/TGQ12*slha" % slhadir )
    if len(files)==0:
        print ( "[createTGQ12] still cannot find TGQ12*slha. Stopping" )
        sys.exit()
    if debug:
        files = files[:2]
    for f in files:
        mgluino,msquark,mN = retrieveMassesFromFileName ( f )
        weights[(mgluino,msquark,mN)]={}
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=f)
        decomp = decompose(model)
        for d in decomp:
            for e in d.elementList:
                txn=str(e).replace("c","q")
                if txn == "[[[q]],[[q]]]": 
                    txn = "T2"
                if txn in [ "[[[q,q]],[[q,q]]]" ]: 
                    txn = "T1"
                if txn in [ "[[[q,q]],[[q]]]", "[[[q]],[[q,q]]]" ]:
                    txn = "TGQ"
                if txn in [ '[[],[[q]]]', '[[],[[q,q]]]', '[[],[]]' ]:
                    # we can savely ignore these
                    continue
                wts = e.weight
                for w in wts:
                    if w.info.sqrts < 10*TeV:
                        continue
                    if not txn in weights[(mgluino,msquark,mN)]:
                        weights[(mgluino,msquark,mN)][txn]=0.
                    weights[(mgluino,msquark,mN)][txn]+=w.value.asNumber(pb)
                    # print ( "decomp", e, txn, w.value )
    if writePickle:
        import pickle
        f=open("weights.pcl","wb")
        pickle.dump(weights,f)
        f.close()
    return weights
                            
def writeHeaders():
    dirs = glob.glob("*Meff*")
    for d in dirs:
        tgq = "%s/TGQ.txt" % d
        tgqf = open ( tgq, "rt" )
        tgqlines = tgqf.readlines()
        tgqf.close()
        out = "%s/TGQ12.txt" % d
        outf = open ( out, "wt" )
        for line in tgqlines:
            line = line.replace("TGQ","TGQ12")
            if line.startswith("constraint:" ):
                line = "constraint: [[[q]],[[q,q]]]+[[[q]],[[q]]]+[[[q,q]],[[q,q]]]\n"
            if line.startswith("efficiencyMap:"):
                break
            outf.write ( line )
        # outf.write ( "comment: FIXME check if i didnt mix up gluino-squark\n" )
        outf.close()

def computeWeightedEfficiencies ( effs, weights ):
    """ given the efficiencies and the weights, compute weighted efficiencies """
    ret = {}
    for dataset,points in effs.items():
        ret[dataset]={}
        # print ( "dataset", dataset )
        for point,es in points.items():
            if point in weights:
                w = weights[point]
                wtot = sum ( w.values() )
                # print ( "point", point )
                #print ( " `- the weights are", w )
                #print ( " `- the effs are", es )
                eavg = 0.
                if wtot > 0.:
                    for txn,e in es.items():
                        wt = 0.
                        if txn in w:
                            wt = w[txn]
                        eavg += e * wt / wtot
                # print ( "point", point, "avg eff is", eavg )
                ret[dataset][point]=eavg
    return ret

def writeBody ( efficiencies ):
    """ write the actual efficiency maps (as opposed to the metadata)
        in */TGQ12.txt """
    for dataset,points in efficiencies.items():
        f=open("%s/TGQ12.txt" % dataset, "at" )
        f.write ( "efficiencyMap: [" )
        for ctr,(point,eff) in enumerate(points.items()):
            print ( "point", point, eff )
            mG,mQ,mlsp = point[0],point[1],point[2]
            massvec="[[%d*GeV,%d*GeV],[%d*GeV,%d*GeV]]" % ( mG, mlsp, mQ, mlsp )
            cm=",\n"
            if ctr == len(points)-1:
                cm=""
            f.write ( "[%s,%f]%s" % ( massvec, eff, cm ) )
        f.write ( "]\n" )

if __name__ == "__main__":
    cloneExclLine()
    #sys.exit()
    #checkNewExclLine()
    weights = getWeights( fromPickle=True, debug=False, writePickle=True )
    #print ( weights )
    effs = getEfficiencies( fromPickle=True, debug=False, writePickle=True )
    weffs = computeWeightedEfficiencies ( effs, weights )
    print ( "weffs", weffs )
    writeHeaders()
    writeBody ( weffs )
