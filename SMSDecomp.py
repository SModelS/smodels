#Do LHE-based decomposition.
#lhefile = LHE file with pythia events
#W = dictionary with event weights
#nevts = number of generated events
def LHEdecomp(lhefile,W,nevts):
    import SMSmethods

    SMSTopList = []
    LHEfile = open(lhefile,"r")
    for iev in range(nevts):
##Read event    
        PList = SMSmethods.getNextEvent(LHEfile) 
##Get mother PDGs:
        momPDG = tuple(SMSmethods.getMom(PList))
#Get event weight list:
        weight = {}
        for k in W.keys(): 
            if W[k].has_key(momPDG):
                weight.update({k : W[k][momPDG]})
            else:
                print "LHEdecomp: Error getting weight"
                return False

#Get event topology    
        SMSTopListEv = SMSmethods.getEventTop(PList, weight)
    
#Add event topology to topology list:  
        for TopEv in SMSTopListEv:  
            SMSmethods.AddToList(TopEv,SMSTopList)
            
    return SMSTopList




#Do SLHA-based decomposition.
#slhafile = file with mass spectrum and branching ratios
#Xsec = dictionary with cross-sections for pair production
#sigcut = minimum sigma*BR to be generated
def SLHAdecomp(slhafile,Xsec,sigcut):
    import sys, os, copy

    workdir = os.getcwd() 
    pyslhadir = workdir + "/pyslha-1.4.3"
    sys.path.append(pyslhadir)
    import pyslha, SMSmethods
    from SMSHelpers import addunit


#Read SLHA file
    res = pyslha.readSLHAFile(slhafile)
    
#Get list of particles with maximum production cross-section above sigcut and maximum cross-sections
    Pdic = {}
    for k in Xsec.keys():
        for k2 in Xsec[k]:
            if Xsec[k][k2] > sigcut:
                for ip in range(2):
                    if not Pdic.has_key(k2[ip]) or Pdic[k2[ip]] < Xsec[k][k2]: 
                        Pdic.update({k2[ip] : Xsec[k][k2]})

#Get list of branching ratios for all particles:
    BRdic = {}
    for k in res[1].keys():
        brs = copy.deepcopy(res[1][abs(k)].decays)
        for i in range(len(brs)):
            brs[i].ids = [-x for x in brs[i].ids]
        BRdic.update({k : res[1][abs(k)].decays, -k : brs})
#Get mass list for all particles
    Massdic = {}
    for k in res[1].keys():
        if k: 
            Massdic.update({k : addunit(res[1][k].mass,'GeV'), -k : addunit(res[1][k].mass,'GeV')})
            
#Loop over all particles and generate all possible 1branch-elements with sigmamax*BR > sigcut
    ElList = []
    WList = []
    for ptc in Pdic.keys():
        NewEl = SMSmethods.BElement()
        NewEl.momID = [ptc,ptc]
        NewEl.masses.append(Massdic[ptc])
        weight = Pdic[ptc]
        ElList.append(NewEl)
        WList.append(weight)
        
    FinalList = []
    WFinal = []
    newel = True    
    while newel:
        newel = False
        NewList = []
        NewWeight = []
        for iel in range(len(ElList)):
            BaseEl = ElList[iel]
            ptc = BaseEl.momID.pop()
            weight = WList[iel]
            
            if len(BRdic[ptc]) == 0:       # Stable final state (LSP)
                BaseEl.momID = BaseEl.momID[0]
                FinalList.append(BaseEl)
                WFinal.append(weight)
                continue
            
            for BR in BRdic[ptc]:
                NewEl = copy.deepcopy(BaseEl)
                vertparts = []
                mass = []
                for x in BR.ids:
                    if x in SMSmethods.Reven:
                        vertparts.append(SMSmethods.Reven[x])
                    elif x in SMSmethods.Rodd:
                        mass.append(Massdic[x])
                        NewEl.momID.append(x)
                    else:
                        print 'SLHAdecomp: unknown particle:',x
                        return False
                    
                NewEl.particles.append(vertparts)
                if len(mass) == 1:
                    NewEl.masses.append(mass[0])
                else:
                    print 'SLHAdecomp: unknown decay (R-parity violation?)'
                    return False
                if weight*BR.br > sigcut:
                    NewList.append(NewEl)
                    NewWeight.append(weight*BR.br)
                    newel = True
                                
        if newel:
            ElList = copy.deepcopy(NewList)
            WList = copy.deepcopy(NewWeight)

#Combine 1branch elements according to production cross-section:
    SMSTopList = []    
    for ptcs in Xsec[Xsec.keys()[0]].keys():
        for iel in range(len(FinalList)):
            for jel in range(len(FinalList)):
                if ptcs[0] == ptcs[1] and jel < iel: continue       #Avoid double counting

                if FinalList[iel].momID == ptcs[0] and FinalList[jel].momID == ptcs[1]:
                    Els = SMSmethods.EElement()        
                    Els.B = [copy.deepcopy(FinalList[iel]),copy.deepcopy(FinalList[jel])]
                    weight = {}
                    for w in Xsec.keys():
                        weight.update({w : Xsec[w][ptcs]*WFinal[iel]*WFinal[jel]/(Pdic[ptcs[0]]*Pdic[ptcs[1]])})
                    Els.weight = weight                    
                    
                    if max(weight.values()) < sigcut: continue
                    
                    Einfo = Els.getEinfo()
                    Top = SMSmethods.GTop()
                    Top.vertnumb = Einfo["vertnumb"]
                    Top.vertparts = Einfo["vertparts"]
                    Top.ElList.append(Els)
                    for ib in range(len(Top.vertnumb)):
                        if len(Top.vertparts[ib]) != Top.vertnumb[ib]:
                            print 'SLHAdecomp: error creating topology'
                            return False
                    
#Add topology to topology list:  
                    SMSmethods.AddToList(Top,SMSTopList)
                    
           
    return SMSTopList        
