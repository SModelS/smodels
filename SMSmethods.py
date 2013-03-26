Rodd={ 
1000021 : "gluino", 1000022: "N1", 1000023 : "N2", 1000025 : "N3", 1000035 : "N4", 1000024 : "C1", 1000037 : "C2", 1000039 : "gravitino",  1000001 : "squark", 1000002 : "squark", 1000003 : "squark", 1000004 : "squark", 2000001 : "squark", 2000002 : "squark", 2000003 : "squark", 2000004 : "squark", 1000005 : "sbottom", 2000005 : "sbottom", 1000006 : "stop", 2000006 : "stop", 1000011 : "slepton", 1000013 : "slepton", 1000015 : "stau", 2000011 : "slepton", 2000013 : "slepton", 2000015 : "stau", 1000012 : "sneutrino", 1000014 : "sneutrino", 1000016 : "sneutrino", 2000012 : "sneutrino", 2000014 : "sneutrino", 2000016 : "sneutrino", -1000021 : "gluino", -1000022: "N1", -1000023 : "N2", -1000025 : "N3", -1000035 : "N4", -1000024 : "C1", -1000037 : "C2", -1000039 : "gravitino",  1000001 : "squark", -1000002 : "squark", -1000003 : "squark", -1000004 : "squark", -2000001 : "squark", -2000002 : "squark", -2000003 : "squark", -2000004 : "squark", -1000005 : "sbottom", -2000005 : "sbottom", -1000006 : "stop", -2000006 : "stop", -1000011 : "slepton", -1000013 : "slepton", -1000015 : "stau", -2000011 : "slepton", -2000013 : "slepton", -2000015 : "stau", -1000012 : "sneutrino", -1000014 : "sneutrino", -1000016 : "sneutrino", -2000012 : "sneutrino", -2000014 : "sneutrino", -2000016 : "sneutrino"  

}


Reven={
 25 : "higgs", 35 : "H0", 36 : "A0", 37 : "H+", -37 : "H-", 23 : "Z", 22 : "photon", 24 : "W+", -24 : "W-", 16 : "nu", -16 : "nu", 15 : "ta-", -15 : "ta+", 14 : "nu", -14 : "nu", 13 : "l-", -13 : "l+", 12 : "nu", -12 : "nu", 11 : "l-", -11 : "l+", 5 : "b", -5 : "b", 6 : "t+", -6 : "t-", 1 : "jet", 2 : "jet", 3 : "jet", 4 : "jet", 21 : "jet", -1 : "jet", -2 : "jet", -3 : "jet", -4 : "jet", -21 : "jet"
 }

PtcDic={
"l" : ["l+","l-"], "ta" : ["ta+","ta-"], "W" : ["W+","W-"], "t" : ["t+","t-"], "L+" : ["l+","ta+"], "L-" : ["l-","ta-"], "L" : ["l+","ta+","l-","ta-"]
}



class DBranch:
    def __init__(self):
        self.vertnumb = 0
        self.vertins = []
        self.ElList = []
    
    
class TElement:
    def __init__(self):
        self.masses = []
        self.particles = []
        self.momID = 0
        
class TMassEntry:
    def __init__(self):
        self.B1masses = []
        self.B2masses = []
        self.weight = []
                

class GTop:
    def __init__(self):
        self.B = [DBranch(),DBranch()]
        self.WeightList = []
    
    def AddElement(self, masses, particles, momID, weight = 0.):
        
#Sanity checks:
        if len(masses) != 2 or len(particles) !=2 or len(momID) !=2:
            print "AddElement: wrong number of branches"
            return False
        for i in range(2):
            if len(masses[i]) != self.B[i].vertnumb:
                print "AddElement: wrong number of masses"
                return False
            if len(particles[i]) != sum(self.B[i].vertins):
                print "AddElement: wrong number of particles"
                return False
            
#Create elements for each branch and add them to ElList:            
        for ibranch in range(2):
            NewEl = TElement()
            NewEl.masses = masses[ibranch]
            NewEl.particles = particles[ibranch]
            NewEl.momID = momID[ibranch]
            self.B[ibranch].ElList.append(NewEl)
     
        self.WeightList.append(weight)
        return True
    
    def OrderBranches(self):
        Branch1 = DBranch()
        Branch2 = DBranch()
        Branch1 = self.B[0]
        Branch2 = self.B[1]
        if (self.B[0].vertnumb < self.B[1].vertnumb) or (self.B[0].vertnumb == self.B[1].vertnumb and sum(self.B[0].vertins) < sum(self.B[1].vertins)):
            self.B = [Branch2,Branch1]


class MParticle:
    def __init__(self):
        self.pdg = 0
        self.status = 0
        self.moms = []
        self.px = 0.
        self.py = 0.
        self.pz = 0.
        self.e = 0.
        self.mass = 0.
        
        
def getNextEvent(filename): 
    
    line = " "
    PartList = []

#Find next event
    while line.find("<event>") == -1: line = filename.readline()
        
#Read event info:
    line = filename.readline()

#Get particles info:            
    line = filename.readline()  
    while line.find("</event>") == -1: 
        if line.find("#")>-1:
          line=line[:line.find('#')]
        if len(line)==0: 
          line=filename.readline()
          continue
        particle = MParticle()
        linep = [float(x) for x in line.split()]
        particle.pdg = int(linep[0])
        particle.status = int(linep[1])
        particle.moms = [int(linep[2]),int(linep[3])]
        particle.px = linep[6]
        particle.py = linep[7]
        particle.pz = linep[8]
        particle.e = linep[9]
        particle.mass = linep[10]
        
        PartList.append(particle)
        line = filename.readline()  


    return PartList


def getEventTop(PList, weight):
    import math, SMSglobals 
    
    leppT_cut = SMSglobals.min_lep_pt 
    jetpT_cut = SMSglobals.min_jet_pt
    taupT_cut = SMSglobals.min_tau_pt
    bpT_cut = SMSglobals.min_b_pt
    pT_cut = SMSglobals.min_pt

    ETop = GTop()
    
    momspdg = [0,0]
    mompos = [0,0]
    imom = 0
    
#First get Mothers:    
    for i in range(len(PList)):
        if PList[i].moms[0] == 1 or PList[i].moms[1] == 1:
            momspdg[imom] = PList[i].pdg
            mompos[imom] = i
            imom +=1
            
#Each mother = different branch. Loop over branchs: 
    for ib in range(2):
        mother = mompos[ib]
        newmom = mompos[ib]
        El = TElement()
        
        
        
        ndaugs = 2
        while ndaugs > 0:
            ndaugs = 0  
            nmoms = 0
            nvertins = 0      
            for i in range(len(PList)):
                if PList[i].moms[0] != mother+1 and PList[i].moms[1] != mother+1: continue
                if abs(PList[i].pdg) in Rodd:
                    newmom = i
                    ndaugs += 1
                    nmoms +=1
                elif abs(PList[i].pdg) in Reven:
                    pname = ptype(PList[i].pdg)
                    pT = math.sqrt(PList[i].px**2 + PList[i].py**2)
                    if pT < pT_cut: continue
                    if EqualParticles(pname,"l") and pT < leppT_cut: continue
                    if EqualParticles(pname,"ta") and pT < taupT_cut: continue
                    if EqualParticles(pname,"jet") and pT < jetpT_cut: continue
                    if EqualParticles(pname,"b") and pT < bpT_cut: continue

                    El.particles.append(pname)
                    nvertins +=1
                    ndaugs += 1
                else:
                    print "getEventTop: Unknown particle!"
                    return False
                
            if ndaugs > 1:
                if nmoms != 1 or ndaugs != nvertins+1:
                    print "getEventTop: R-parity violated!"
                    return False

                El.masses.append(PList[mother].mass)
                mother = newmom
                ETop.B[ib].vertins.append(nvertins)
                ETop.B[ib].vertnumb += 1
            elif ndaugs == 1:         #Compressed case->glue mother and daughter
                SMSglobals.ncomp[ib] +=1
                if mother == mompos[ib]:   # If compression happens at 1st decay
                    mompos[ib] = newmom    # redifine first mother
                    momspdg[ib] = PList[newmom].pdg
                mother  = newmom
                                 
        ETop.B[ib].ElList.append(El)    
        ETop.B[ib].ElList[0].momID = momspdg[ib]


    ETop.WeightList.append(weight)
#Sort branches
    Branch1 = ETop.B[0]
    Branch2 = ETop.B[1]
    if (ETop.B[0].vertnumb < ETop.B[1].vertnumb) or (ETop.B[0].vertnumb == ETop.B[1].vertnumb and sum(ETop.B[0].vertins) < sum(ETop.B[1].vertins)):
        ETop.B = [Branch2,Branch1]
        
       
    return ETop    
        
        
def ptype(pdg):
  p=int(pdg)
  if p in Rodd: return Rodd[p]
  if p in Reven: 
    return Reven[p]
  else:
      return False
  
  
def EqualTops(Top1,Top2):
    if (Top1.B[0].vertnumb == Top2.B[0].vertnumb and Top1.B[1].vertnumb == Top2.B[1].vertnumb):
        if Top1.B[0].vertins != Top2.B[0].vertins: return False 
        if Top1.B[1].vertins != Top2.B[1].vertins: return False
        return True
    elif (Top1.B[0].vertnumb == Top2.B[1].vertnumb and Top1.B[1].vertnumb == Top2.B[0].vertnumb):        
        if Top1.B[0].vertins != Top2.B[1].vertins: return False 
        if Top1.B[1].vertins != Top2.B[0].vertins: return False
        return True
    else:    
        return False

       
def EqualEls(El1,El2):
    if len(El1.particles) != len(El2.particles): return False
    if len(El1.masses) != len(El2.masses): return False
    if El1.particles != El2.particles: return False
    if not EqualMasses(El1.masses,El2.masses): return False
    return True

def SimEls(El1,El2):
    if len(El1.particles) != len(El2.particles): return False
    if len(El1.masses) != len(El2.masses): return False
    if not EqualParticles(El1.particles,El2.particles): return False
    if not EqualMasses(El1.masses,El2.masses): return False
    return True


def EqualMasses(mass1,mass2):
    from SMSglobals import massequiv,maxgap
    mass10 = mass1
    mass20 = mass2
    mass1v = []
    mass2v = []
    if mass1 == mass2: return True
    if type(mass1) == type(list()):
        if type(mass1[0]) == type(list()):            
                mass1 = mass1[0] + mass1[1]
                mass2 = mass2[0] + mass2[1]
        if type(mass1[0]) == type(list()) or type(mass2[0]) == type(list()):
            print "EqualMasses: wrong input format"
            return False
        
        mass1v = mass1
        mass2v = mass2
    else:
        mass1v.append(mass1)
        mass2v.append(mass2) 
        
    if len(mass1v) != len(mass2v):
        print "EqualMasses: wrong mass comparison",mass1v,mass2v,mass10,mass20
        return False
 
    for i in range(len(mass1v)):
        massA = mass1v[i]
        massB = mass2v[i]   
        if massA != massB:
            if 2.*abs(massA-massB)/(massA+massB) > massequiv: return False
            if abs(massA-massB) > maxgap: return False

    return True

#Compares 2 particle names or 2 name arrays. Allows for dictionary labels
#(Ex: L = l, l+ = l, l = l-,...)
def EqualParticles(ptype1,ptype2):
       
    ptype1v = []
    ptype2v = []
    if type(ptype1) == type(list()):
        if type(ptype1[0]) == type(list()):            
            ptype1 = ptype1[0] + ptype1[1]
            ptype2 = ptype2[0] + ptype2[1]
        if type(ptype1[0]) == type(list()) or type(ptype2[0]) == type(list()):
            print "EqualParticles: wrong input format"
            return False
        ptype1v = ptype1
        ptype2v = ptype2
    else:
        ptype1v.append(ptype1)
        ptype2v.append(ptype2) 

    if len(ptype1v) != len(ptype2v):
        return False
        
    for i in range(len(ptype1v)):
        if PtcDic.has_key(ptype1v[i]):
            ptypeA = PtcDic[ptype1v[i]]        
        else:
            ptypeA = [ptype1v[i]]
        if PtcDic.has_key(ptype2v[i]):
            ptypeB = PtcDic[ptype2v[i]]        
        else:
            ptypeB = [ptype2v[i]]

        if set(ptypeA) & set(ptypeB) == set([]): return False

    return True



def AddToList(SMSTop,SMSTopList):
    
    equaltops = -1
    equalels = -1
    for i in range(len(SMSTopList)):
        if EqualTops(SMSTop,SMSTopList[i]):
            equaltops = i
            for j in range(len(SMSTopList[i].B[0].ElList)):
                MEl_B1 = SMSTopList[i].B[0].ElList[j]
                MEl_B2 = SMSTopList[i].B[1].ElList[j]
                El_B1 = SMSTop.B[0].ElList[0]
                El_B2 = SMSTop.B[1].ElList[0]
                if (EqualEls(MEl_B1,El_B1) and EqualEls(MEl_B2,El_B2)):
                    equalels = j
                    ebr = [0,1]
                elif (EqualEls(MEl_B1,El_B2) and EqualEls(MEl_B2,El_B1)):
                    equalels = j
                    ebr = [1,0]
#If element exists, add weight:
                if equalels == j:
                    if len(SMSTopList[i].WeightList[j]) != len(SMSTop.WeightList[0]):
                        print "Wrong number of weights"
                    else:    
                        for iw in range(len(SMSTop.WeightList[0])):
                            SMSTopList[i].WeightList[j][iw] += SMSTop.WeightList[0][iw]
#For the "mothers" keep the smallest set of PDG mother IDs (not used in the analysis, only relevant to set a standard):
                    if min(abs(El_B1.momID),abs(El_B2.momID)) < min(abs(MEl_B1.momID),abs(MEl_B2.momID)):
                        for im in range(2): SMSTopList[i].B[im].ElList[j].momID = SMSTop.B[ebr[im]].ElList[0].momID
#For the masses, keep the highest sum:
                    MsumA = sum(MEl_B1.masses) + sum(MEl_B2.masses)
                    MsumB = sum(El_B1.masses) + sum(El_B2.masses)
                    if MsumB >= MsumA:
                        for im in range(2): SMSTopList[i].B[im].ElList[j].masses = [x for x in SMSTop.B[ebr[im]].ElList[0].masses]
                        
                        
#If topology and/or element does not exist, add:
    if equaltops == -1:
        SMSTopList.append(SMSTop)
    elif equalels == -1:
        allmass = [SMSTop.B[0].ElList[0].masses,SMSTop.B[1].ElList[0].masses]
        allparts = [SMSTop.B[0].ElList[0].particles,SMSTop.B[1].ElList[0].particles]
        allmom = [SMSTop.B[0].ElList[0].momID,SMSTop.B[1].ElList[0].momID]
        if not SMSTopList[equaltops].AddElement(allmass, allparts, allmom, SMSTop.WeightList[0]):
            print "Error adding element"
            print '\n'
            
            
def AddToAnalyses(SMSTopList):
    from SMSglobals import AnalysesRes    
    
#Find matching topologies:
    for iA in range(len(AnalysesRes)):
        for iE in range(len(SMSTopList)):  
            if not EqualTops(AnalysesRes[iA].Top,SMSTopList[iE]): continue
    
            SMSTop = SMSTopList[iE]
            
#Loop over analyses elements:            
            for i in range(len(AnalysesRes[iA].StrList)):
                Rptclist = eltostr(AnalysesRes[iA].StrList[i])
                MassList = []
                WeightList = []
#Collect all elements matching analyses element and order masses according
#to analyses branch order. If both orderings match, tag element with negative
#weights and add both possibilities to masslist
                for j in range(len(SMSTop.WeightList)):                
                    weight = SMSTop.WeightList[j]
                    El1 = SMSTop.B[0].ElList[j]
                    El2 = SMSTop.B[1].ElList[j] 
                    ptclist = [El1.particles,El2.particles]
                    masslist = [El1.masses,El2.masses]
                    pmatch = False        
                    
                    
                    if EqualParticles(ptclist,Rptclist):
                        pmatch = True
                        massnew = [masslist[0],masslist[1]]
                        weightnew = [x for x in weight]
                        if EqualParticles(ptclist,[Rptclist[1],Rptclist[0]]) and not EqualMasses(masslist[0],masslist[1]):
                            weightnew = [-x for x in weight]                    
                    elif EqualParticles(ptclist,[Rptclist[1],Rptclist[0]]):    
                        pmatch = True
                        massnew = [masslist[1],masslist[0]]
                        weightnew = [x for x in weight]
                    if pmatch:
                        MassList.append(massnew)
                        WeightList.append(weightnew)
                        if sum(weightnew) < 0.:
                            MassList.append([massnew[1],massnew[0]])
                            WeightList.append(weightnew)
                        
                    
                                
                  
#Now put all elements with ambiguous ordering at the bottom of the list 
#(to be added for last):
                for j in range(len(WeightList)):
                    if sum(WeightList[j]) < 0.:
                        Rmass = MassList.pop(j)
                        Rweight = WeightList.pop(j)
                        MassList.append(Rmass)
                        WeightList.append(Rweight)
                    

#Finally include all matching elements to AnalysesRes:
                for j in range(len(WeightList)):
                    match = -1
                    for k in range(len(AnalysesRes[iA].ElList[i].MassList)):
                        if (EqualMasses(MassList[j],AnalysesRes[iA].ElList[i].MassList[k])):
                            match = k
                            break                        
                        
                    if match >= 0:
                        for iw in range(len(AnalysesRes[iA].ElList[i].WeightList[match])):
                            AnalysesRes[iA].ElList[i].WeightList[match][iw] += abs(WeightList[j][iw])
                        MsumB = sum(MassList[j][0]) + sum(MassList[j][1])
                        MsumA = sum(AnalysesRes[iA].ElList[i].MassList[match][0]) + sum(AnalysesRes[iA].ElList[i].MassList[match][1])
                        if  MsumB >= MsumA:
                            AnalysesRes[iA].ElList[i].MassList[match] = [x for x in MassList[j]]     
                                
                    else:    #Mass not found, add to list
                        mass = [x for x in MassList[j]]
                        weight = [abs(x) for x in WeightList[j]]
                        AnalysesRes[iA].ElList[i].MassList.append(mass)
                        AnalysesRes[iA].ElList[i].WeightList.append(weight)


#Evaluate theoretical predictions for the analysis constraint\condition:
def Aeval(Analysis,instr):
    from SMSglobals import AnalysesRes
    
#First find equivalent topology in AnalysesRes:
    itop = -1
    for i in range(len(AnalysesRes)):
        if EqualTops(Analysis.Top,AnalysesRes[i].Top):
            itop = i
            break
    if itop == -1: return False
    
    ARes = AnalysesRes[itop]    
    
    
#First remove all blanks:
    outstr = instr.replace(" ","")
#Get ordered list of elements:
    El = []
    iels = 0
    while "[[" in outstr:  #String has element list                
        st = outstr[outstr.find("[["):outstr.find("]]")+2] #Get duplet        
        ptclist = eltostr(st)     # Get particle list
#Syntax checks:
        for ib in range(2):
            for ip in range(len(ptclist[ib])):
                if not ptclist[ib][ip] in Reven.values() and not PtcDic.has_key(ptclist[ib][ip]):
                    print "Ceval: Unknown particle ",ptclist[ip]
                    return False
        outstr = outstr.replace(st,"El["+str(iels)+"]")  # Replace element
        El.append(ptclist)
        iels +=1
        
        
#Get list of expressions (separated by commas):
    outstrv = outstr.rsplit(",")
#Find elements in AnalysesRes corresponding to El[i]:
    ieq =[]
    for i in range(len(El)):
        for j in range(len(ARes.StrList)):
            Rptc = eltostr(ARes.StrList[j])
            if El[i] == Rptc:
                ieq.append(j)
                break
#Loop through masslist and evaluate with weights:
    Allres = []
    for i in range(len(ARes.ElList[0].MassList)):
        resw = []
        for iw in range(len(ARes.ElList[0].WeightList[0])):
            for j in range(len(El)):            
                El[j] = ARes.ElList[ieq[j]].WeightList[i][iw]
            eout = [Ceval(x,El) for x in outstrv]
            resw.append(eout)
        Allres.append([ARes.ElList[0].MassList[i],resw])
        
    return Allres      
                

#Get upper limit on sigma*BR for a specific array of masses from plot
def GetPlotLimit(masses,plot):
    sigmax = 0.
    return sigmax
    
    
    
    
    
    

    
def eltostr(invar):
    if type(invar) == type(list()):
        st = str(invar).replace("'","")
        st = st.replace(" ","")
        return st
    elif type(invar) == type(str()):
        st = invar.replace(" ","")
        st = st[st.find("[["):st.find("]]")+2]
        st_B1 = st[0:st.find("],[")]    
        st_B2 = st[st.find("],[")+3:]
        newstr = [st_B1,st_B2]
        ptclist = []   
        for k in range(2):    
            newstr[k] = newstr[k].replace("]","")
            newstr[k] = newstr[k].replace("[","")
            newstr[k] = newstr[k].replace(" ","")
            ptclist.append(newstr[k].rsplit(","))
        return ptclist
    else:
        print "ElToStr: Wrong input"
        return False
        
            
            
class GTopAnalyses:
     
    def __init__(self):
        self.Top = GTop()
        self.ElList = []
        self.StrList = []

#Generates a common mass list for all the masses appearing in the topology
#Useful for evaluating constraints and conditions        
    def CreateAll(self):

#Collect all masses:
        AllMasses = []
        for i in range(len(self.ElList)):
            for j in range(len(self.ElList[i].MassList)):
                AllMasses.append(str(self.ElList[i].MassList[j]))
                wlength = len(self.ElList[i].WeightList[j])
#Remove duplicates:
        AllMasses = set(AllMasses)
        AllMasses = [eval(x) for x in AllMasses]

#Remove Equivalent masses:
        for i in range(len(AllMasses)-1):
            massA = AllMasses[i]
            if massA == []: continue
            for j in range(i+1,len(AllMasses)):
                massB = AllMasses[j]
                if massB == []: continue
                if EqualMasses(massA,massB):
                    MsumA = sum(massA[0]) + sum(massA[1])
                    MsumB = sum(massB[0]) + sum(massB[1])
                    if  MsumB >= MsumA:
                        AllMasses[i] = []
                    else:    
                        AllMasses[j] = []
                    
        AllMassesNew = []        
        for i in range(len(AllMasses)):
            if AllMasses[i] != []: AllMassesNew.append(AllMasses[i])
        AllMasses = [x for x in AllMassesNew]

        for i in range(len(self.ElList)):
            NewWeightList = []
            for j in range(len(AllMasses)):
                weight = [0.]*wlength
                for imass in range(len(self.ElList[i].MassList)):
                    if EqualMasses(AllMasses[j],self.ElList[i].MassList[imass]):
                        for iw in range(wlength):
                            weight[iw] += self.ElList[i].WeightList[imass][iw]
                NewWeightList.append(weight)
            

            self.ElList[i].MassList = [x for x in AllMasses]   
            self.ElList[i].WeightList = [x for x in NewWeightList]
                        
            
        
       
        
class TElementAnalyses:                  
    def __init__(self):
        self.MassList = []
        self.WeightList = []
            

class EAnalysis:    
    def __init__(self):
        self.label = ""
        self.Top = GTop()
        self.sqrts = 8.
        self.lum = 0.
        self.constraints = {}
        self.plots = {}
        self.plotpath = ""


#Given the constraints dictionary, automatically fill AnalysesRes with
# the element list for all the elemenents appearing in the constraints 
# and conditions, skipping repeated or equivalent ones.
    def FillElements(self):
        from SMSglobals import AnalysesRes
       
        ListOfStrs = []
        vertnumb = [self.Top.B[0].vertnumb,self.Top.B[1].vertnumb]
        vertins = [self.Top.B[0].vertins,self.Top.B[1].vertins]
# Syntax check:        
        for k in range(2):
            if len(vertins[k]) != vertnumb[k]:
                print "FillElements: Wrong syntax"
                return False
        
        inelements = self.constraints.items()

        for iii in range(len(inelements)):
            for ii in range(2):    

                con = inelements[iii][ii].replace(" ","")
                while "[" in con:  #String has element list                
                    st = con[con.find("[["):con.find("]]")+2] #Get duplet
                    con = con.replace(st,"")  # Remove element duplet
                    ptclist = eltostr(st)     # Get particle list
#Syntax checks:
                    for k in range(2):
                        if len(ptclist[k]) != sum(vertins[k]):
                            print "FillElements: Wrong syntax"
                            return False
                        for ip in range(len(ptclist[k])):
                            if not ptclist[k][ip] in Reven.values() and not PtcDic.has_key(ptclist[k][ip]):
                                print "FillElements: Unknown particle ",ptclist[k][ip]
                                return False        
                    ListOfStrs.append(st)
                    
#Find correct AnalysesRes topology to insert the strings:
        itop = -1
        for i in range(len(AnalysesRes)):
            ATop = AnalysesRes[i].Top
            if EqualTops(ATop,self.Top):
                itop = i
                break
            
        if itop == -1:
            NewAnalyses = GTopAnalyses()
            for i in range(2):
               NewAnalyses.Top.B[i].vertnumb = vertnumb[i]
               NewAnalyses.Top.B[i].vertins = vertins[i]
            AnalysesRes.append(NewAnalyses)
            itop = len(AnalysesRes)-1
            
            
#Now find if element already exists:
        ListOfStrs = set(ListOfStrs)
        while len(ListOfStrs) > 0:
            strA = ListOfStrs.pop()        #Check for same branch order
#            strB = eltostr(strA)           
#            strB = [strB[1],strB[0]]
#            strB = eltostr(strB)
            Allstr = AnalysesRes[itop].StrList
#            if (strA in Allstr) or (strB in Allstr): continue #Skip element
            if (strA in Allstr): continue #Skip element

            NewEl = TElementAnalyses()
            AnalysesRes[itop].ElList.append(NewEl)
            AnalysesRes[itop].StrList.append(strA)
            
            
        
        
#Defines similar function when comparing two list of numbers.
#If any of the elements differ by < 10%, returns True
def similar(els):
    for i in range(len(els)):
        for j in range(i+1,len(els)):
            if els[i] != els[j]: 
                if 2.*abs(els[i]-els[j])/abs(els[i]+els[j]) > 0.1: return False  
    return True          
    
#Routine to evaluate the analyses conditions and constraints.
#Flexible version of eval to allow for additional operators, 
#such as ~ (= similar)
def Ceval(instring,El):
    
    run = instring.replace(" ","")  #Remove blanks
    if "~" in run:
        simels = run.split("~")
        run = 'similar(' + str(simels) + ')'
        run = run.replace("'","")
    return eval(run)

  
    






        
        


def runpythia(slhafile,nevts = 10000, Sqrts = 7):
    import sys
    
    pyslhadir = "/home/lessa/SMS_pythia/pyslha-1.4.3"
    sys.path.append(pyslhadir)
    import PMSSM


    workdir = "/home/lessa/SMS_pythia"
    PMSSM.installdir="%s/Pythia_LHE" % workdir
    PMSSM.etcdir="%s/Pythia_LHE/etc" % workdir
    PMSSM.logdir="%s/Pythia_LHE/log" % workdir

    PMSSM.verbose=True
    PMSSM.basedir = "%s/Pythia_LHE" % workdir
    PMSSM.datadir = "%s/Pythia_LHE/data" % workdir

#run pythia
    print "running pythia"
    D = PMSSM.runPythiaLHE ( nevts, slhafile, " ", sqrts=Sqrts)
    print "done running"
    return D
