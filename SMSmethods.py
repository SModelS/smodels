#---------------Dictionaries:
Rodd={ 
1000021 : "gluino", 1000022: "N1", 1000023 : "N2", 1000025 : "N3", 1000035 : "N4", 1000024 : "C1", 1000037 : "C2", 1000039 : "gravitino",  1000001 : "squark", 1000002 : "squark", 1000003 : "squark", 1000004 : "squark", 2000001 : "squark", 2000002 : "squark", 2000003 : "squark", 2000004 : "squark", 1000005 : "sbottom", 2000005 : "sbottom", 1000006 : "stop", 2000006 : "stop", 1000011 : "slepton", 1000013 : "slepton", 1000015 : "stau", 2000011 : "slepton", 2000013 : "slepton", 2000015 : "stau", 1000012 : "sneutrino", 1000014 : "sneutrino", 1000016 : "sneutrino", 2000012 : "sneutrino", 2000014 : "sneutrino", 2000016 : "sneutrino", -1000021 : "gluino", -1000022: "N1", -1000023 : "N2", -1000025 : "N3", -1000035 : "N4", -1000024 : "C1", -1000037 : "C2", -1000039 : "gravitino",  1000001 : "squark", -1000002 : "squark", -1000003 : "squark", -1000004 : "squark", -2000001 : "squark", -2000002 : "squark", -2000003 : "squark", -2000004 : "squark", -1000005 : "sbottom", -2000005 : "sbottom", -1000006 : "stop", -2000006 : "stop", -1000011 : "slepton", -1000013 : "slepton", -1000015 : "stau", -2000011 : "slepton", -2000013 : "slepton", -2000015 : "stau", -1000012 : "sneutrino", -1000014 : "sneutrino", -1000016 : "sneutrino", -2000012 : "sneutrino", -2000014 : "sneutrino", -2000016 : "sneutrino"  

}


Reven={
 25 : "higgs", 35 : "H0", 36 : "A0", 37 : "H+", -37 : "H-", 23 : "Z", 22 : "photon", 24 : "W+", -24 : "W-", 16 : "nu", -16 : "nu", 15 : "ta-", -15 : "ta+", 14 : "nu", -14 : "nu", 13 : "l-", -13 : "l+", 12 : "nu", -12 : "nu", 11 : "l-", -11 : "l+", 5 : "b", -5 : "b", 6 : "t+", -6 : "t-", 1 : "jet", 2 : "jet", 3 : "jet", 4 : "jet", 21 : "jet", -1 : "jet", -2 : "jet", -3 : "jet", -4 : "jet", -21 : "jet"
 }

PtcDic={
"l" : ["l+","l-"], "ta" : ["ta+","ta-"], "W" : ["W+","W-"], "t" : ["t+","t-"], "L+" : ["l+","ta+"], "L-" : ["l-","ta-"], "L" : ["l+","ta+","l-","ta-"]
}


#--------------------Classes:
    
    
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
      
class TElement:
    def __init__(self):
        self.masses = []
        self.particles = []
        self.momID = 0
        
class DBranch:
    def __init__(self):
        self.vertnumb = 0
        self.vertparts = []
        self.ElList = []

class GTop:
    def __init__(self):
        self.B = [DBranch(),DBranch()]
        self.WeightList = []

#Adds element to ElLists in branches
#OBS: input must be given with the correct branch ordering!  
    def AddElement(self, NewElement, weight = {}):
        
#Sanity checks:
        for i in range(2):
            masses = NewElement[i].masses
            particles = NewElement[i].particles
            if len(masses) != self.B[i].vertnumb:
                print "AddElement: wrong number of masses"
                return False
            if (self.B[i].vertparts[self.B[i].vertnumb-1] == 0 and len(particles) != self.B[i].vertnumb-1) or (self.B[i].vertparts[self.B[i].vertnumb-1] != 0 and len(particles) != self.B[i].vertnumb): 
                print "AddElement: wrong number of particles"
                return False
            else:
                for ipt in range(len(particles)):
                    if len(particles[ipt]) != self.B[i].vertparts[ipt]:
                        print "AddElement: wrong number of particles"
                        return False

            
#Create elements for each branch and add them to ElList:            
        for ibranch in range(2):
            self.B[ibranch].ElList.append(NewElement[ibranch])
     
        neweight = {}
        neweight.update(weight)
        self.WeightList.append(neweight)
        return True
    
class EAnalysis:    
    def __init__(self):
        self.label = ""
        self.Top = GTop()
        self.sqrts = 8.
        self.lum = 0.
        self.results = {}
        self.plots = {}
        self.run = ""
        self.masscomp = 0.001


#Given the constraints dictionary, automatically fill the element list with the
#elements corresponding to the strings in the dictionary, skipping repeated ones
    def GenerateElements(self):
        import sys
       
        ListOfStrs = []
        vertnumb = [self.Top.B[0].vertnumb,self.Top.B[1].vertnumb]
        vertparts = [self.Top.B[0].vertparts,self.Top.B[1].vertparts]
#Syntax check:        
        for k in range(2):
            if len(vertparts[k]) != vertnumb[k]:
                print "GenerateElements: Wrong syntax"
                return False
        
#Get all element strings:        
        inelements = self.results.items()
        
        for iii in range(len(inelements)):
            for ii in range(2):    

                con = inelements[iii][ii].replace(" ","")
                while "[" in con:  #String has element                
                    st = con[con.find("[[["):con.find("]]]")+3] #Get duplet
                    con = con.replace(st,"")  # Remove element duplet
                    ptclist = eltostr(st)     # Get particle list
#Syntax checks:
                    for ib in range(2):
                        for ipt in range(len(ptclist[ib])):
                            if len(ptclist[ib][ipt]) != vertparts[ib][ipt]:
                                print "GenerateElements: Wrong syntax2"
                                return False
                            for ptc in ptclist[ib][ipt]:
                                if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                                    print "GenerateElements: Unknown particle ",ptc
                                    return False        
                    ListOfStrs.append(st)
                    
                    
#Remove repeated elements:
        ListOfStrs = set(ListOfStrs)
#Now add all elements to element list        
        while len(ListOfStrs) > 0:
            ptclist = eltostr(ListOfStrs.pop())            
            for ib in range(2):
                NewEl = TElement()
                NewEl.particles = ptclist[ib]
                self.Top.B[ib].ElList.append(NewEl)
                self.Top.WeightList.append([])
    
    

#Check if the plots listed in results exist and 
#store bin information in masscomp for combining masses
    def GetPlots(self):
        import SMSResults
        
        self.masscomp = {}
        run = self.run        
        if run == "": run = None  #If run has not been defined, use latest
        for res in self.results.keys():
            if not self.plots.has_key(res):
                print "GetPlots: Plot for result",res,"in Analysis",self.label,"not found"
                topo = ""
                ana = []
            else:
                topo = self.plots[res][0]
                analyses = self.plots[res][1]
                
            binx = 0.
            biny = 0.
            xmin = 0.
            ymin = 0.
            for ana in analyses:
                if not SMSResults.exists(ana,topo,run):
                    print "GetPlots: Histogram for ",topo," in ",ana," for run ",run," not found"
                    continue
#                if binx == 0.:
#                    binx = SMSResults.getXBinWidth(ana,topo,run)
#                    biny = SMSResults.getYBinWidth(ana,topo,run)
#                    xmin = SMSResults.getMinX(ana,topo,run)
#                    ymin = SMSResults.getMinY(ana,topo,run)
#                else:    
#                    binx = min(SMSResults.getXBinWidth(ana,topo,run),binx)
#                    biny = min(SMSResults.getYBinWidth(ana,topo,run),biny)
#                    xmin = min(SMSResults.getMinX(ana,topo,run),xmin)
#                    ymin = min(SMSResults.getMinY(ana,topo,run),ymin)
                    
            if binx == 0.:  #If no histogram has been found, use default values
                binx = 20.
                biny = 20.    
                    
            self.masscomp.update({res: [binx,biny,xmin,ymin]}) #Store bin information
            

        


#--------------------Methods:
    
#Event reader. Returns a list of MParticles
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

#Reads particle list and returns an array with the mother
#PDG codes:
def getMom(PList):

    momspdg = []
    imom = 0
    for i in range(len(PList)):
        if PList[i].moms[0] == 1 or PList[i].moms[1] == 1:
            momspdg.append(PList[i].pdg)
            imom += 1
    if momspdg[0] > momspdg[1]:
        momspdg[0], momspdg[1] = momspdg[1], momspdg[0]
            
    if imom != 2:
        print "getMom: Number of mother particles != 2"
        return False
    else:
        return momspdg

#Reads particle list (PList) from event and generates Topology (GTop)
#with only one element in ElList, corresponding to the event
#If the DoCompress and/or DoInvisible flags are on, also generate
#compressed topologies with there are small mass gaps and/or neutrinos
#emitted in the last step of the cascade ("effective LSP")
def getEventTop(PList, weight = {}):
    import SMSglobals 
        
    ETopList = []
    ETop = GTop()
    
    momspdg = [0,0]
    mompos = [0,0]
    imom = 0
    nptcs = 0      #Particle counter just for sanity checks
    
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
        
        nptcs += 1
        
        ndaugs = 2
        while ndaugs > 0:
            ndaugs = 0  
            nmoms = 0
            nvertparts = 0
            El.particles.append([])
            for i in range(len(PList)):
                if PList[i].moms[0] != mother+1 and PList[i].moms[1] != mother+1: continue
                if abs(PList[i].pdg) in Rodd:
                    newmom = i
                    ndaugs += 1
                    nmoms += 1
                    nptcs += 1
                elif abs(PList[i].pdg) in Reven:
                    pname = ptype(PList[i].pdg)
                    El.particles[ETop.B[ib].vertnumb].append(pname)
                    nvertparts +=1
                    ndaugs += 1
                    nptcs += 1
                else:                    
                    print "getEventTop: Unknown particle!"
                    return False

            El.masses.append(PList[mother].mass)
            mother = newmom
            ETop.B[ib].vertparts.append(nvertparts)
            ETop.B[ib].vertnumb += 1


        if El.particles[ETop.B[ib].vertnumb-1] == []:
            El.particles.pop()     #Remove last empty insertion if LSP is stable

        ETop.B[ib].ElList.append(El)    
        ETop.B[ib].ElList[0].momID = momspdg[ib]

    neweight = {}
    neweight.update(weight)
    ETop.WeightList.append(neweight)

#Check to see if all particles (minus initial state) have been included 
#(Only for later. Right now higgsses decays are not 
#included in the topology, so ignore warning for now)
#    if nptcs != len(PList)-2:
#        for ipt in PList:
#            print ipt.pdg,ipt.mass
#        print "nptcs = ",nptcs    
#        print "getEventTop: Error reading event!"
#        return False
#    else:
    ETopList.append(ETop)
        
        
#Check for compressed topologies:
    if SMSglobals.DoCompress:
        ETopComp = GTop()
        neweight = {}
        neweight.update(weight)
        ETopComp.WeightList.append(neweight)
        Comp = False
#Loop over branches        
        for ib in range(2):
            ETopComp.B[ib].ElList.append(TElement())            
            for i in range(len(ETop.B[ib].ElList[0].masses)):
                massA = ETop.B[ib].ElList[0].masses[i]
                if i < len(ETop.B[ib].ElList[0].masses)-1:
                    massB = ETop.B[ib].ElList[0].masses[i+1]
                    if abs(massA-massB) < SMSglobals.minmassgap:
                        Comp  = True
                        continue           #If compressed, skip this vertex
              
                ETopComp.B[ib].ElList[0].masses.append(massA) 
                vertparts = ETop.B[ib].vertparts
                ptcs = ETop.B[ib].ElList[0].particles[sum(vertparts[:i]):sum(vertparts[:i+1])]
                for ptc in ptcs:
                    ETopComp.B[ib].ElList[0].particles.append(ptc)
                ETopComp.B[ib].vertnumb += 1
                ETopComp.B[ib].vertparts.append(len(ptcs))               
        if Comp:
            SMSglobals.nComp += 1
            ETopList.append(ETopComp)

                
            
            
#Check for effective LSPs (LSP + neutrino = LSP'):                   
    if SMSglobals.DoInvisible:
        ETopIComp = GTop()
        neweight = {}
        neweight.update(weight)
        ETopIComp.WeightList.append(neweight)
        if SMSglobals.DoCompress:  #Do invisible compression+mass compression
            ETopIComp = ETopComp
        else:    
            ETopIComp = ETop        #or just invisble compression
        Icomp = False
#Loop over branches        
        for ib in range(2):
            nnu = 0
            ptcs = ETopIComp.B[ib].ElList[0].particles
            i = len(ptcs)-1
            nnu = 0
            while i >= 0 and ptcs[i] == 'nu':   #Count number of nus emitted at the end of the cascade decay
                nnu += 1
                i -= 1
            ilsp = ETopIComp.B[ib].vertnumb    #Position of effective LSP
            vertin = ETopIComp.B[ib].vertparts[ilsp-1]
            while nnu >= vertin and ilsp >= 1:                            
                ilsp -= 1
                vertin += ETopIComp.B[ib].vertparts[ilsp-1]
                
   
            if ilsp < ETopIComp.B[ib].vertnumb-1:
                Icomp = True    
                ETopIComp.B[ib].vertnumb = ilsp + 1
                ETopIComp.B[ib].vertparts = ETopIComp.B[ib].vertparts[:ilsp]
                ETopIComp.B[ib].vertparts.append(0)
                ETopIComp.B[ib].ElList[0].particles = ptcs[:sum(ETopIComp.B[ib].vertparts)]
                ETopIComp.B[ib].ElList[0].masses = ETopIComp.B[ib].ElList[0].masses[:ilsp+1]
                
        if Icomp:
            SMSglobals.nInvis += 1
            ETopList.append(ETopIComp)
            
    return ETopList           
                
                
     
                        

#Converts pdg number to particle name according to the dictionaries Rodd
# and Reven        
def ptype(pdg):
  p=int(pdg)
  if p in Rodd: return Rodd[p]
  if p in Reven: 
    return Reven[p]
  else:
      return False
  
  
#Compares two global topologies. Returns true if they have the same
#number of vertices and particles, independent of branch ordering.
def EqualTops(Top1,Top2):
    if (EqualBranches(Top1.B[0],Top2.B[0]) and EqualBranches(Top1.B[1],Top2.B[1])) or (EqualBranches(Top1.B[0],Top2.B[1]) and EqualBranches(Top1.B[1],Top2.B[0])): return True
    else:
        return False
    
#Compares two branches. Returns true if they have the same
#number of vertices and particles
def EqualBranches(B1,B2):
    if B1.vertparts != B2.vertparts: return False 
    if B1.vertnumb != B2.vertnumb: return False
    return True
    

#Compare two elements or array of elements.
#If all masses and particles are equal, returns True,
#otherwise returns False       
def EqualEls(El1,El2):
    
    if len(El1) != len(El2): return False
    
#If it is an array of elements:    
    if type(El1) == type(list()) and type(El1[0]) == type(TElement()):
        if type(El2) != type(list()): return False
        if type(El2[0]) != type(TElement()): return False                
        for i in range(len(El1)):
            if El1[i].particles != El2[i].particles: return False
            if El1[i].masses != El2[i].masses: return False
            
#If it is a single element:            
    else:
        if El1.particles != El2.particles: return False
        if El1.masses != El2.masses: return False

    return True

#Compare two elements or array of elements.
#If all masses and particles are similar, returns True,
#otherwise returns False       
def SimEls(El1,El2):
    
    if len(El1) != len(El2): return False
    
#If it is an array of elements:    
    if type(El1) == type(list()) and type(El1[0]) == type(TElement()):
        if type(El2) != type(list()): return False
        if type(El2[0]) != type(TElement()): return False                
        for i in range(len(El1)):
            if not SimParticles(El1[i].particles,El2[i].particles): return False
            if El1[i].masses != El2[i].masses: return False
            
#If it is a single element:            
    else:
        if not SimParticles(El1.particles,El2.particles): return False
        if El1.masses != El2.masses: return False

    return True


#Compares 2 particle names or 2 nested name arrays. Allows for dictionary labels
#(Ex: L = l, l+ = l, l = l-,...)
def SimParticles(ptype1,ptype2):

    ptype1v = [ptype1]
    ptype2v = [ptype2]
 
#First flatten nested arrays: 
    isList = True
    while isList:
        newptype1v = []
        newptype2v = []
        if len(ptype1v) != len(ptype2v): return False        
        for i in range(len(ptype1v)):
            if type(ptype1v[i]) == type(list()):
                if len(ptype1v[i]) != len(ptype2v[i]): return False
                for j in range(len(ptype1v[i])):
                    newptype1v.append(ptype1v[i][j])
                    newptype2v.append(ptype2v[i][j])
            else:
                newptype1v.append(ptype1v[i])
                newptype2v.append(ptype2v[i])
                
        ptype1v = newptype1v
        ptype2v = newptype2v
        isList = False
        for i in range(len(ptype1v)):
            if type(ptype1v[i]) == type(list()): isList = True
            if type(ptype2v[i]) == type(list()): isList = True

        
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



#Check if elements in SMSTop matches an entry in SMSTopList. If it does, add weight.
#If the same topology exists, but not the same element, add element.
#If neither element nor topology exist, add the new topology and all its elements
def AddToList(SMSTop,SMSTopList):
    
    for inew in range(len(SMSTop.B[0].ElList)):
        NewEl_a = [SMSTop.B[0].ElList[inew],SMSTop.B[1].ElList[inew]]  #Check both orderings
        NewEl_b = [SMSTop.B[1].ElList[inew],SMSTop.B[0].ElList[inew]]                
        equaltops = -1
        equalels = -1
        i = -1
        while (equaltops < 0 or equalels < 0) and i < len(SMSTopList)-1:
            i += 1
            if EqualTops(SMSTop,SMSTopList[i]):    #First look for matching topology
                equaltops = i
            else: continue 
                
            for j in range(len(SMSTopList[i].B[0].ElList)):  #Search for matching element
                OldEl = [SMSTopList[i].B[0].ElList[j],SMSTopList[i].B[1].ElList[j]]
                if EqualEls(OldEl,NewEl_a):
                    equalels = j
                    NewEl = NewEl_a
                    break                   
                elif EqualEls(OldEl,NewEl_b):
                    equalels = j
                    NewEl = NewEl_b
                    break
                    
                    
#If element exists, add weight:
        if equalels >= 0:
            if len(SMSTopList[equaltops].WeightList[equalels]) != len(SMSTop.WeightList[inew]):
                print "Wrong number of weights"
            else:
                for w in SMSTop.WeightList[inew].keys():
                    SMSTopList[equaltops].WeightList[equalels][w] += SMSTop.WeightList[inew][w]
                    
    
                    
#When combining elements, keep the smallest set of PDG mother IDs (not used in the analysis, only relevant to set a standard):
                if min(abs(NewEl[0].momID),abs(NewEl[0].momID)) < min(abs(OldEl[0].momID),abs(OldEl[1].momID)):
                    for ib in range(2): SMSTopList[equaltops].B[ib].ElList[equalels].momID = NewEl[ib].momID
                        
                        
#If topology and/or element does not exist, add:
    if equaltops == -1:
        SMSTopList.append(SMSTop)
    elif equalels == -1:
        if EqualBranches(SMSTop.B[0],SMSTopList[equaltops].B[0]):
            NewEl = [SMSTop.B[0].ElList[inew],SMSTop.B[1].ElList[inew]]
        else:
            NewEl = [SMSTop.B[1].ElList[inew],SMSTop.B[0].ElList[inew]]
        if not SMSTopList[equaltops].AddElement(NewEl,SMSTop.WeightList[inew]):
            print "Error adding element"
            print '\n'
            

#Loop over all elements in SMSTopList and add the weight to the 
#matching elements in Analysis.
def AddToAnalysis(SMSTopList,Analysis):
    
    for itop in range(len(SMSTopList)):
        NewTop = SMSTopList[itop]    
#Check if topologies match:
        if not EqualTops(NewTop,Analysis.Top): continue
        
#Loop over (event) element list:
        for iel in range(len(NewTop.B[0].ElList)):
#Loop over analysis elements:
            for jel in range(len(Analysis.Top.B[0].ElList)):
                NewEl_a = [NewTop.B[0].ElList[iel],NewTop.B[1].ElList[iel]]
                NewEl_b = [NewTop.B[1].ElList[iel],NewTop.B[0].ElList[iel]]
                oldptcs = [Analysis.Top.B[0].ElList[jel].particles,Analysis.Top.B[1].ElList[jel].particles]
                newptcs_a = [NewEl_a[0].particles,NewEl_a[1].particles]
                newptcs_b = [NewEl_b[0].particles,NewEl_b[1].particles]
                
                
                
#Check if particles match (independent of branch order)
                ptcmatch = SimParticles(newptcs_a,oldptcs) or SimParticles(newptcs_b,oldptcs)
                if not ptcmatch: continue

                weight = {}
                weight.update(NewTop.WeightList[iel])

                
#Loop over nested mass list in Analysis element:                
                massmatch = False
                for imass in range(len(Analysis.Top.B[0].ElList[jel].masses)):
                    OldEl = [TElement(),TElement()]
                    for ib in range(2):
                        OldEl[ib].particles = oldptcs[ib]
                        OldEl[ib].masses = Analysis.Top.B[ib].ElList[jel].masses[imass]
                        
                        
#Check if elements match (with identical masses) and  add weight
                    if SimEls(NewEl_a,OldEl) or SimEls(NewEl_b,OldEl):
                        massmatch = True
                        for w in weight.keys():
                            Analysis.Top.WeightList[jel][imass][w] += weight[w]

#If particles match, but not masses, add element to nested 
# mass list and nested weight list
#(check for both branch orderings. If both match, add both mass orderings):
                if not massmatch:
                    if SimParticles(newptcs_a,oldptcs):
                        for ib in range(2):
                            Analysis.Top.B[ib].ElList[jel].masses.append(NewEl_a[ib].masses)                            
                        Analysis.Top.WeightList[jel].append(weight)
                        massmatch = True
                    if SimParticles(newptcs_b,oldptcs) and (not massmatch or NewEl_b[0].masses != NewEl_b[1].masses):
                        for ib in range(2):
                            Analysis.Top.B[ib].ElList[jel].masses.append(NewEl_b[ib].masses)
                        Analysis.Top.WeightList[jel].append(weight)
           


#Replace mass array by its binned equivalent using the information
#in binpars:
def BinMass(mass,binpars=[20.,20.,0.,0.]):
    
    newmass = [[],[]]
    for ib in range(2):
        for imass in range(len(mass[ib])):
            if imass < len(mass[ib])-1:  #Parameters for mx (mother and intermediate masses)
                binw = binpars[0]
                mmin = binpars[2]
            else:                        #Parameters for LSP
                binw = binpars[1]
                mmin = binpars[3]
            Nbins = int((mass[ib][imass]-mmin)/binw)
            newmass[ib].append(Nbins*binw + mmin + binw/2.)
            
    return newmass       
        
                          
                        
#For a list of equivalent masses, compute an average mass (or mass array)
#using the defined method. 
#harmonic = harmonic mean
#mean = standard mean
def MassAvg(equivin, method = "mean"):
    import numpy

    N = len(equivin)
    if N == 0:
        print "MassAvg: Empty array"
        return False
    if N == 1: return equivin[0]

    if type(equivin[0]) != type(list()):
        equivinBr = [equivin]
#In case the input has 2 branches of different sizes, average
#each one individually
    elif len(equivin[0]) == 2 and type(equivin[0][0]) == type(list()):
        if len(equivin[0][0]) != len(equivin[0][1]):
            equivinBr = [[],[]]
            for mass in equivin:
                equivinBr[0].append(mass[0])
                equivinBr[1].append(mass[1])
        else:
            equivinBr = [equivin]

    massout = []    
    for ib in range(len(equivinBr)):            
        equivmasses = numpy.array(equivinBr[ib])  #Generate numpy array
    
#Sanity checks:    
        for mass in equivmasses.flat:
            if mass == 0.:
                print "MassAvg: Zero mass!"
                return False
            if mass < 0.:
                print "MassAvg: Negative mass!"
                return False

        if method == "mean":
            massavg = equivmasses[0]
        elif method == "harmonic":
            massavg = 1./equivmasses[0]
        else:
            print "MassAvg: Unknown method"
            return False

        for imass in range(1,N):
            mass = equivmasses[imass]
            if mass.shape != massavg.shape:        #Sanity check
                print "MassAvg: Wrong input"
                return False
            if method == "mean":
                massavg = massavg + mass
            elif method == "harmonic":
                massavg = massavg + 1./mass
            
        if method == "mean":
            massavg = massavg/float(N)
        elif method == "harmonic":
            massavg = float(N)/massavg
    
        if massavg.shape != equivmasses[0].shape:
            print "MassAvg: Error computing average"
            return False
        
        massout.append(massavg.tolist())
    
    if len(massout) == 1:    
        return massout[0]
    else:
        return massout
    
                    

#Evaluate theoretical predictions for the analysis result and conditions:
def EvalRes(res,Analysis):
    import copy
    
    if not Analysis.plots.has_key(res) or not Analysis.results.has_key(res) or not Analysis.masscomp.has_key(res):
        print "EvalRes: Wrong analysis input"
        return False
    
    binpars = Analysis.masscomp[res] #Parameters for mass binning
    
    
#Create new element list with all masses replaced by
#their binned values
    NewTop = copy.deepcopy(Analysis.Top)
    for iel in range(len(NewTop.B[0].ElList)):      
        masslist = []
        weightlist = []      
        while len(NewTop.B[0].ElList[iel].masses) > 0:
            mass = [NewTop.B[0].ElList[iel].masses.pop(0),NewTop.B[1].ElList[iel].masses.pop(0)]
            binmass = BinMass(mass,binpars)
            masslist.append(binmass)
            weightlist.append(NewTop.WeightList[iel].pop(0))
            
#Now combine all equal masses:
        newmasslist = []
        neweightlist = []
        while len(masslist) > 0:
            mass = masslist.pop(0)
            weight = weightlist.pop(0)
            while masslist.count(mass) > 0:
                w1 = weight
                w2 = weightlist.pop(masslist.index(mass))
                masslist.remove(mass)
                for k in weight.keys():
                    weight.update({k: w1[k] + w2[k]})
            newmasslist.append(mass)
            neweightlist.append(weight)
            
            
        for imass in range(len(newmasslist)):
            NewTop.B[0].ElList[iel].masses.append(newmasslist[imass][0])
            NewTop.B[1].ElList[iel].masses.append(newmasslist[imass][1])
            NewTop.WeightList[iel].append(neweightlist[imass])
            
            


#Create mass list with all mass elements in analysis:
    Allmasses = []
    for iel in range(len(NewTop.B[0].ElList)):
        for imass in range(len(NewTop.B[0].ElList[iel].masses)):
            mass = [NewTop.B[0].ElList[iel].masses[imass],NewTop.B[1].ElList[iel].masses[imass]]
            Allmasses.append(str(mass))
      
    if len(Allmasses) == 0: 
        return [{'mass' : [],'result' : 0., 'conditions' : []}]  #Empty result
    
#Remove repeated entries:
    Allmasses = set(Allmasses)
    Allmasses = [eval(x) for x in Allmasses]
    

#Get weight format and generate a zero weight dictionary:
    zeroweight = {}
    for iel in range(len(Analysis.Top.WeightList)):
        if len(Analysis.Top.WeightList[iel]) > 0:
            for wk in Analysis.Top.WeightList[iel][0].keys():                
                zeroweight.update({wk : 0.})
            break    
        
        
#Replace mass list for each element by a common mass list with
#zero weights when necessary
    for iel in range(len(NewTop.B[0].ElList)):
        neweightlist = []
        for imass in range(len(Allmasses)):
            neweightlist.append({})
            neweightlist[imass].update(zeroweight)
            mass = Allmasses[imass]
            for jmass in range(len(NewTop.B[0].ElList[iel].masses)):
                if mass == [NewTop.B[0].ElList[iel].masses[jmass],NewTop.B[1].ElList[iel].masses[jmass]]:
                    neweightlist[imass].update(NewTop.WeightList[iel][jmass])
                    break
        NewTop.B[0].ElList[iel].masses = []        
        NewTop.B[1].ElList[iel].masses = []        
        NewTop.WeightList[iel] = []
        for imass in range(len(Allmasses)):
            NewTop.B[0].ElList[iel].masses.append(Allmasses[imass][0])
            NewTop.B[1].ElList[iel].masses.append(Allmasses[imass][1])
            NewTop.WeightList[iel].append(neweightlist[imass])



#Evaluate result:
    result = EvalRes_aux(res,NewTop)
#Evaluate conditions:
    conditions = EvalRes_aux(Analysis.results[res],NewTop)
   
    output = []
    for imass in range(len(Allmasses)):        
        output.append({'mass' : Allmasses[imass], 'result' : result[imass], 'conditions' : conditions[imass]})
    return output    
    

#Evaluates string expression in instr using the elements and weights
#stored in InTop
def EvalRes_aux(instr,InTop):
        
    outstr = instr.replace(" ","")
#Get ordered list of elements:
    El = []
    iels = 0
    while "[[[" in outstr:  #String has element                
        st = outstr[outstr.find("[[["):outstr.find("]]]")+3] #Get duplet        
        ptclist = eltostr(st)     # Get particle list
#Syntax checks:
        for ib in range(2):
            for ptcL in ptclist[ib]:
                for ptc in ptcL:
                    if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                        print "EvalRes: Unknown particle ",ptc
                        return False
        outstr = outstr.replace(st,"El["+str(iels)+"]")  # Replace element
        El.append(ptclist)   #Store elements
        iels +=1
        
#Get list of expressions (separated by commas):
    outstrv = outstr.rsplit(",")
    
#Find elements in InTop corresponding to elements in El:
    ieq =[]
    for i in range(len(El)):
        for j in range(len(InTop.B[0].ElList)):
            AEl = [InTop.B[0].ElList[j].particles,InTop.B[1].ElList[j].particles]
            if El[i] == AEl:
                ieq.append(j)
                break
            
    if len(ieq) != len(El):
        print "EvalRes_aux: Error evaluating analysis"
        return False
    
#Loop through common masslist and replace element by its weight:
    Allres = []
    Allmasses = InTop.B[0].ElList[0].masses
    for imass in range(len(Allmasses)):
        resw = {}
        for w in InTop.WeightList[0][imass].keys():
            for j in range(len(El)):
                if len(InTop.WeightList[ieq[j]]) != len(InTop.WeightList[0]):
                    print "EvalRes_aux: Mismatch in Analysis.ElList"
                    return False
                El[j] = InTop.WeightList[ieq[j]][imass][w]
                
            eout = [Ceval(x,El) for x in outstrv]
            if len(eout) == 1: eout = eout[0]
            resw.update({w : eout})
        Allres.append([resw])
        
    return Allres      
                


#Converts an element string to a nested particle list or vice-versa    
def eltostr(invar):
    if type(invar) == type(list()):
        st = str(invar).replace("'","")
        st = st.replace(" ","")
        return st
    elif type(invar) == type(str()):
        st = invar.replace(" ","")
        st = st[st.find("[[["):st.find("]]]")+3]
        st_B = []
        st_B.append(st[2:st.find("]],[[")+1])
        st_B.append(st[st.find("]],[[")+4:st.find("]]]")+1])

        ptclist = [[],[]]
        for ib in range(2):
            while "[" in st_B[ib] or "]" in st_B[ib]:
                ptcs = st_B[ib][st_B[ib].find("[")+1:st_B[ib].find("],[")].split(",")
                
#Syntax check:                
                for ptc in ptcs:
                    if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                        print "eltostr: Unknown particle:",ptc
                        return False
                    
                ptclist[ib].append(ptcs)
                sptcs = str(ptcs).replace("'","")
                sptcs = str(sptcs).replace(" ","")
                st_B[ib] = st_B[ib].replace(sptcs,"",1)

        return ptclist
        

            
            
        
        
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

  
    
