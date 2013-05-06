from SMSHelpers import addunit, rmvunit

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
    def __str__(self):
        return "particle pdg %d p=(%f,%f,%f,m=%f) st %d" % ( self.pdg, self.px, self.py, self.pz, self.mass, self.status )
      
class TElement:
    def __init__(self):
        self.masses = []
        self.particles = []
        self.momID = 0

    def __str__ ( self ):
        ret="[Element] particles=["
        for p in self.particles: 
          ret+"%s " % type(p)
        ret+="] masses=["
        for m in self.masses: ret+="%s " % rmvunit(m,"GeV")
        ret+="]"
        return ret
        
class DBranch:
    def __init__(self):
        self.vertnumb = 0
        self.vertparts = []
        self.ElList = []

    def __str__(self):
        ret="V#%d VP%s " % ( self.vertnumb, self.vertparts )
        for i in self.ElList: ret+="  el %s" % i
        return ret
        
class GTop:
    def __init__(self):
        self.B = [DBranch(),DBranch()]
        self.WeightList = []

    def __str__(self):
        return "[GTop] B1 "+str(self.B[0])+", B2 "+str(self.B[1])

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
        self.masscomp = 0.2


#Given the constraints dictionary, automatically fill the element list with the
#elements corresponding to the strings in the dictionary, skipping repeated ones
    def GenerateElements(self):
       
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
    
    

#Check if the plots listed in results exist
    def GetPlots(self,verbose=True):
        import SMSResults

        run = self.run        
        if run == "": run = None  #If run has not been defined, use latest
        for res in self.results.keys():
            if not self.plots.has_key(res):
                if verbose: print "SMSmethods.py: GetPlots: Plot for result",res,"in Analysis",self.label,"not found"
                topo = ""
                ana = []
            else:
                topo = self.plots[res][0]
                analyses = self.plots[res][1]
                
            for ana in analyses:
                if not SMSResults.exists(ana,topo,run):
                    if verbose: print "SMSmethods.py: GetPlots: Histogram for ",topo," in ",ana," for run ",run," not found"
            

        


#--------------------Methods:
    
#Event reader. Returns a list of MParticles
def getNextEvent(filename): 
    """ reads LHE file 'filename', returns a list of MParticles """
    
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

def getEventTop(PList, weight = {}):
    """ Reads particle list (PList) from event (as it is created by .getNextEvent)
    and generates Topology (GTop) with only one element in ElList, corresponding
    to the event. If the DoCompress and/or DoInvisible flags are on, also generate
    compressed topologies with small mass gaps and/or neutrinos emitted
    in the last step of the cascade ("effective LSP"). """
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

            El.masses.append(addunit(PList[mother].mass,'GeV'))
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
    added = True
    
#Keep compressing the topologies generated so far until no new compressions can happen:
    while added:        
        added = False
#Check for mass compressed topologies        
        if SMSglobals.DoCompress:
            for Top in ETopList:
                ETopComp = False
                ETopComp = MassCompTop(Top,SMSglobals.minmassgap)
                if ETopComp:
                    exists = False
                    for Topp in ETopList:
                        if EqualTops(Topp,ETopComp): exists = True 
                    if not exists:   #Avoid double counting (conservative)
                        SMSglobals.nComp += 1
                        ETopList.append(ETopComp)
                        added = True
            
#Check for invisible compressed topologies 
#(look for effective LSP, such as LSP + neutrino = LSP')          
        if SMSglobals.DoInvisible:
            for Top in ETopList:
                ETopInComp = False
                ETopInComp = InvCompTop(Top)
                if ETopInComp:
                    exists = False
                    for Topp in ETopList:
                        if EqualTops(Topp,ETopInComp): exists = True
                    if not exists:   #Avoid double counting (conservative)
                        SMSglobals.nInvis += 1
                        ETopList.append(ETopInComp)
                        added = True
    return ETopList


#If two masses in InTop are degenerate, return compressed topology
def MassCompTop(InTop,mingap):
    import copy    
       
    ETopComp = copy.deepcopy(InTop)  
#Loop over branches        
    for ib in range(2):
        if ETopComp.B[ib].vertnumb < 2: continue
#Remove all external particles between compressed masses          
        for ivert in range(ETopComp.B[ib].vertnumb-1):
            massA = ETopComp.B[ib].ElList[0].masses[ivert]
            massB = ETopComp.B[ib].ElList[0].masses[ivert+1]
            if abs(massA-massB) < mingap:
                ETopComp.B[ib].ElList[0].particles[ivert] = []
                ETopComp.B[ib].vertparts[ivert] = 0
                
#Remove all vertices and masses with zero particle emissions:
        while ETopComp.B[ib].vertparts.count(0) > 1:
            ivert = ETopComp.B[ib].vertparts.index(0)
            ETopComp.B[ib].vertnumb -= 1
            massA = ETopComp.B[ib].vertparts.pop(ivert) 
            massA = ETopComp.B[ib].ElList[0].masses.pop(ivert)
            massA = ETopComp.B[ib].ElList[0].particles.pop(ivert)    
                
            
    if not EqualTops(ETopComp,InTop):
        return ETopComp
    else:
        return False


            
            
#If InTop has an effective LSPs (LSP + neutrino = LSP'), return compressed topology
def InvCompTop(InTop):
    import copy    
    
      
    ETopComp = copy.deepcopy(InTop)
#Loop over branches        
    for ib in range(2):
        if ETopComp.B[ib].vertnumb < 2: continue
#Remove all external neutrinos        
        for ivert in range(ETopComp.B[ib].vertnumb):
            if ETopComp.B[ib].vertparts[ivert] > 0:
                ptcs = ETopComp.B[ib].ElList[0].particles[ivert]
                while ptcs.count('nu') > 0: ptcs.remove('nu')   #Delete neutrinos                
                ETopComp.B[ib].ElList[0].particles[ivert] = ptcs
                ETopComp.B[ib].vertparts[ivert] = len(ptcs)
                
                
#First first non-empty vertex at the end of the branch
        inv  = ETopComp.B[ib].vertnumb-1
        while inv > 0 and ETopComp.B[ib].vertparts[inv-1] == 0: inv -= 1
#Remove empty vertices at the end of the branch:
        ETopComp.B[ib].vertnumb = inv + 1
        ETopComp.B[ib].vertparts = InTop.B[ib].vertparts[0:inv]
        ETopComp.B[ib].vertparts.append(0)
        ETopComp.B[ib].ElList[0].particles = InTop.B[ib].ElList[0].particles[0:inv]
        ETopComp.B[ib].ElList[0].masses = InTop.B[ib].ElList[0].masses[0:inv+1]
        

    if not EqualTops(ETopComp,InTop):
        return ETopComp
    else:
        return False
          
         
                
                
     
                        

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
#If order = False, test both branch orderings (for an element doublet only) 
def EqualEls(El1,El2,order=True):
    
    if len(El1) != len(El2): return False
    
#If it is an array of elements:    
    if type(El1) == type(list()) and type(El1[0]) == type(TElement()):
        if type(El2) != type(list()): return False
        if type(El2[0]) != type(TElement()): return False
        if len(El1) == 2 and not order:
            ptcsA = [El2[0].particles,El2[1].particles]
            massA = [El2[0].masses,El2[1].masses]
            ptcs = [El1[0].particles,El1[1].particles]
            mass = [El1[0].masses,El1[1].masses]
            ptcs_b = [El1[1].particles,El1[0].particles]
            mass_b = [El1[1].masses,El1[0].masses]
            if ptcsA == ptcs and mass == massA: 
                return True
            elif ptcsA == ptcs_b and mass_b == massA: 
                return True
            else:
                return False
        else:       
            for i in range(len(El1)):
                if El1[i].particles != El2[i].particles: return False
                if El1[i].masses != El2[i].masses: return False
           
#If it is a single element:            
    else:
        if El1.particles != El2.particles: return False
        if El1.masses != El2.masses: return False

    return True

#Compare two elements or array of elements.
#If particles are similar and all masses equal, returns True,
#otherwise returns False
#If order = False, test both branch orderings (for an element doublet only)
def SimEls(El1,El2,order=True):
    
    if len(El1) != len(El2): return False
    
#If it is an array of elements:    
    if type(El1) == type(list()) and type(El1[0]) == type(TElement()):
        if type(El2) != type(list()): return False
        if type(El2[0]) != type(TElement()): return False
        
        if len(El1) == 2 and not order:
            ptcsA = [El2[0].particles,El2[1].particles]
            massA = [El2[0].masses,El2[1].masses]
            ptcs = [El1[0].particles,El1[1].particles]
            mass = [El1[0].masses,El1[1].masses]
            ptcs_b = [El1[1].particles,El1[0].particles]
            mass_b = [El1[1].masses,El1[0].masses]
            if SimParticles(ptcsA,ptcs) and mass == massA: 
                return True
            elif SimParticles(ptcsA,ptcs_b) and mass_b == massA: 
                return True
            else:
                return False
        else:
            for i in range(len(El1)):
                if not SimParticles(El1[i].particles,El2[i].particles): return False
                if El1[i].masses != El2[i].masses: return False
        
#If it is a single element:            
    else:
        if not SimParticles(El1.particles,El2.particles) or  El1.masses == El2.masses: return False
        
    return True


#Compares 2 particle names or 2 nested name arrays. Allows for dictionary labels
#(Ex: L = l, l+ = l, l = l-,...)
def SimParticles(ptype1,ptype2):

    ptype1v = [ptype1]
    ptype2v = [ptype2]
 
#First flatten nested arrays: 
    isNested = True
    while isNested:
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
        isNested = False
        for i in range(len(ptype1v)):
            if type(ptype1v[i]) == type(list()): isNested = True
            if type(ptype2v[i]) == type(list()): isNested = True

        
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
    """ Check if elements in SMSTop matches an entry in SMSTopList. If it does,
    add weight.  If the same topology exists, but not the same element, add
    element.  If neither element nor topology exist, add the new topology and
    all its elements """
    
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
                w1 = SMSTopList[equaltops].WeightList[equalels]
                w2 = SMSTop.WeightList[inew]                
                SMSTopList[equaltops].WeightList[equalels] = sumweights([w1,w2])
                    
    
                    
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
    import copy
    
    for itop in range(len(SMSTopList)):
        NewTop = SMSTopList[itop]    
#Check if topologies match:
        if not EqualTops(NewTop,Analysis.Top): continue
        
#Loop over (event) element list:
        for iel in range(len(NewTop.B[0].ElList)):
#Loop over analysis elements:
            for jel in range(len(Analysis.Top.B[0].ElList)):
                NewEl = [NewTop.B[0].ElList[iel],NewTop.B[1].ElList[iel]]
                weight = copy.deepcopy(NewTop.WeightList[iel])
                newptcs = [NewTop.B[0].ElList[iel].particles,NewTop.B[1].ElList[iel].particles]
                oldptcs = [Analysis.Top.B[0].ElList[jel].particles,Analysis.Top.B[1].ElList[jel].particles]
#Loop over nested mass list in Analysis element:
                added = False
                for imass in range(len(Analysis.Top.B[0].ElList[jel].masses)):
                    OldEl = [TElement(),TElement()]
                    for ib in range(2):
                        OldEl[ib].particles = Analysis.Top.B[ib].ElList[jel].particles
                        OldEl[ib].masses = Analysis.Top.B[ib].ElList[jel].masses[imass]
                        
                                            
#Check if elements match (with identical masses) for any branch ordering
                    if SimEls(NewEl,OldEl,order=False):
                        Analysis.Top.WeightList[jel][imass] = sumweights([Analysis.Top.WeightList[jel][imass],weight])
                        added = True
                        break   #To avoid double counting only add event to one mass combination
                    
                if not added:
#Check for both branch orderings, but only add one (if matches) to avoid double counting                    
                    ptcmatch = False
                    if SimParticles([newptcs[0],newptcs[1]],oldptcs):
                        ptcmatch = 1
                    elif SimParticles([newptcs[1],newptcs[0]],oldptcs):
                        ptcmatch = 2
                    if ptcmatch:
                        for ib in range(2):
                            Analysis.Top.B[ib].ElList[jel].masses.append(NewEl[(ib+ptcmatch-1)%2].masses)
                        Analysis.Top.WeightList[jel].append(weight)


#Definition of distance between two mass arrays
#If Dana is defined, use maximum distance in all analyses
def MassDist(mass1,mass2):
    """ definition of distance between two mass arrays """
    import SMSglobals, SMSgetlimit
    
    Dana = SMSglobals.DistAnalyses  #List of analyses to be used

#Get upper bounds for each mass:
    xmass1 = SMSgetlimit.GetPlotLimit(mass1,Dana[0],Dana[1],complain=True)
    xmass2 = SMSgetlimit.GetPlotLimit(mass2,Dana[0],Dana[1],complain=True)
    if xmass1==None or xmass1==False:
      print "[SMSmethods.MassDist] no limit for plot 1"
      return None
    if xmass2==None or xmass2==False:
      print "[SMSmethods.MassDist] no limit for plot 2"
      return None
    
    d = -1.
    for iana in range(len(xmass1)):
        x1 = rmvunit(xmass1[iana][1],'fb')
        x2 = rmvunit(xmass2[iana][1],'fb')
        if type(x1) == type(str()) or type(x1) == type(str()): continue  #Skip analysis error messages
        if x1 and x2:
            if xmass1[iana][0] != xmass2[iana][0]:    #Check if analysis label match
                print "MassDist: Error getting upper limit"
                return None            
            
            newd = 2.*abs(x1-x2)/(x1+x2)   #Relative distance in "upper limit space"
            d = max(d,newd)

            
    if d < 0.: return None   #Skip masses without an upper limit
    
   #If masses differ by more than 100%, do not define distance
    if abs(mass1[0][0]-mass2[0][0])/(mass1[0][0]+mass2[0][0]) > 0.5: 
        return None
    
    return d

#Definition of distance two clusters
#MD = square matrix of distances
def ClusterDist(cluster1,cluster2,MD):
    d = 0.
    if type(cluster1) != type(set()) or type(cluster2) != type(set()):
        print "ClusterDist: unknown format input"
        return False
        
    for ic in cluster1:
        for jc in cluster2:
            if MD[ic][jc] == None: return None
            d = max(d,MD[ic][jc])
    return d
    

#Test if a mass array is "good"
# = have similar branch masses if branch topologies are equal
# = have similar mother and LSP masses if branch topologies are different
#If it is, return an equivalent array with equal masses (= mass avg)
def GoodMass(mass,Distfunc,dmin):
    
    if mass[0] == mass[1]: return mass
    if len(mass[0]) == len(mass[1]):
        mass1 = [mass[0],mass[0]]
        mass2 = [mass[1],mass[1]]
        MD = Distfunc(mass1,mass2)
        if MD == None or MD > dmin: 
            return False
        else:
            return MassAvg([mass1,mass2],"harmonic")
    else:
        mass1 = mass
        mass2 = mass
        mass1[1][0] = mass1[0][0]     #Force mothers and daughters to be equal in each branch
        mass1[1][len(mass1)-1] = mass1[0][len(mass1)-1]
        mass2[0][0] = mass2[1][0]
        mass2[0][len(mass2)-1] = mass2[1][len(mass2)-1]
        MD = Distfunc(mass1,mass2)
        if MD == None or MD > dmin: 
            return False
        else:
            return MassAvg([mass1,mass2],"harmonic")
    


#Cluster algorithm (generic for any type of object, as long as the distance function is given):
def DoCluster(objlist,Distfunc,dmin):
    import copy
    
    MD = []
#Compute distance matrix
    for i in range(len(objlist)):
        line = []
        for j in range(len(objlist)):
            if j >= i:
                line.append(Distfunc(objlist[i],objlist[j]))
            else:
                line.append(addunit(0.,'GeV'))
        MD.append(line)
        
    for i in range(len(objlist)):
        for j in range(len(objlist)):
            if j < i: MD[i][j] = MD[j][i]

#Begin clustering
    ClusterList = []
    for i in range(len(objlist)):
        cluster = set([])
        for j in range(len(objlist)):
            if MD[i][j] == None: continue
            if MD[i][j] <= dmin: cluster.add(j) 
        if not cluster in ClusterList: ClusterList.append(cluster)   #Zero level clusters (individual masses)


    FinalCluster = []
    newClusters = [0]
    while len(newClusters) > 0:
        newClusters = []    
        for cluster in ClusterList:
            split = False
            if len(cluster) > 2:
                for i in cluster:
                    ClDist = ClusterDist(set([i]),cluster,MD)
                    if  ClDist == None or ClDist > dmin:
                        newcluster = copy.deepcopy(cluster)
                        newcluster.remove(i)
                        split = True
                        if not newcluster in newClusters:
                            newClusters.append(newcluster)
                    
            if not split and not cluster in FinalCluster: FinalCluster.append(cluster)            
                    
        ClusterList = newClusters
        if len(ClusterList) > 1000:    #Check for oversized list of cluster (too time consuming)
            print "DoCluster: Error clustering. ClusterList >",len(ClusterList)
            return False
        
                
#Clean up clusters
    FinalCluster = FinalCluster + ClusterList
    i = 0
    for i in range(len(FinalCluster)):
        clusterA = FinalCluster[i]
        for j in range(len(FinalCluster)):
            clusterB = FinalCluster[j]
            if i != j and clusterB.issubset(clusterA):
                FinalCluster[j] = set([])
            
    while FinalCluster.count(set([])) > 0: FinalCluster.remove(set([]))
    
    return FinalCluster       
                          
                        
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
            if rmvunit(mass,'GeV') == 0.:
                print "MassAvg: Zero mass!"
                return False
            if rmvunit(mass,'GeV') < 0.:
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

#Sum a list of weights
def sumweights(wlist):
    
    neweight = {}
    for wk in wlist[0].keys(): neweight.update({wk : addunit(0.,'fb')})
    for wk in wlist[0].keys():
        wsum = addunit(0.,'fb')
        for weight in wlist: wsum = wsum + weight[wk]
        neweight.update({wk : wsum})
        
    return neweight
            
                    

#Evaluate theoretical predictions for the analysis result and conditions:
def EvalRes(res,Analysis,uselimits = False):
    import SMSglobals, copy

    output = []
    if not Analysis.plots.has_key(res) or not Analysis.results.has_key(res):
        print "EvalRes: Wrong analysis input"
        return False
    
#Get minimum distance parameter (it can be analysis dependent)
    dmin = Analysis.masscomp
#List of analyses and anlysis itself (necessary to compute distances)    
    analyses = Analysis.plots[res]
    SMSglobals.DistAnalyses = [analyses,Analysis]
    
#Create a mass list with all masses appearing in the analysis which have similar branch masses:
    Goodmasses = []
    Top = copy.deepcopy(Analysis.Top)
    for iel in range(len(Top.B[0].ElList)):        
        for imass in range(len(Top.B[0].ElList[iel].masses)):
            mass = [Top.B[0].ElList[iel].masses[imass],Top.B[1].ElList[iel].masses[imass]]
            gmass = GoodMass(mass,MassDist,dmin)
            if gmass:
                 Top.B[0].ElList[iel].masses[imass] = gmass[0]
                 Top.B[1].ElList[iel].masses[imass] = gmass[1]
                 if not gmass in Goodmasses: Goodmasses.append(gmass)
                 
#Cluster masses:
    MCluster = DoCluster(Goodmasses,MassDist,dmin)
   
#Loop over clusters to evaluate constraints and conditions inside each cluster
    for cluster in MCluster:
#Get masses in cluster
        masscluster = []
        for ic in cluster: masscluster.append(Goodmasses[ic])
        
#Shrink Topology elements to cluster:
        NewTop = ClusterTop(Top,masscluster)

#Now NewTop contains only elements with a common mass (replaced by the average mass)
#Evaluate result inside cluster
        result = Eval_cluster(res,NewTop)
#Evaluate conditions
        conditions = Eval_cluster(Analysis.results[res],NewTop)
        
#Save cluster result   
        mavg = [NewTop.B[0].ElList[0].masses,NewTop.B[1].ElList[0].masses]
        
#Check if average mass is inside the cluster (exp. limit for average mass ~ exp. limit for individual masses):
        davg = -1.
        for mass in masscluster:
            davg = max(davg,MassDist(mass,mavg))
        if davg == -1. or davg > dmin:
            print "EvalRes: Wrong clustering"
            continue
                
                    
        output.append({'mass' : mavg, 'result' : result, 'conditions' : conditions})
        
    return output    



#Given a Topology and the clustered masses, return a new topology 
#with only the elements belonging to the cluster
#(all masses are replaced by the average mass and elements with
#equivalent masses have their weights combined)
def ClusterTop(Top,masscluster):
    import copy

#Compute average mass in cluster
    mavg = MassAvg(masscluster,"harmonic")
        
#Keep only elements which belong to the cluster
    NewTop = copy.deepcopy(Top)
    NewTop.B[0].ElList = []
    NewTop.B[1].ElList = []
    NewTop.WeightList = []
    for iel in range(len(Top.B[0].ElList)):
        for imass in range(len(Top.B[0].ElList[iel].masses)):
            mass = [Top.B[0].ElList[iel].masses[imass],Top.B[1].ElList[iel].masses[imass]]
            ptc = [Top.B[0].ElList[iel].particles,Top.B[1].ElList[iel].particles]
            weight = Top.WeightList[iel][imass]

#If mass is in cluster, add element to NewTop:                
            if mass in masscluster:
                Elm = [TElement(),TElement()]
                for ib in range(2):
                    Elm[ib].masses = mavg[ib]
                    Elm[ib].particles = ptc[ib]
                match = False    
                for iel2 in range(len(NewTop.B[0].ElList)):    
                    ptcB = [NewTop.B[0].ElList[iel2].particles,NewTop.B[1].ElList[iel2].particles]
                    if ptcB == ptc:
                        match = True
                        oldweight = NewTop.WeightList[iel2]
                        NewTop.WeightList[iel2] = sumweights([oldweight,weight])
                        break
                    
                if not match:
                    NewTop.AddElement(Elm,weight)
                        
    return NewTop                        


#Evaluates string expression in instr using the elements and weights
#stored in InTop
def Eval_cluster(instr,InTop):
        
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

#Generate zeroweight entry
    zeroweight = {}  
    for wk in InTop.WeightList[0].keys():
        zeroweight.update({wk : addunit(0.,'fb')})
        
#Find elements in InTop corresponding to elements in El and fill Elw with the respective weights: 
    Elw = []       
    for i in range(len(El)):
        Elw.append(zeroweight)
        for j in range(len(InTop.B[0].ElList)):
            AEl = [InTop.B[0].ElList[j].particles,InTop.B[1].ElList[j].particles]
            if El[i] == AEl:                
                Elw[i] = InTop.WeightList[j]
                break

#Evaluate the instr expression (condition or constraint) for each weight entry:
    result = {}        
    Els = []
    if len(Elw) > 0:    
        for w in Elw[0].keys():
            Els = [weight[w] for weight in Elw]            
            eout = [Ceval(x,Els) for x in outstrv]        
            if len(eout) == 1: eout = eout[0]
            result.update({w : eout})
    else:
        eout = [Ceval(x,Els) for x in outstrv]
        result = eout
        
        
    return result      


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

  

