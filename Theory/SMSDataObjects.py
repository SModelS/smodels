#Converts pdg number to particle name according to the dictionaries Rodd
# and Reven    
def ptype(pdg):
  p=int(pdg)
  if p in Rodd: return Rodd[p]
  if p in Reven: 
    return Reven[p]
  else:
    return False
  
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

#---------------Dictionaries:
Rodd={ 
1000021 : "gluino", 1000022: "N1", 1000023 : "N2", 1000025 : "N3", 1000035 : "N4", 1000024 : "C1", 1000037 : "C2", 1000039 : "gravitino",  1000001 : "squark", 1000002 : "squark", 1000003 : "squark", 1000004 : "squark", 2000001 : "squark", 2000002 : "squark", 2000003 : "squark", 2000004 : "squark", 1000005 : "sbottom", 2000005 : "sbottom", 1000006 : "stop", 2000006 : "stop", 1000011 : "slepton", 1000013 : "slepton", 1000015 : "stau", 2000011 : "slepton", 2000013 : "slepton", 2000015 : "stau", 1000012 : "sneutrino", 1000014 : "sneutrino", 1000016 : "sneutrino", 2000012 : "sneutrino", 2000014 : "sneutrino", 2000016 : "sneutrino", -1000021 : "gluino", -1000022: "N1", -1000023 : "N2", -1000025 : "N3", -1000035 : "N4", -1000024 : "C1", -1000037 : "C2", -1000039 : "gravitino",  -1000001 : "squark", -1000002 : "squark", -1000003 : "squark", -1000004 : "squark", -2000001 : "squark", -2000002 : "squark", -2000003 : "squark", -2000004 : "squark", -1000005 : "sbottom", -2000005 : "sbottom", -1000006 : "stop", -2000006 : "stop", -1000011 : "slepton", -1000013 : "slepton", -1000015 : "stau", -2000011 : "slepton", -2000013 : "slepton", -2000015 : "stau", -1000012 : "sneutrino", -1000014 : "sneutrino", -1000016 : "sneutrino", -2000012 : "sneutrino", -2000014 : "sneutrino", -2000016 : "sneutrino"  

}


Reven={
 25 : "higgs", -25: "higgs", 35 : "H0", -35 : "H0", 36 : "A0", -36 : "A0", 37 : "H+", -37 : "H-", 23 : "Z", -23 : "Z", 22 : "photon", -22 : "photon", 24 : "W+", -24 : "W-", 16 : "nu", -16 : "nu", 15 : "ta-", -15 : "ta+", 14 : "nu", -14 : "nu", 13 : "mu-", -13 : "mu+", 12 : "nu", -12 : "nu", 11 : "e-", -11 : "e+", 5 : "b", -5 : "b", 6 : "t+", -6 : "t-", 1 : "jet", 2 : "jet", 3 : "jet", 4 : "jet", 21 : "jet", -1 : "jet", -2 : "jet", -3 : "jet", -4 : "jet", -21 : "jet"
 }

PtcDic={ 
"e" : ["e+","e-"], "mu" : ["mu+", "mu-"], "ta" : ["ta+","ta-"], "l+" : ["e+","mu+"],"l-" : ["e-","mu-"],"l" : ["e-","mu-","e+","mu+"], "W" : ["W+","W-"], "t" : ["t+","t-"], "L+" : ["e+","mu+","ta+"], "L-" : ["e-","mu-","ta-"], "L" : ["e+","mu+","ta+","e-","mu-","ta-"]
}
class BElement:
  def __init__(self):
    self.masses = []
    self.particles = []
    self.momID = 0

  def __str__ ( self ):
    from Experiment.SMSUnits import rmvunit
    ret="particles=%s masses=%s" % \
       ( self.particles, [ rmvunit(x,"GeV") for x in self.masses ] )
    return ret
  
class EElement:
  def __init__(self):
    self.B = []
    self.weight = []

#Get global topology info from element structure  
  def getEinfo(self):
    vertnumb = []
    vertparts = []  
    for el in self.B:
      vertnumb.append(len(el.masses))
      vertparts.append([len(x) for x in el.particles])
      if len(vertparts[len(vertparts)-1]) == vertnumb[len(vertnumb)-1]-1:
        vertparts[len(vertparts)-1].append(0)  #Append 0 for stable LSP
    return {"vertnumb" : vertnumb, "vertparts" : vertparts}
    
  def __str__ ( self ):
    ret="Branch #1={{"+str(self.B[0])+"}}, Branch #2={{"+str(self.B[1])+"}}"
    return ret
    
class GTop:
  """ global topology. contains a list of elements, and
    the number of vertices, and .. FIXME andre? """
    
  def __init__(self):    
    self.vertnumb = []
    self.vertparts = []
    self.ElList = []

  def leadingElement ( self ):
    """ often, a topology carries only one element, so
      we have a special accessor for this """
    if len(self.ElList)==0: return None
    return self.ElList[0]

  def elements ( self ):
    return self.ElList

  def __str__(self):
    ret="number of vertices=%s number of vertex particles=%s" % \
        ( self.vertnumb, self.vertparts )
    return ret

  def checkConsistency ( self, verbose=False ):
    """ the number of vertices and insertions per vertex is 
      redundant information in a topology, so we can perform
      an internal consistency check """
    for element in self.ElList:
      info=element.getEinfo()
      if self.vertnumb!=info["vertnumb"]: 
        if verbose: print "[SMSmethods.py] inconsistent topology!!!"
        return False
      if self.vertparts!=info["vertparts"]: 
        if verbose: print "[SMSmethods.py] inconsistent topology!!!"
        return False
    if verbose: print "[SMSmethods.py] topology is consistent."
    return True

#Adds Eelement to ElList
#OBS: NewElement must have the correct branch ordering!
  def AddElement(self, NewElement):

#First get global topology info from NewElement:
    Einfo = NewElement.getEinfo()  
#Sanity checks:
    if Einfo["vertnumb"] != self.vertnumb or Einfo["vertparts"] != self.vertparts:
      print "AddElement: wrong element topology"
      return False      
#Append element to ElList:    
    self.ElList.append(NewElement)
    return True
  
   
class EAnalysis:  
  def __init__(self):
    self.label = ""
    self.Top = GTop()
    self.sqrts = 0
    self.lum = 0
    self.results = {}
    self.plots = {}
    self.run = ""
    self.masscomp = 0.2


#Given the constraints dictionary, automatically fill the element list with the
#elements corresponding to the strings in the dictionary, skipping repeated ones
  def GenerateElements(self):
     
    ListOfStrs = []
    vertnumb = self.Top.vertnumb
    vertparts = self.Top.vertparts
#Syntax check:    
    for k in range(len(vertnumb)):
      if len(vertparts[k]) != vertnumb[k]:
        print "GenerateElements: Inconsistent data: ninsertions=%d len(insertions)=%d for ``%s''." % ( vertnumb[k], len(vertparts[k]), self.Top )
        return False
    
#Get all element strings:    
    inelements = self.results.items()
    
    for iii in range(len(inelements)):
      for ii in range(2):  
        con = inelements[iii][ii].replace(" ","")
        while "[" in con:  #String has element        
          st = con[con.find("[[["):con.find("]]]")+3] #Get duplet
          con = con.replace(st,"")  # Remove element duplet
          ptclist = eltostr(st)   # Get particle list
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
      NewEl = EElement()
      NewEl.B = [BElement(),BElement()]      
      for ib in range(2):
        NewEl.B[ib].particles = ptclist[ib]
      self.Top.ElList.append(NewEl)
  
  

#Check if the plots listed in results exist
  def GetPlots(self,verbose=True):
    from Experiment import SMSResults

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
