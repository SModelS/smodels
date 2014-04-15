#!/usr/bin/env python

"""
.. module:: Theory.Topology
   :synopsis: missing
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
    
"""
from ParticleNames import simParticles
from Tools.PhysicsUnits import addunit
from Theory import ClusterTools, TheoryPrediction, CrossSection
from AuxiliaryFunctions import getelements, eltonum, Ceval
import copy
import logging
from Theory.element import AElement, CElement, EElement, MassWeight
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Topology(object):
    """
    Abstract class. DO NOT INSTANTIATE!
    """


    def __init__(self):
        """
        Constructor
        """
        self.vertnumb = []
        self.vertparts = []
        self.ElList = []


    def __eq__ ( self, other ):
        return self.isEqual ( other )
    
    
    def isEqual ( self, Top2, order=False ):
        """
        is this topology equal to Top2?
        Returns true if they have the same number of vertices and particles.
        If order=False and each topology has two branches, ignore branch ordering.
        """
        if order or len(self.vertnumb) != 2 or len(Top2.vertnumb) != 2:
            if self.vertnumb != Top2.vertnumb: return False
            if self.vertparts != Top2.vertparts: return False
            return True
        else:
            x1 = [self.vertnumb[0],self.vertparts[0]]
            x2 = [self.vertnumb[1],self.vertparts[1]]
            xA = [Top2.vertnumb[0],Top2.vertparts[0]]
            xB = [Top2.vertnumb[1],Top2.vertparts[1]]
            if x1 == xA and x2 == xB: return True
            if x1 == xB and x2 == xA: return True
            return False
        
    
    def leadingElement(self):
        """
        often, a topology carries only one element, so
          we have a special accessor for this
        """
        if len(self.ElList) == 0:
            return None        
        elif len(self.ElList) > 1:
            logger.warning("ElList has " + len(self.ElList) + "elements!")
        
        return self.ElList[0]
    
    
    def checkConsistency ( self, verbose=False ):
        """
        the number of vertices and insertions per vertex is
        redundant information in a topology, so we can perform
        an internal consistency check
        """
        for element in self.ElList:
            info=element.getEinfo()
            if self.vertnumb!=info["vertnumb"]:
                if verbose: print "[SMSDataObjects.py] inconsistent topology!!!"
                return False
            if self.vertparts!=info["vertparts"]:
                if verbose: print "[SMSDataObjects.py] inconsistent topology!!!"
                return False
        if verbose: print "[SMSDataObjects.py] topology is consistent."
        return True
    
    
    def describe(self):
        """ a lengthy description """
        ret="number of vertices=%s number of vertex particles=%s number of elements=%d" % \
              ( self.vertnumb, self.vertparts, len(self.ElList) )
        return ret
    
    
    def elements ( self ):
        return self.ElList


class GTop(Topology):
    
    
    def addElement(self, NewElement):
        """
        missing
        """
        #First get global topology info from NewElement:
        Einfo = NewElement.getEinfo()
        #Sanity checks:
        if Einfo["vertnumb"] != self.vertnumb or Einfo["vertparts"] != self.vertparts:
            print "[GTop.addElement] wrong element topology"
            return False
        #Append element to ElList:
        self.ElList.append(NewElement)
        return True
    
    
    def massCompressedTopology ( self, mingap ):
        """
        if two masses in this topology are degenerate, create
        a compressed copy of this topology
        """
        from Tools.PhysicsUnits import rmvunit
        mingap=rmvunit ( mingap, "GeV" )
        ETopComp = copy.deepcopy(self)
        #Loop over branches
        for ib in range(len(ETopComp.vertnumb)):
            if ETopComp.vertnumb[ib] < 2: continue
            #Remove all external particles between compressed masses
            for ivert in range(ETopComp.vertnumb[ib]-1):
                massA = rmvunit(ETopComp.ElList[0].B[ib].masses[ivert],"GeV")
                massB = rmvunit(ETopComp.ElList[0].B[ib].masses[ivert+1],"GeV")
                if abs(massA-massB) < mingap:
                    ETopComp.ElList[0].B[ib].particles[ivert] = []
                    ETopComp.vertparts[ib][ivert] = 0
      
            #Remove all vertices and masses with zero particle emissions:
            while ETopComp.vertparts[ib].count(0) > 1:
                ivert = ETopComp.vertparts[ib].index(0)
                ETopComp.vertnumb[ib] -= 1
                massA = ETopComp.vertparts[ib].pop(ivert)
                massA = ETopComp.ElList[0].B[ib].masses.pop(ivert)
                massA = ETopComp.ElList[0].B[ib].particles.pop(ivert)
      
        if not ETopComp.isEqual(self):
            return ETopComp
        else:
            return False
      
      
    def invisibleCompressedTopology ( self ):
        """
        missing
        """
        ETopComp = copy.deepcopy(self)
        #Loop over branches
        for ib in range(len(ETopComp.vertnumb)):
            if ETopComp.vertnumb[ib] < 2: continue
            #Remove all external neutrinos
            for ivert in range(ETopComp.vertnumb[ib]):
                if ETopComp.vertparts[ib][ivert] > 0:
                    ptcs = ETopComp.ElList[0].B[ib].particles[ivert]
                    while ptcs.count('nu') > 0: ptcs.remove('nu')   #Delete neutrinos
                    ETopComp.ElList[0].B[ib].particles[ivert] = ptcs
                    ETopComp.vertparts[ib][ivert] = len(ptcs)
                    #First first non-empty vertex at the end of the branch
            inv  = ETopComp.vertnumb[ib]-1
            while inv > 0 and ETopComp.vertparts[ib][inv-1] == 0: inv -= 1
            #Remove empty vertices at the end of the branch:
            ETopComp.vertnumb[ib] = inv + 1
            ETopComp.vertparts[ib] = self.vertparts[ib][0:inv]
            ETopComp.vertparts[ib].append(0)
            ETopComp.ElList[0].B[ib].particles = self.ElList[0].B[ib].particles[0:inv]
            ETopComp.ElList[0].B[ib].masses = self.ElList[0].B[ib].masses[0:inv+1]
            
        if not ETopComp.isEqual(self):
            return ETopComp
        else:
            return False
        
        
class ATop(Topology):
    """
    analysis global topology. contains a list of analysis elements (AElements),
    the number of vertices and the number of particle insertions in each vertex
    """

    def addEventElement(self,NewElement,sqrts):
        """
        Adds an event element (EElement) to the corresponding analysis elements (AElements) in ElList\
        The event element and analysis elements DO NOT need to have the same branch ordering\
        sqrts = analysis center of mass energy (necessary for adding only the relevant theoretical cross-sections)
        """
        
        
        #Get cross-section information    
        XsecsInfo=None
        try:
          XsecsInfo = CrossSection.XSectionInfo  #Check if cross-section information has been defined
        except:
          pass
        if not XsecsInfo:
          XsecsInfo = CrossSection.XSecInfoList()   #If not, define default cross-sections
          CrossSection.XSectionInfo = XsecsInfo
          log = logging.getLogger(__name__)
          log.warning ( "Cross-section information not found. Using default values" )

#Restrict weight to the cross-section labels corresponding to sqrts
        zeroweight = {}     
        for xsec in XsecsInfo.xsecs:    #If not only add the weight labels corresponding to the analysis sqrts
          if sqrts == xsec.sqrts: zeroweight[xsec.label] = addunit(0., 'fb')          
        if not zeroweight: return False    #Skip analyses with unwanted sqrts
        
    
        if type(NewElement) != type(EElement()):
            print "[SMSDataObjects.py] wrong input! Must be an EElement object"
            return False
    
        newparticles_a = [NewElement.B[0].particles,NewElement.B[1].particles]
        newparticles_b = [NewElement.B[1].particles,NewElement.B[0].particles]
    
        for OldElement in self.ElList:
            if type(OldElement) != type(AElement()):
                print "[SMSDataObjects.py] wrong input! Elements in ATop must be an AElement object"
                return False    
            oldparticles = OldElement.getParticleList()
            #Check if particles match
            if not simParticles(newparticles_a,oldparticles) and not simParticles(newparticles_b,oldparticles): continue
            #Format the new weight to the analysis-dependent format (remove weights which do not match the analysis format and add zero to missing weights)
            neweight = copy.deepcopy(NewElement.weight)
            for key in zeroweight.keys()+neweight.keys():
                if not neweight.has_key(key): neweight[key] = addunit(0.,'fb')
                if not zeroweight.has_key(key): neweight.pop(key)      
                #Check if masses match
            added = False
            OldEl = EElement(OldElement.ParticleStr)  #Create temporary EElement for easy comparison
            for massweight in OldElement.MassWeightList:
                OldEl.B[0].masses = massweight.mass[0]
                OldEl.B[1].masses = massweight.mass[1]
                OldEl.weight = massweight.weight
                if NewElement.isSimilar(OldEl,order=False):
                    massweight.weight = ClusterTools.sumweights([OldEl.weight,neweight])
                    added = True
                    break   #To avoid double counting only add the event weight to one mass combination
                #If no identical mass was found, add entry to mass-weight dictionary
            if not added:
                newmass = [NewElement.B[0].masses,NewElement.B[1].masses]
                if not NewElement.isSimilar(OldEl,order=True,igmass=True): newmass = [newmass[1],newmass[0]] #Check for correct branch ordering
                newmassweight = MassWeight(newmass,neweight)
                OldElement.MassWeightList.append(newmassweight)    
    
        return True
    
    
class CTop(Topology):
    """
    cluster global topology. contains a list of cluster elements (CElements),
    the number of vertices and the number of particle insertions in each vertex and the cluster mass
    """    
    
    def __init__(self, Top=None, masscluster=None):
        self.clustermass = None
        super()
        
        if Top and masscluster:
            if type(Top) != type(ATop()):
                print "[SMSDataObjects.CTop] Input must be an analysis topology (ATop) !"
                return False
            if type(masscluster) != type([]):
                print "[SMSDataObjects.CTop] Input masscluster must be a list !"
                return False
            #Compute average mass in cluster
            self.clustermass = ClusterTools.MassAvg(masscluster,"harmonic")
            self.vertnumb = Top.vertnumb
            self.vertparts = Top.vertparts
            #Add elements to cluster
            for El in Top.ElList: self.addAnalysisElement(El,masscluster)
            
                 
    def addAnalysisElement(self, NewElement, masscluster):
        """
        Adds an analysis element (AElement) to the corresponding cluster elements (CElements) in ElList \
        Both the particles and branch orderings must be identical!
        """

        if type(NewElement) != type(AElement()):
            print "[SMSDataObjects.py] wrong input! Must be an AElement object"
            return False

        #Consistency checks:
        for OldEl in self.ElList:
            if type(OldEl) != type(CElement()):
                print "[SMSDataObjects.py] wrong input! Elements in ATop must be an AElement object"
                return False
            if OldEl.getEinfo() != NewElement.getEinfo(): return False
        if self.clustermass != ClusterTools.MassAvg(masscluster,"harmonic"):
            print "[SMSDataObjects.py] wrong masscluster input!"
            return False

        newparticles = NewElement.getParticleList()

        for massweight in NewElement.MassWeightList:
            mass = massweight.mass
            weight = massweight.weight
            if not mass in masscluster: continue
            #If mass is in cluster, add element to Cluster Topology
            match = False
            for OldEl in self.ElList:
                oldparticles = OldEl.getParticleList()
                if simParticles(oldparticles,newparticles,useDict=False):
                    match = True
                    oldweight = OldEl.Weight
                    OldEl.Weight = ClusterTools.sumweights([oldweight,weight])
                    break

            if not match:
                NewEl = CElement(NewElement.ParticleStr,weight)        
                self.ElList.append(NewEl)
        return True
    

    def evaluateCluster(self,results):
        """
        Evaluates the constraints and conditions in results using the
        respective theoretical cross section predictions for each element in the cluster topology. 

        :type results: a dictionary with the constraints and the conditions \
        of the experiment results

        :returns: an XSecPredictionForCluster object
        """
        #To store the result:
        ClusterResult = TheoryPrediction.XSecPredictionForCluster()

        #Get constraints and conditions:
        consts = results.keys()
        if len(consts) > 1:
            print "evaluateCluster: Analysis contains more than one entry"
            return False    
        conds = results.values()[0]
        if ";" in conds or "Csim" in conds or "Cgtr" in conds or "cSim" in conds or "cGtr" in conds:
            conds = conds.rsplit(";")
        else:
            conds = conds.rsplit(",")

        #Get a list of all elements appearing in results:
        allEl = set(getelements(consts) + getelements(conds))

        #Generate zeroweight and list of numerical elements with zero weights
        zeroweight = {}
        for wk in self.ElList[0].Weight.keys(): zeroweight[wk]=addunit(0.,'fb')
        nEll = [zeroweight]*len(allEl)  
        #Build a dictionary to map the relevant elements to its respective numerical element and the element to its weight
        thdic = {}
        iel = 0
        for el in allEl:
            ptcsA = CElement(el).getParticleList()
            nel = "nEl["+str(iel)+"]"
            thdic[el] = nel
            for El in self.ElList:
                ptcsB = El.getParticleList()
                if simParticles(ptcsA,ptcsB,useDict=False): nEll[iel] = El.Weight
            iel += 1

        #Replace string elements by their respective numerical element (nEl):
        consts_num = {}
        conds_num = {}
        for const in consts: consts_num[const] = eltonum(const,thdic)
        for cond in conds: conds_num[cond] = eltonum(cond,thdic)

        #Loop over weights and
        const_res = {}
        cond_res = {}
        #Evaluate each constraint
        for ckey in consts_num.keys():
            const = consts_num[ckey]
            res = {}
            for weight in zeroweight.keys():
                nEl = []
                for el in nEll: nEl.append(el[weight])  #Select weight
                res[weight] = Ceval(const,nEl)
            const_res[ckey] = res
            #Evaluate each condition   
            for ckey in conds_num.keys():
                cond = conds_num[ckey]
                res = {}
                for weight in zeroweight.keys():
                    nEl = []
                    for el in nEll: nEl.append(el[weight])  #Select weight
                    res[weight] = Ceval(cond,nEl)
                cond_res[ckey] = res   

        ClusterResult.conditions_dic = cond_res
        ClusterResult.result_dic = const_res

        return ClusterResult
    
    
class TopologyList(object):
    """
    Implements a list of topologies, knows how to correctly add a topology.
    """
    
    def __init__(self, topos=[]):
        """
        If topos are given, we add all of them sequentially.
        """
        self.topos = []
        for topo in topos:
            self.add ( topo )


    def __len__ ( self ): 
        return len(self.topos)
    

    def __getitem__ ( self, n ):
        return self.topos[n]


    def addList ( self, List ):
        for topo in List: self.add ( topo )


    def __str__ ( self ):
        s="TopologyList:\n" 
        for topo in self.topos:
            s+=str(topo)+"\n"
        return s


    def describe( self ):
        s="TopologyList:\n" 
        for topo in self.topos:
            s+=str(topo)+"\n"
        return s


    def add ( self, topo ):
        """
        Check if elements in topo matches an entry in self.topos. If it does,
        add weight.  If the same topology exists, but not the same element, add
        element.  If neither element nor topology exist, add the new topology and
        all its elements 

        :type topo: GTop    
        """
        for (inew,element) in enumerate(topo.ElList): ## range(len(topo.ElList)):
            if len(element.B)<2:
                print "[SMSDataObjects.TopologyList] error: assumed at least two branches"
                continue
            NewEl_a = element
            NewEl_b = copy.deepcopy(NewEl_a)
            NewEl_b.B[1] = element.B[0]
            NewEl_b.B[0] = element.B[1]   #Check both orderings
            equaltops = -1
            equalels = -1
            i = -1
            while (equaltops < 0 or equalels < 0) and i < len(self.topos)-1:
                i += 1
                if topo.isEqual(self.topos[i],order=False):  #First look for matching topology
                    equaltops = i
                else: continue

            for j in range(len(self.topos[i].ElList)):  #Search for matching element
                OldEl = self.topos[i].ElList[j]
                if OldEl.isEqual(NewEl_a):
                    equalels = j
                    NewEl = NewEl_a
                    break
                elif OldEl.isEqual(NewEl_b):
                    equalels = j
                    NewEl = NewEl_b
                    break

            #If element exists, add weight:
            if equalels >= 0:
                if len(OldEl.weight) != len(NewEl.weight):
                    print "Wrong number of weights"
                else:
                    w1 = OldEl.weight
                    w2 = NewEl.weight
                    self.topos[equaltops].ElList[equalels].weight = ClusterTools.sumweights([w1,w2])             
                    #When combining elements, keep the smallest set of PDG mother IDs (not used in the analysis, only relevant to set a standard):
                    if min(abs(NewEl.B[0].momID),abs(NewEl.B[1].momID)) < min(abs(OldEl.B[0].momID),abs(OldEl.B[1].momID)):
                        for ib in range(2): self.topos[equaltops].ElList[equalels].B[ib].momID = NewEl.B[ib].momID
                
        #If topology and/or element does not exist, add:
        if equaltops == -1:
            self.topos.append(topo)
        elif equalels == -1:
            if topo.isEqual(self.topos[equaltops],order=True):
                NewEl = NewEl_a
            else:
                NewEl = NewEl_b
            if not self.topos[equaltops].addElement(NewEl):
                print "Error adding element"
                print '\n'
