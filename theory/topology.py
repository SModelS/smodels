#!/usr/bin/env python

"""
.. module:: theory.Topology
   :synopsis: missing
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
    
"""
from theory import crossSection
import logging
from theory.element import Element
from printer import Printer
from prettytable import PrettyTable
from tools import SMSPrettyPrinter

logger = logging.getLogger(__name__)


class Topology(object):
    """
    Abstract class. DO NOT INSTANTIATE!
    """


    def __init__(self,elements=None):
        """
        Constructor. If elements is defined, create the topology from it. 
        If elements it is a list, all elements must share a common global topology.
        """
        self.vertnumb = []
        self.vertparts = []
        self.ElList = []
        
        if elements:
            if type(elements) == type(Element()): self.addElement(elements)              
            elif type(elements) == type([]):
                for element in elements: self.addElement(element)


    def __eq__(self,other):
        return self.isEqual(other)
    
    def __ne__(self,other):
        return not self.isEqual(other)
    
    
    def isEqual(self, other, order=False):
        """
        is this topology equal to other?
        Returns true if they have the same number of vertices and particles.
        If order=False and each topology has two branches, ignore branch ordering.
        """
        
        if type(self) != type(other): return False
        if order or len(self.vertnumb) != 2 or len(other.vertnumb) != 2:
            if self.vertnumb != other.vertnumb: return False
            if self.vertparts != other.vertparts: return False
            return True
        else:
            x1 = [self.vertnumb[0],self.vertparts[0]]
            x2 = [self.vertnumb[1],self.vertparts[1]]
            xA = [other.vertnumb[0],other.vertparts[0]]
            xB = [other.vertnumb[1],other.vertparts[1]]
            if x1 == xA and x2 == xB: return True
            if x1 == xB and x2 == xA: return True
            return False
    
    
    def checkConsistency ( self, verbose=False ):
        """
        the number of vertices and insertions per vertex is
        redundant information in a topology, so we can perform
        an internal consistency check
        """
        for element in self.ElList:
            info=element.getEinfo()
            if self.vertnumb!=info["vertnumb"]:
                if verbose: print "[topology.py] inconsistent topology!!!"
                return False
            if self.vertparts!=info["vertparts"]:
                if verbose: print "[topology.py] inconsistent topology!!!"
                return False
        if verbose: print "[topology.py] topology is consistent."
        return True
    
    def getTinfo(self):
        """
        Returns a dictionary with the topology number of vertices and vertparts
        """
        
        return {'vertnumb' : self.vertnumb, 'vertparts' : self.vertparts}
    
    def describe(self):
        """ a lengthy description """
        ret="number of vertices=%s number of vertex particles=%s number of elements=%d" % \
              ( self.vertnumb, self.vertparts, len(self.ElList) )
        return ret
    
    
    def elements ( self ):
        return self.ElList

     
    def addElement(self, newelement):
        """
        Adds an Element object to the ElList.
        For all the pre-existing elements which match the new element, add weight.
        If no pre-existing elements match the new one, add it to the list.
        If order=False, try both branch orderings        
        """
                
#If the topology info has not been set yet, set it using the element info
        if not self.vertparts: self.vertparts = newelement.getEinfo()["vertparts"]
        if not self.vertnumb: self.vertnumb = newelement.getEinfo()["vertnumb"]
        
#Check if newelement matches the topology structure:
        info = newelement.getEinfo()
        info_b = newelement.switchBranches().getEinfo()
        if info != self.getTinfo() and info_b != self.getTinfo():
            logger.warning('Element to be added does not match topology')
            return False
        
              
        added = False
#Include element to ElList:
        for element in self.ElList:
            if element == newelement:
                added = True
                element.weight.combineWith(newelement.weight)
                if element.getMothers() != newelement.getMothers():  #When combining elements with different mothers, erase mother info
                    element.branches[0].momID = None
                    element.branches[1].momID = None
                if element.getDaughters() != newelement.getDaughters(): #When combining elements with different daughters, erase daughter info
                    element.branches[0].daughterID = None
                    element.branches[1].daughterID = None
                
                
        if added: return True
#If element has not been found add to list (OBS: if both branch orderings are good, add the original one)
        tryelements = [newelement,newelement.switchBranches()]  #Check both branch orderings
        for newel in tryelements:
            info = newel.getEinfo()
            if info["vertnumb"] != self.vertnumb or info["vertparts"] != self.vertparts: continue
            self.ElList.append(newel)
            return True
        
#If element does not match topology info, return False
        return False
    
    def getTotalWeight(self):
        """ Returns the sum of all elements weights
        """
        if len(self.ElList) == 0: return None
        
        sumw = crossSection.XSectionList()        
        for element in self.ElList: sumw.combineWith(element.weight)
           
        return sumw    
 
    
class TopologyList(Printer):
    """
    Implements a list of topologies, knows how to correctly add a topology.
    """
    
    def __init__(self, topos=[]):
        """
        If topos are given, we add all of them sequentially.
        """
        self.topos = []
        for topo in topos: self.add(topo)

    def __len__ (self): 
        return len(self.topos)    

    def __getitem__(self,index):
        return self.topos[index]
    
    def __iter__ (self):
        return iter(self.topos)
    
    def __str__(self):
        s="TopologyList:\n" 
        for topo in self.topos:
            s+=str(topo)+"\n"
        return s

    def addList(self,List):
        for topo in List: self.add(topo)

    def describe(self):
        s="TopologyList:\n" 
        for topo in self.topos:
            s+=str(topo)+"\n"
        return s


    def add(self,newtopo):
        """
        Check if elements in newtopo matches an entry in self.topos. If it does,
        add weight.  If the same topology exists, but not the same element, add
        element.  If neither element nor topology exist, add the new topology and
        all its elements 

        :type topo: Topology    
        """
        
        topmatch = False        
        for itopo,topo in enumerate(self.topos):
            if topo == newtopo: topmatch = itopo
#If no pre-existing topologies match, append it to list of topologies
        if topmatch is False:
            self.topos.append(newtopo)
        else:
            for newelement in newtopo.ElList: self.topos[topmatch].addElement(newelement)
            
    def getTotalWeight(self):
        """ Returns the sum of all topologies total weights
        """
        sumw = crossSection.XSectionList()
        for topo in self:
            topoweight = topo.getTotalWeight()
            if topoweight: sumw.combineWith(topoweight)       
        return sumw
    
    def getElements(self):
        """ Returns a list with all the elements in all the topologies """
        
        elements = []
        for top in self.topos: elements.extend(top.ElList)
        return elements
    
    def prepareData(self):
        output = ""
        
        printer=SMSPrettyPrinter.SMSPrettyPrinter()
        EvTop_table = PrettyTable(["Topology","#Vertices", "#Insertions", "#Elements", "Sum of weights"])
        EvElement_table = PrettyTable(["Topology","Element","Particles B[0]","Particles B[1]", "Masses B[0]","Masses B[1]","Element Weight"])
        
        eltot = 0
        # totweight = []
        # Print Results:
        # for i in range(len(SMSTopList)):
        for i,topo in enumerate(self):
            sumw = topo.getTotalWeight().getDictionary()
            EvTop_table.add_row([i,topo.vertnumb,topo.vertparts,len(topo.ElList),SMSPrettyPrinter.wrap(printer.pformat(sumw),width=30)])
            eltot += len(topo.ElList)
        
        #Print element list for Topology[i]:  
            if i == 0:
                for j,el in enumerate(topo.ElList):
                    if el.getParticles() != [[['b','b']],[['b','b']]]:
                        continue
                    EvElement_table.add_row([i,j,el.getParticles()[0],el.getParticles()[1],SMSPrettyPrinter.wrap(printer.pformat(el.getMasses()[0]),width=25),SMSPrettyPrinter.wrap(printer.pformat(el.getMasses()[1]),width=25),SMSPrettyPrinter.wrap(printer.pformat(el.weight.getDictionary()),width=30)])
                EvElement_table.add_row(["---","---","---","---","---","---","---"])  
                     
        output += "Number of Global topologies = " + str(len(self)) + "\n"     
        output += str(EvTop_table) + "\n"
        output += "Total Number of Elements = " + str(eltot) + "\n"
        output += "Total weight = " + str(self.getTotalWeight()) + "\n"
        # output += EvElement_table + "\n"
        
        return output
        
