"""
.. module:: auxiliaryFunctions
   :synopsis: A collection of functions used to evaluate fuzzy the conditions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory import crossSection
from smodels.tools.physicsUnits import pb, GeV, fb
import numpy as np
from scipy import stats
from collections import Iterable
import copy
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
from smodels.experiment.finalStateParticles import finalStates
import networkx as nx

#Get all finalStateLabels
finalStateLabels = finalStates.getValuesFor('label')

def massPosition(mass, txdata):
    """ Give mass position in upper limit space.    
    Use the analysis experimental limit data. 
    :param txdata: TxNameData object holding the data and interpolation   
    """
    xmass = txdata.getValueFor(mass)
    if type(xmass) != type(1.*pb):
        return None

    return xmass.asNumber(fb)


def distance(xmass1, xmass2):
    """
    Define distance between two mass positions in upper limit space.
    The distance is defined as d = 2*|xmass1-xmass2|/(xmass1+xmass2).
    
    
    :parameter xmass1: upper limit value (in fb) for the mass1
    :parameter xmass2: upper limit value (in fb) for the mass2
    :returns: relative mass distance in upper limit space     
    """
    if xmass1 is None or xmass2 is None:
        return None
    distanceValue = 2.*abs(xmass1 - xmass2) / (xmass1 + xmass2)
    if distanceValue < 0.:
        # Skip masses without an upper limit
        return None
    
    return distanceValue


def massAvg(massList, method='weighted', weights=None):
    """
    Compute the average mass of massList according to method.

    If method=weighted but weights were not properly defined,
    switch method to harmonic.    
    If massList contains a zero mass, switch method to mean.
    
    :parameter method: possible values: harmonic, mean, weighted
    :parameter weights: weights of elements (only for weighted average)
    
    """
    if not massList:
        return massList
    if massList.count(massList[0]) == len(massList):
        return massList[0]

    if method == 'weighted' and (not weights or len(weights) != len(massList)):
        method = 'harmonic'
    flatList = [ mass / GeV for mass in _flattenList(massList)]
    if method == 'harmonic' and 0. in flatList:
        method = 'mean'

    for mass in massList:
        if len(mass) != len(massList[0]) \
                or len(mass[0]) != len(massList[0][0]) \
                or len(mass[1]) != len(massList[0][1]):
            logger.error('Mass shape mismatch in mass list:\n' + str(mass) +
                         ' and ' + str(massList[0]))
            raise SModelSError()

    avgmass = copy.deepcopy(massList[0])
    for ib, branch in enumerate(massList[0]):
        for ival in enumerate(branch):
            vals = [ float(mass[ib][ival[0]] / GeV) for mass in massList]
            if method == 'mean':
                avg = np.mean(vals)
            elif method == 'harmonic':
                avg = stats.hmean(vals)
            elif method == 'weighted':
                weights = [ float(weight) for weight in weights ]
                avg = np.average(vals,weights=weights)                
            avgmass[ib][ival[0]] = float(avg)*GeV

    return avgmass


def cSim(*weights):
    """
    Define the auxiliar similar function.
    
    Return the maximum relative difference between any element weights of the
    list, normalized to [0,1].
    
    :returns: XSectionList object with the values for each label.
    
    """
    for weight in weights:
        if type(weight) != type(crossSection.XSectionList()):
            logger.error("Trying to evaluate non-xsection objects")
            raise SModelSError()

    # Make sure both xsec lists have the same entries (add zero xsecs for the
    # missing entries)
    infoList = []
    for weight in weights:
        for info in weight.getInfo():
            if not info in infoList:
                infoList.append(info)
    zeros = crossSection.XSectionList(infoList)
    for weight in weights:
        weight.combineWith(zeros)

    # Evaluate the inequality for each cross section info
    result = crossSection.XSectionList()
    for info in infoList:
        res = 0.
        xsecRes = crossSection.XSection()
        xsecRes.info = info
        for weightA in weights:
            for weightB in weights:
                a = weightA.getXsecsFor(info.label)[0].value.asNumber ( fb )
                b = weightB.getXsecsFor(info.label)[0].value.asNumber ( fb )
                if a + b == 0.:
                    continue
                res = max(res, abs(a - b) / abs(a + b))
        xsecRes.value = res
        result.add(xsecRes)

    return result


def cGtr(weightA, weightB):
    """
    Define the auxiliary greater function.
    
    Return a number between 0 and 1 depending on how much it is violated
    (0 = A > B, 1 = A << B).
    
    :returns: XSectioList object with the values for each label.
    
    """
    if type(weightA) != type(crossSection.XSectionList()) or \
            type(weightB) != type(crossSection.XSectionList()):
        logger.error("Trying to evaluate non-xsection objects")
        raise SModelSError()

    # Make sure both xsec lists have the same entries (add zero xsecs for the
    # missing entries)
    infoList = weightA.getInfo()
    for info in weightB.getInfo():
        if not info in infoList:
            infoList.append(info)
    if not infoList:
        # If there are no cross sections, can not evaluate
        return 'N/A'
    zeros = crossSection.XSectionList(infoList)
    weightA.combineWith(zeros)
    weightB.combineWith(zeros)

    # Evaluate the inequality for each cross section info
    result = crossSection.XSectionList()
    for info in infoList:
        a = weightA.getXsecsFor(info.label)[0].value / fb
        b = weightB.getXsecsFor(info.label)[0].value / fb
        xsecRes = crossSection.XSection()
        xsecRes.info = info
        if a + b == 0.:
            xsecRes.value = 'N/A'
        else:
            xsecRes.value = (abs(a - b) - (a - b)) / (2.*(a + b))
        result.add(xsecRes)

    return result


def _flattenList(inlist, dims=None):
    """
    Flatten a multi-dimensional nested list.
    
    Output ordering: [first level objects, second level objects, ...].    
    
    If dims == [], include dimensions of nested list to it. This is useful when
    comparing lists).
    
    """
    flat = []
    for item in inlist:
        if isinstance(item, Iterable) and not isinstance(item, str ):
            if not dims is None:
                dims.append(len(item))
            for x in _flattenList(item, dims):
                flat.append(x)
        else:
            flat.append(item)
    return flat
    
def index_bisect(inlist, el):
    """
    Return the index where to insert item el in inlist.
    inlist is assumed to be sorted and a comparison function (lt or cmp)
    must exist for el and the other elements of the list.
    If el already appears in the list, inlist.insert(el) will
    insert just before the leftmost el already there.  
    """

    lo = 0    
    hi = len(inlist)
    while lo < hi:
        mid = (lo+hi)//2
        if inlist[mid] < el: lo = mid+1
        else: hi = mid
    return lo


def elementsInStr(instring,removeQuotes=True):
    """
    Parse instring and return a list of elements appearing in instring.
    instring can also be a list of strings.
    
    :param instring: string containing elements (e.g. "[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]")
    :param removeQuotes: If True, it will remove the quotes from the particle labels.
                         Set to False, if one wants to run eval on the output.
    
    :returns: list of elements appearing in instring in string format
    
    """
    outstr = ""
    if isinstance(instring,str):
        outstr = instring
    elif isinstance(instring,list):
        for st in instring:
            if not isinstance(st,str):
                logger.error("Input must be a string or a list of strings")
                raise SModelSError()
            # Combine list of strings in a single string
            outstr += st
    else:
        raise SModelSError ( "syntax error in constraint/condition: ``%s''." \
              "Check your constraints and conditions in your database." % str(instring) )

    elements = []
    outstr = outstr.replace(" ", "")
    if removeQuotes:
        outstr = outstr.replace("'", "")
    elStr = ""
    nc = 0
    # Parse the string and looks for matching ['s and ]'s, when the matching is
    # complete, store element
    for c in outstr:
        delta = 0
        if c == '[':
            delta = -1
        elif c == ']':
            delta = 1
        nc += delta
        if nc != 0:
            elStr += c
        if nc == 0 and delta != 0:
            elements.append(elStr + c)
            elStr = ""
            # Syntax checks
            ptclist = elements[-1].replace(']', ',').replace('[', ',').\
                    split(',')       
            for ptc in ptclist:
                ptc = ptc.replace("'","")
                if not ptc:
                    continue          
                if not ptc in finalStateLabels:
                    raise SModelSError("Unknown particle. Add " + ptc + " to finalStateParticles.py")

    # Check if there are not unmatched ['s and/or ]'s in the string
    if nc != 0:
        raise SModelSError("Wrong input (incomplete elements?) " + instring)

    return elements
           
def stringToTree(info,finalState=None):
    """
    Creates a Tree (DiGraph object) object from a string in bracket notation.
    
    :parameter info: string describing the element in bracket notation
                     (e.g. [[[e+],[jet]],[[e-],[jet]]])
                     
    :parameter finalState: list containing the final state labels for each branch
                     (e.g. ['MET', 'HSCP'] or ['MET','MET']). If not defined it will
                     be assumed to be MET for all branches
    :return: Tree (nerworkX DiGraph object)
                         
    """


    elements = elementsInStr(info,removeQuotes=False)
    if not elements or len(elements) > 1:
        nel = 0
        if elements:
            nel = len(elements)
        logger.error("Malformed input string. Number of elements "
                      "is %d (expected 1) in: ``%s''", nel, info)
        raise SModelSError()
    
    branches = eval(elements[0].replace('[*]',"['InclusiveBranch']"))
    if not branches:
        logger.error("Malformed input string. Number of "
                      "branches is %d (expected 2) in: ``%s''",
                      len(branches), info)
        raise SModelSError()

    if finalState:
        if not isinstance(finalState,list) or len(finalState) != len(branches):
            logger.error("Number of final states (%i) does not match number of branches (%i)" 
                               %(len(finalState),len(branches)))
            raise SModelSError()
    else:
        finalState = ['MET']*len(branches)

    from smodels.theory.branch import InclusiveBranch
    #Create map with all required particle objects for building the graph:
    fstateDict = {'anyOdd' : finalStates.getParticlesWith(label='anyOdd')[0],
                  'PV' : finalStates.getParticlesWith(label='PV')[0],
                 'InclusiveBranch' : InclusiveBranch()}
    #Add final state to the end of each branch:    
    for ib,branch in enumerate(branches):
        if not branch:
            branches[ib] = finalState[ib]  #The PV decays directly to the final state
        else:
            branch.append(finalState[ib])
    for ptc in _flattenList(branches) + finalState:
        if not ptc in fstateDict:
            particle = finalStates.getParticlesWith(label=ptc)
            if not particle or len(particle) != 1:
                logger.error("Error retrieving particle %s from finalStateParticles. Is this particle uniquely defined?")
                raise SModelSError()
            particle = particle[0]
            fstateDict[ptc] = particle

    #Create element rooted tree (starting with primary vertex = zero node):
    tree = vertexListToTree(branches,particleDict=fstateDict,parentNode=0,parentParticle=fstateDict['PV'])

    return tree


def branchListToTree(bList,particleDict,parentNode,parentParticle):
    """
    Creates a Tree (DiGraph object) from a nested list.
    The list should correspond to a series of sequential decays
    of 'anyOdd' particles (e.g.  [[a,b],[c]] =  X -> a+b+Y, Y->c) 
    If one of the entries of the list is a string, it is assumed
    to correspond to a final state to the previous odd particle
    (e.g. [[a,b],[c],MET] =  X -> a+b+Y, Y->c+MET, [[a,b],MET,[c]] =  X -> a+b+Y+MET, Y->c).
    A branch with a single particle should be represented by [X]
    (e.g. [HSCP] -> PV->HSCP).
    
    :parameter bList: List object containing a list of vertices or vertices and particles
                     (e.g. [[e+],[jet]],[[e+],[jet],MET])
    :parameter particleDict: Dictionary used to map strings to Particle objects.
    :parameter parentNode: Name of the parent node (integer). First node in the sequential decay
    :parameter parentParticle: Particle object corresponding to the first particle in the sequential decay.
                     
                     
    :return: Tree (nerworkX DiGraph object)
                         
    """    
    
    if not isinstance(bList,list):
        logger.error("Malformed bracket expression. Expecting a list and received %s" %str(bList))
        raise SModelSError()

    #Create tree
    T = nx.DiGraph()
    #The first entry in bList should correspond to the parent decay
    node = parentNode
    momNode = parentNode
    #Create a sequential Tree (each entry in bList is added in a sequence of levels)
    for vList in bList:        
        if isinstance(vList,list): #entry correspond to a vertex:
            #Each entry corresponding to an odd-particle in the vertex            
            if node != parentNode:
                T.add_node(node,particle=particleDict['anyOdd'])
                T.add_edge(momNode,node)
            else:
                T.add_node(node,particle=parentParticle)
            vTree = vertexListToTree(vList,particleDict,parentNode=node,
                                     parentParticle=T.nodes[node]['particle'])
            T = nx.compose(T,vTree)
            momNode = node
        elif isinstance(vList,str):  #entry correspond to an odd particle (e.g. a Z2-odd final state)
            if not vList in particleDict:
                logger.error("Particle %s not found in dictionary" %vList)
                raise SModelSError()
            ptc = particleDict[vList]
            T.add_node(node,particle=ptc)
            if node != parentNode:
                T.add_edge(momNode,node)            
        #Every sequential node added corresponds to going one level down the tree
        #Add digit everytime you go down one level:
        node *= 10
    return T
        

def vertexListToTree(vList,particleDict,parentNode,parentParticle):
    """
    Creates a Tree (DiGraph object) from a list.
    The list should correspond to a vertex (list of daughters from the same parent)
    (e.g.  [a,b] =  parentParticle -> a+b) .
    If one of the entries of the list is another list, it is assumed to correspond
    to a new (nested) branch (e.g. [ [[a,b]] , c ] -> parent -> X(->a+b) + c)
    
    :parameter vList: List object containing a list of daughters belonging to the same decay vertex
                     (e.g. [e+,jet,MET])
    :parameter particleDict: Dictionary used to map strings to Particle objects.
    :parameter parentNode: Name of a parent node (integer)
    :parameter parentParticle: Particle object corresponding to the parent of the vertex
                     
    :return: Tree (nerworkX DiGraph object)
                         
    """        
    
    if not isinstance(vList,list):
        logger.error("Malformed bracket expression. Expecting a list and received %s" %str(vList))
        raise SModelSError()
                     
                     
    #Create tree
    T = nx.DiGraph()
    #Always include the parent node in the graph:
    momNode = parentNode
    T.add_node(parentNode,particle=parentParticle)
    #A vertex is always one level below the parent. Add digit everytime you go down one level:
    node = 10*parentNode
    #Create a parallel Tree (each entry in vList is added at the same level/same parent)
    for p in vList:
        node += 1 # Add unit everytime you add a sibling (node at the same tree level)
        if isinstance(p,str):  #Entry in vertex is a Z2-even particle
            if not p in particleDict:
                logger.error("Particle %s not found in dictionary" %p)
                raise SModelSError()
            ptc = particleDict[p]
            T.add_node(node,particle=ptc)
            T.add_edge(momNode,node)
        elif isinstance(p,list): #Entry in vertex is another branch
            T.add_node(node,particle=particleDict['anyOdd'])
            T.add_edge(momNode,node)            
            branchTree = branchListToTree(p,particleDict,parentNode=node,
                                          parentParticle=T.nodes[node]['particle'])
            T = nx.compose(T,branchTree)
    return T
            
def fromTreeToList(T,node=None):
    """
    Convert a Tree to a nested list with the Z2 even
    final states. The Z2 odd final states (e.g. 'MET', 'HSCP') are
    not included.
    
    :param T: DiGraph object
    :param node: Name of node. If None it will start with the primary node
    
    :return: Nested list with Z2-even particle objects (e.g. [[[e-,mu],[L]],[[jet]]])
    """

    if not isinstance(T,nx.DiGraph):
        raise SModelSError("Input must be a DiGraph object.")
    
    #Get the primary vertex:
    if node is None:
        node = getTreeRoot(T)
        
    dList = []
    children = list(T[node])
    if children:
        for n in children:
            if list(T[n]):
                dList.append(fromTreeToList(T,n))
            else:
                if T.nodes[n]['particle'].Z2parity == 'even':
                    dList.append(T.nodes[n]['particle'])

    return dList


def getTopologyName(T,node=None):
    """
    Recursively compute the topology (canonical) name for the Tree T.
    Returns the name in integer form.
    
    :param T: Tree (networkX DiGraph object)
    :param node: Node to get the name for. If None, it will use the root
    
    :return: Integer representing the Tree
    """
    
    if not isinstance(T,nx.DiGraph):
        raise SModelSError("Input must be a DiGraph object.")
    
    if node is None:
        node = getTreeRoot(T)

    children = list(T[node])
    if not children:
        return 10
    else:
        tp = sorted([getTopologyName(T,n) for n in children])
        tpStr = '1'+"".join(str(c) for c in tp)+'0'
        return int(tpStr)
    
def getNodeLevelDict(T):
    """
    Return a dictionary with the nodes in each level of the tree T
    (e.g. {0 : [PV], 1 : [A,B],...})
    
    :param T: DiGraph object
    :return: Dictionary with the levels as keys and a list of nodes in each level as values.     
    """         

    if not isinstance(T,nx.DiGraph):
        raise SModelSError("Input must be a DiGraph object.")

    root = getTreeRoot(T)
    levelNodes = {}
    for n in T.nodes():
        d = nx.shortest_path_length(T,root,n) #Distance to root
        if not d in levelNodes:
            levelNodes[d] = [n]
        else:
            levelNodes[d].append(n)
            
    return levelNodes

def getTreeRoot(T):
    """
    Get the root node (primary vertex) of a tree T.
    
    :param T: DiGraph object
    :return: node
    """

    if not isinstance(T,nx.DiGraph):
        raise SModelSError("Input must be a DiGraph object.")

    root = [n for n in T.nodes() if T.in_degree[n] == 0]
    if len(root) != 1:
        raise SModelSError("Malformed Tree, %i root(s) have been found." %len(root))
    
    return root[0]
    
            
    
def sortParticleList(ptcList):
    """
    sorts a list of particle or particle list objects by their label
    :param ptcList: list to be sorted containing particle or particle list objects
    :return: sorted list of particles
    """
    
    newPtcList = sorted(ptcList, key=lambda x: x.label) 
    return newPtcList    


def getValuesForObj(obj, attribute):
    """
    Loops over all attributes in the object and in its attributes
    and returns a list of values for the desired attribute:
    
    :param obj: Any object with a __dict__ attribute
    :param attribute: String for the desired attribute
    
    :return: List with unique attribute values. If the attribute is not found, returns empty list.
    """
    
    values = []
    try:
        objDict = obj.__dict__.items()
    except:
        return values        
    
      
    for attr,value in objDict:
        if attribute == attr:
            values += [value]
        elif isinstance(value,Iterable):
            values += [getValuesForObj(v,attribute) for v in value]
        else:
            values += getValuesForObj(value,attribute)
    
    values =  list(filter(lambda a: (not isinstance(a,list)) or a != [], values))
    values = _flattenList(values)
    uniqueValues = [v for n,v in enumerate(values) if v not in values[:n]]
    
    return uniqueValues


def getAttributesFrom(obj):
    """
    Loops over all attributes in the object and return a list
    of the attributes.
    
    :param obj: Any object with a __dict__ attribute
    
    :return: List with unique attribute labels.
    """
    
    attributes = []
    try:
        objDict = obj.__dict__.items()
    except:
        return attributes        
    
      
    for attr,value in objDict:
        attributes.append(attr)
        if isinstance(value,list):
            attributes += [getAttributesFrom(v) for v in value]
        elif isinstance(value,dict):
            attributes += [getAttributesFrom(v) for v in value.values()]
        else:
            attributes += getAttributesFrom(value)
    
    attributes =  list(filter(lambda a: a != [], attributes))    
    
    return list(set(_flattenList(attributes)))

