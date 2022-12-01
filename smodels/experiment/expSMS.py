"""
.. module:: tree
   :synopsis: This is a class for describing Simplified Model Topologies
              used for describing experimental results BSM models.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.base.particleNode import ParticleNode,InclusiveParticleNode
from smodels.base.genericSMS import GenericSMS
from smodels.experiment.expAuxiliaryFuncs import maximal_matching, bracketToProcessStr
from smodels.base.physicsUnits import fb
from itertools import product
from collections import OrderedDict

class ExpSMS(GenericSMS):
    """
    A class for describing Simplified Model Topologies generated by the decompostion
    of full BSM models.
    """

    def __init__(self):
        """
        Initialize basic attributes.
        """

        GenericSMS.__init__(self)

    @classmethod
    def from_string(cls, stringSMS, model=None, finalState=None,
                     intermediateState=None):
        """
        Converts a string describing an SMS to a SMS object. It accepts
        the (old) bracket notation or the process notation. For the old notation the
        optional arguments finalState and intermediateState can also be defined.
        If the argument model is defined, the particle labels will be converted to
        Particle objects from the model. Otherwise the nodes will hold the particle strings.

        :param stringSMS: The process in string format
                              (e.g. '(PV > gluino(1),squark(2)), (gluino(1) >
                               MET,jet,jet), (squark(2) > HSCP,u)' or [[['jet','jet']],[['u']]]).
                               The particle labels should match the particles in the Model
                               (if Model != None).
        :parameter model: The model (Model object) to be used when converting particle labels to
                          particle objects. If None, the nodes will only store the particle labels.
        :parameter finalState: (optional) list containing the final state labels for each branch
                               (e.g. ['MET', 'HSCP'] or ['MET','MET'])
        :parameter intermediateState: (optional) nested list containing intermediate state labels
                                         for each branch  (e.g. [['gluino'], ['gluino']])
        """

        # First check if string is in old format:
        if '[' in stringSMS and ']' in stringSMS:
            procString = bracketToProcessStr(stringSMS, finalState=finalState,
                                               intermediateState=intermediateState)
        elif '>' in stringSMS and 'PV' in stringSMS:
            procString = stringSMS
        else:
            raise SModelSError("Could not recognize string format for element (%s)" % stringSMS)


        decays = procString.replace(" ", "").split("),(")
        decays[0] = decays[0][1:]  # Remove trailing parenthesis
        decays[-1] = decays[-1][:-1]  # Remove remaining parenthesis

        # Split decays into mother and daughters tuples:
        decayParticles = []
        allParticles = []
        for dec in decays:
            momStr = dec.split('>')[0].strip()
            daughtersStr = [p.strip() for p in dec.split('>')[1].split(',')]
            decayParticles.append((momStr,daughtersStr))
            allParticles += [momStr]+daughtersStr

        # Build a mapping (particle string > nodeIndex) for all unstable particles:
        nodesDict = {}
        for particle in allParticles:
            # If the particle has a node assigned, store it
            if '(' in particle and ')' in particle:  # Particles should always have a unique numbering
                nodeIndex = eval(particle.split('(')[1].split(')')[0])
            elif particle == 'PV':
                nodeIndex = 0  # PV is always node zero
            else:
                continue  # Stable particles will have their nodes defined later
            nodesDict[particle] = nodeIndex

        # Store highest node index defined so far:
        maxNode = max(nodesDict.values())

        # Make sure particles have unique nodes:
        if len(set(list(nodesDict.values()))) != len(list(nodesDict.values())):
            raise SModelSError("Input string has non unique nodes: %s" % nodesDict)

        # Add the particles with undefined nodes to nodesDict
        # and create successors dict
        successorsStr = {}
        for mom,daughters in decayParticles:
            successorsStr[mom] = []
            for ptc in daughters:
                if ptc in nodesDict:
                    successorsStr[mom].append(ptc)
                    continue
                maxNode = maxNode+1
                # Enumerate the unstable daughter labels, so they are unique
                ptcLabel = '%s(%i)' %(ptc,maxNode)
                nodesDict[ptcLabel] = maxNode
                successorsStr[mom].append(ptcLabel)

        # Sort successors according to the node index
        sortedMoms = sorted(successorsStr.keys(), key = lambda momStr : nodesDict[momStr])
        sortedSuccessors = OrderedDict()
        for mom in sortedMoms:
            sortedDaughters = sorted(successorsStr[mom], key = lambda dStr : nodesDict[dStr])
            sortedSuccessors[mom] = sortedDaughters

        # Convert strings to node objects:
        nodesObjDict = {}
        for ptcStr,nodeIndex in nodesDict.items():
            # Check for inclusive node
            label = ptcStr.split('(')[0]
            if label.lower() == 'inclusive' or label.lower() == 'inclusivenode':
                node = InclusiveParticleNode()
            else:
                # Check for inclusive lists:
                inclusiveList = False
                if label[0] == '*':
                    label = label[1:]
                    inclusiveList = True
                if model is None:
                    particle = label
                else:
                    particle = model.getParticle(label=label)
                node = ParticleNode(particle=particle,
                                    inclusiveList=inclusiveList)

            # Change key from node index to node object
            nodesObjDict[nodeIndex] = node

        # Add all nodes following the nodeIndex order
        newSMS = ExpSMS()
        for nodeIndex in sorted(nodesObjDict.keys()):
            node = nodesObjDict[nodeIndex]
            newSMS.add_node(node,nodeIndex=nodeIndex)

        # Finally add all edges according to successorsDict:
        for momStr, daughtersStr in sortedSuccessors.items():
            if not daughtersStr:
                continue
            momIndex = nodesDict[momStr]
            daughterIndices = [nodesDict[d] for d in daughtersStr]
            newSMS.add_edges_from(product([momIndex],daughterIndices))

        return newSMS

    def __eq__(self, other):
        """
        SMS equality based on the compareTo method.

        :parameter other: TheorySMS object

        :return: True if objects are equivalent.
        """

        match = self.matchesTo(other)

        return (match is not None)

    def __hash__(self):
        return object.__hash__(self)

    def matchesTo(self, other):
        """
        Check if self matches other.

        :param other: TheorySMS or ExpSMS object to be compared to

        :return: None if objects do not match or a copy of self,
                 but with the nodes from other.
        """

        if not isinstance(other,GenericSMS):
            raise SModelSError("Can not compare ExpSMS and %s" %str(type(other)))

        mapDict = self.computeMatchingDict(other,self.rootIndex,
                                           other.rootIndex)

        if mapDict is None:
            return None


        # Invert mapping dictionary {self -> other} -> {other -> self}
        # following the ordering of self
        invMapDict = OrderedDict()
        for n1 in self.nodeIndices:
            if not n1 in mapDict:  # For InclusiveNodes the match is partial
                continue
            n2 = mapDict[n1]
            # For inclusiveLists, set the firt match
            if isinstance(n2,dict):
                n2 = list(n2.keys())[0]
            invMapDict[n2] = n1


        # Get max node number from self
        maxNode = max(invMapDict.values())
        # Get unmatched nodes from other (in case of InclusiveNodes or InclusiveLists)
        missingNodes = set(other.nodeIndices).difference(set(invMapDict.keys()))
        # Add missing nodes to invMapDict:
        for n2 in missingNodes:
            maxNode += 1
            invMapDict[n2] = maxNode

        # Make a new tree from other
        matchedTree = other.copy()
        # Relabel nodes following the numbering and ordering of self:
        matchedTree.relabelNodeIndices(invMapDict)
        # Sort according to node order from self:
        nodesOrder = self.nodeIndices
        matchedTree.sortAccordingTo(nodesOrder)

        return matchedTree

    def computeMatchingDict(self, other, n1, n2):
        """
        Compare the subtrees with n1 and n2 as roots and
        return a dictionary with the node matchings {n1 : n2,...}.
        It uses the node comparison to define semantically equivalent nodes.

        :param other: TheorySMS or ExpSMS object to be compared to self.
        :param n1: Node index belonging to self
        :param n2: Node index belonging to other

        :return: None (subtrees differ) or a dictionary
                 with the mapping of the nodes and their daughters
                ({n1 : n2, d1 : d2, ...}).
        """

        if n1 is None:
            n1 = self.rootIndex
        if n2 is None:
            n2 = other.rootIndex

        # print('comparing',n1,n2,self.indexToNode(n1),other.indexToNode(n2))

        # Compare node canonical names
        cName1 = self.nodeCanonName(n1)
        cName2 = other.nodeCanonName(n2)
        if cName1 != cName2:
            return None

        # Get node objects
        node1 = self.indexToNode(n1)
        node2 = other.indexToNode(n2)

        # Compare nodes directly (canon name and particle content)
        if node1 != node2:
            return None

        # Check for equality of daughters
        # In the case of inclusive nodes, jump to final states:
        if node1.isInclusive or node2.isInclusive:
            daughters1 = self.getFinalStates(n1)
            daughters2 = other.getFinalStates(n2)
        else:
            daughters1 = self.daughterIndices(n1)
            daughters2 = other.daughterIndices(n2)

        # print('daughters1 = ',list(zip(daughters1,self.indexToNode(daughters1))))
        # print('daughters2 = ',list(zip(daughters2,other.indexToNode(daughters2))))

        if len(daughters1) == len(daughters2) == 0:
            return {n1: n2}


        # Define left and right nodes in order to compute matching:
        left_nodes = daughters1[:]
        right_nodes = daughters2[:]
        edges = {}
        for d1 in left_nodes:
            for d2 in right_nodes:
                mapDict = self.computeMatchingDict(other,d1, d2)
                if mapDict is not None:
                    if d1 not in edges:
                        edges[d1] = {}
                    edges[d1].update({d2: mapDict})

            # If node had no matches (was not added to the graph),
            # and it is not an inclusiveList, we already know n1 and n2 differs
            if d1 not in edges:
                if not self.indexToNode(d1).inclusiveList:
                    return None

        # Remove nodes and edges from inclusiveList nodes
        # and add them to the final map
        finalMap = {}
        for d1 in left_nodes[:]:
            if not self.indexToNode(d1).inclusiveList:
                continue
            left_nodes.remove(d1)
            if not d1 in edges:
                continue
            for d2 in edges[d1]:
                right_nodes.remove(d2)
            finalMap[d1] = edges.pop(d1)

        for d2 in right_nodes[:]:
            if not other.indexToNode(d2).inclusiveList:
                continue
            right_nodes.remove(d2)
            for d1 in list(edges.keys()):
                if d2 in edges[d1]:
                    left_nodes.remove(d1)
                    finalMap[d1] = edges.pop(d1)


        # print('Left=',left_nodes)
        # print('Right=',right_nodes)
        # print('edges=',edges)

        # If number of left and right nodes do not match,
        # there will be no full matching:
        if len(left_nodes) != len(right_nodes):
            return None
        # Compute the maximal matching
        # (mapping where each node1 is connected to a single node2)
        mapDict = maximal_matching(left_nodes, right_nodes, edges)

        # print('finalMap before:\n',finalMap)
        # print('Maximal matching:\n',mapDict)


        # Check if the match was successful.
        left_matches = len(set(mapDict.keys()))
        right_matches = len(set(mapDict.values()))
        if left_matches != len(left_nodes):
            return None
        if right_matches != len(right_nodes):
            return None

        for d1, d2 in list(mapDict.items()):
            daughtersMap = edges[d1][d2]
            finalMap.update(daughtersMap)
        finalMap[n1] = n2

        # print('   returning for %i = %i' %(n1,n2),finalMap)
        return finalMap

    def identicalTo(self,other):
        
        if self.canonName != other.canonName:
            return False
        
        for nodeIndex in self.dfsIndexIterator():
            cName1 = self.nodeCanonName(nodeIndex)
            cName2 = other.nodeCanonName(nodeIndex)
            if cName1 != cName2:
                return False
            node1 = self.indexToNode(nodeIndex)
            node2 = other.indexToNode(nodeIndex)
            if node1.particle is not node2.particle:
                return False
        return True

    def copy(self, emptyNodes=False):
        """
        Returns a shallow copy of self.

        :param emptyNodes: If True, does not copy any of the nodes from self.

        :return: TheorySMS object
        """

        newSMS = ExpSMS()
        if not emptyNodes:
            nodesObjDict = {n : node for n,node in zip(self.nodeIndices,self.nodes)}
            newSMS.copyTreeFrom(self, nodesObjDict)

        return newSMS

    def compareNodes(self,other,nodeIndex1,nodeIndex2):
        """
        Convenience function for defining how nodes are compared
        within the SMS. For ExpSMS the nodes are sorted according to their
        inclusiveness (InclsuiveNode or inclusiveList), 
        the particle content and finally their string representation.

        :param other: ExpSMS object (if other=self compare subtrees of the same SMS).
        :param nodeIndex1: Index of first node
        :param nodeIndex2: Index of second node

        :return: 1 if node1 > node2, -1 if node1 < node2, 0 if node1 == node2.
        """

        # Comparison parameters:
        node1 = self.indexToNode(nodeIndex1)
        node2 = other.indexToNode(nodeIndex2)
        
        # First sort nodes by inclusiveness:
        cmp1 = (not node1.isInclusive, not node1.inclusiveList)
        cmp2 = (not node2.isInclusive, not node2.inclusiveList)
        if cmp1 != cmp2:
            if cmp1 > cmp2:
                return 1
            else:
                return -1

        # If both are inclusive, compare particles:
        cmp = node1.compareTo(node2)
        if cmp != 0:
            return cmp
        
        # If particles are the same object, return 0
        if node1.particle is node2.particle:
            return 0
        # If they are not the same object, but node2 in node1(MultiParticle):
        # define node1 > node2
        elif node1.particle.contains(node2.particle):
            return 1
        # If they are not the same object, but node1 in node2(MultiParticle):
        # define node2 > node1
        elif node2.particle.contains(node1.particle):
            return -1

        # If everything else fails, compare strings:
        cmp1 = str(node1)
        cmp2 = str(node2)
        if cmp1 != cmp2:
            if cmp1 > cmp2:
                return 1
            else:
                return -1        
        
        # If nodes are identical, return 0
        return 0