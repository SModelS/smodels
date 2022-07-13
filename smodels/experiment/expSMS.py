"""
.. module:: tree
   :synopsis: This is a class for describing Simplified Model Topologies
              used for decomposing BSM models.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.theory.particleNode import ParticleNode,InclusiveParticleNode
from smodels.tools.genericSMS import GenericSMS
from smodels.experiment.expAuxiliaryFuncs import bracketToProcessStr, maximal_matching
from itertools import product
from smodels.tools.physicsUnits import fb
from collections import OrderedDict, deque

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
    def from_string(cls, stringSMS, model, finalState=None,
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
        for dec in decays:
            momStr = dec.split('>')[0].strip()
            daughtersStr = [p.strip() for p in dec.split('>')[1].split(',')]
            decayParticles.append((momStr,daughtersStr))

        # Build a mapping (particle string > nodeIndex) for all unstable particles:
        nodesDict = {}
        maxNode = 0
        for mom,daughters in decayParticles:
            if mom == 'PV':
                nodeIndex = 0  # PV is always node zero
            elif '(' in mom and ')' in mom:  # Unstable particles should always have a unique numbering
                nodeIndex = eval(mom.split('(')[1].split(')')[0])
            else:
                continue  # Stable particles will have their nodes defined later
            nodesDict[mom] = nodeIndex
            maxNode = max(maxNode,nodeIndex)  # Store highest node index

        # Make sure particles have unique nodes:
        if len(set(list(nodesDict.values()))) != len(list(nodesDict.values())):
            raise SModelSError("Input string has non unique nodes: %s" % nodesDict)

        # Add the stable particles to nodesDict
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
        for ptcStr in nodesDict:
            # Check for inclusive node
            # (in this case, the daughters have to be included as final states of the node)
            label = ptcStr.split('(')[0]
            if label.lower() == 'inclusivenode':
                node = InclusiveParticleNode()
                daughters = []
                if ptcStr in sortedSuccessors:
                    daughters = [dStr.split('(')[0]
                                 for dStr in sortedSuccessors[ptcStr]]
                daughters= [model.getParticle(label=d) for d in daughters[:]]
                node.finalStates = sorted(daughters)
            else:
                particle = model.getParticle(label=label)
                node = ParticleNode(particle=particle)

            # Change key from node index to node object
            nodesObjDict[ptcStr] = node

        # Remove daughters from inclusiveNodes from dicts:
        for ptcStr in nodesDict:
            label = ptcStr.split('(')[0]
            if label.lower() != 'inclusivenode':
                continue
            if ptcStr not in sortedSuccessors:
                continue
            daughters = sortedSuccessors[ptcStr]
            for dStr in daughters:
                if dStr in nodesObjDict:
                    nodesObjDict.pop(dStr)
            sortedSuccessors[ptcStr] = []

        # Create SMS, add all unstable nodes following the sorted order
        # and keep track of generated indices
        newSMS = cls()
        smsDict = {}
        for momStr in sortedSuccessors:
            nodeIndex = newSMS.add_node(nodesObjDict[momStr])
            smsDict[momStr] = nodeIndex
        # Add stable nodes and keep track of SMS indices
        for ptcStr in nodesObjDict:
            if ptcStr in smsDict:
                continue  # Node has been added
            nodeIndex = newSMS.add_node(nodesObjDict[ptcStr])
            smsDict[ptcStr] = nodeIndex

        # Finally add all edges according to successorsDict:
        for momStr, daughtersStr in sortedSuccessors.items():
            if not daughtersStr:
                continue
            momIndex = smsDict[momStr]
            daughterIndices = [smsDict[d] for d in daughtersStr]
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

    def matchesTo(self, other):
        """
        Check if self matches other.

        :param other: TheorySMS or ExpSMS object to be compared to

        :return: None if objects do not match or a copy of self,
                 but with the nodes from other.
        """

        if not isinstance(other,GenericSMS):
            raise SModelSError("Can not compare ExpSMS and %s" %str(type(other)))

        mapDict = self.compareSubTrees(other)

        if mapDict is None:
            return None


        # Remove unmatched nodes (which happens in case of InclusiveNodes)
        match = {n1 : n2 for n1, n2 in mapDict.items() if n2 is not None}

        # Create nodes dict with nodes from self:
        nodesDict = {n1 : node for n1,node in zip(self.nodeIndices,self.nodes)}

        # Replace relevant nodes with nodes from other:
        for n1, n2 in match.items():
            node1 = self.indexToNode(n1) # Node from self
            if node1.isInclusive or node1.inclusiveList:
                continue  # Keep node from self
            else:
                node2 = other.indexToNode(n2)  # Node from other
                nodesDict[n1] = node2  # index from self

        # Make a new tree from other
        matchedTree = other.copy()
        # Copy tree structure (node indices, node order,....)
        # from self to other and use nodesDict to set the
        # corresponding node objects
        matchedTree.copyTreeFrom(self,nodesDict)

        # Remove missing nodes (in case there are inclusive matchings)
        for n1 in matchedTree.nodeIndices:
            if n1 not in match:
                matchedTree.remove_node(n1)

        return matchedTree

    def compareSubTrees(self, other, T1_node=None, T2_node=None):
        """
        Compare the subtrees with T1_node and T2_node as roots.
        The comparison is made according to their names, particle content and final states.
        It uses the node comparison to define semantically equivalent nodes.

        :param other: TheorySMS or ExpSMS object to be compared to self.
        :param T1_node: Node index belonging to self.T1
        :param T2_node: Node index belonging to self.T2

        :return: matcDict is None (subtrees differ) or a dictionary
                 with the mapping of the nodes daughters ({nodeINd : d2}).
        """

        if T1_node is None:
            T1_node = self.rootIndex
        if T2_node is None:
            T2_node = other.rootIndex

        # print('comparing',T1_node,T2_node,self.indexToNode(T1_node),other.indexToNode(T2_node))

        # Compare node canonical names
        cName1 = self.nodeCanonName(T1_node)
        cName2 = other.nodeCanonName(T2_node)
        if cName1 != cName2:
            return None


        # Compare nodes directly (canon name and particle content)
        node1 = self.indexToNode(T1_node)
        node2 = other.indexToNode(T2_node)
        if node1.isInclusive:
            other.getFinalStates(T2_node)  # Make sure final states are defined
        if node2.isInclusive:
            self.getFinalStates(T1_node)  # Make sure final states are defined

        cmp = node1.compareTo(node2)
        # print('  node comp=',cmp)
        if cmp != 0:
            return None

        # For inclusive nodes always return True (once nodes are equal)
        if node1.isInclusive or node2.isInclusive:
            # print('  equal nodes (inclusive)')
            return {T1_node: T2_node}

        # Check for equality of daughters
        daughters1 = self.daughterIndices(T1_node)
        daughters2 = other.daughterIndices(T2_node)
        if len(daughters1) == len(daughters2) == 0:
            # print('  equal nodes (leaf)')
            return {T1_node: T2_node}


        # Define left and right nodes in order to compute matching:
        left_nodes = daughters1[:]
        right_nodes = daughters2[:]
        edges = {}
        for d1 in left_nodes:
            for d2 in right_nodes:
                mapDict = self.compareSubTrees(other,d1, d2)
                if mapDict is not None:
                    if d1 not in edges:
                        edges[d1] = {}
                    edges[d1].update({d2: mapDict})

            # If node had no matches (was not added to the graph),
            # we already know T1_node and T2_node differs
            if d1 not in edges:
                return None

        # Remove nodes and edges for inclusiveList nodes
        # and add them to the final map
        finalMap = {}
        for d1 in daughters1[:]:
            if not self.indexToNode(d1).inclusiveList:
                continue
            if not d1 in edges:
                continue
            left_nodes.remove(d1)
            for d2 in edges[d1]:
                right_nodes.remove(d2)
            finalMap[d1] = edges.pop(d1)


        # print('Left=',left_nodes)
        # print('Right=',right_nodes)
        # print('edges=',edges)
        # Compute the maximal matching
        # (mapping where each node1 is connected to a single node2)
        mapDict = maximal_matching(left_nodes, right_nodes, edges)

        # print('Maximal matching:\n',mapDict)

        # Check if the match was successful.
        # Consider a successful match if all nodes in daughters1 were matched
        # or if all nodes in daughers2 were matched and the unmatched nodes
        # in daughters1 match an InclusiveNode.
        matched = False
        # Matched left nodes:
        left_matches = len(set(mapDict.keys()))
        right_matches = len(set(mapDict.values()))
        if left_matches != len(left_nodes):
            return None
        if right_matches != len(right_nodes):
            return None

        finalMap.update(mapDict)
        for d1, d2 in list(finalMap.items()):
            daughtersMap = edges[d1][d2]
            finalMap.update(daughtersMap)
        finalMap[T1_node] = T2_node

        # print('   returning for %i = %i' %(T1_node,T2_node),finalMap)
        return finalMap

    def copy(self):
        """
        Returns a shallow copy of self.

        : return: TheorySMS object
        """

        newSMS = ExpSMS()
        newSMS._successors.update({n: daughters[:]
                                   for n, daughters in self._successors.items()})
        newSMS._predecessors = {k: v for k, v in self._predecessors.items()}
        newSMS._nodesMapping = {n: node for n, node in self._nodesMapping.items()}
        newSMS._rootIndex = self._rootIndex
        newSMS._canonName = self._canonName

        return newSMS



