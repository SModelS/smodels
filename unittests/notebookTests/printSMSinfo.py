import os
import sys



def printInfo(el):

    nodes_and_indices = list(zip(el.nodes,el.nodeIndices))
    nodes_and_indices = [(str(node),inode) for node,inode in nodes_and_indices[:]]
    nodes_and_indices = sorted(nodes_and_indices,key = lambda pt: pt[1])
    edges = sorted([(str(mom),str(daughter)) 
             for mom,daughter in el.edges])
             
    nodeStr = '['
    for inode,node in enumerate(nodes_and_indices):
        if (inode != 0) and inode % 2 == 0:
            nodeStr += '\n                                          '
        if inode < len(nodes_and_indices) -1:
            nodeStr += str(node)+', '
        else:
            nodeStr += str(node)+']'
            
    print('NODES:')
    print(nodeStr)
    
    edgeStr = '['
    for iedge,edge in enumerate(edges):
        if (iedge != 0) and iedge % 2 == 0:
            edgeStr += '\n                              '
        if iedge < len(edges) -1:
            edgeStr += str(edge)+', '
        else:
            edgeStr += str(edge)+']'
            
    print('EDGES:')
    print(edgeStr)
    
