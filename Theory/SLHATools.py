#!/usr/bin/env python

"""
.. module:: SLHATools
    :synopsis: A collection of tools needed for use and manipulation of SLHA files
    
.. moduleauthor:: Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>
    
"""

""" A collection of tools needed for use and manipulation of SLHA files """

import pyslha2
import tempfile
import os

def createSLHAFile(topo, masses = None, filename = None, branching = None, totalwidth = None):
   """ Creates an SLHA File for a certain Tx name and certain masses.

     :param topo: Tx name
     :type topo: str
     :param masses: a dictionary {pid: mass} where pid is an integer and mass is \
       an integer or float, all masses not included in the dictionary are set to 100000\
       (by default masses from initial SLHA file are taken). 
     :param filename: by default a unique random filename will be generated.
     :type filenmame: str
     :param branching: a dictionary from a dictionary {pidmom: {"piddaugther1,piddaugther2,...": branching ratio, ...}...}\
       pidmom is an integer, pids of daugther particles are given as a string seperated by ',',\
       the branching ratio is a float or integer.
     :param totalwidth: a dictionary {pid: total width} where pid is an integer and total width is \
       an integer or float (by default total width from initial SLHA file is taken).
     :returns: the filename in string format
     """

   slha = pyslha2.readSLHAFile('../slha/%s.slha' %topo)

   if masses:
      for pid in masses:
         slha[0]['MASS'].entries[pid] = masses[pid]

      for pid in slha[0]['MASS'].entries:
         if not masses.has_key(pid):
            slha[0]['MASS'].entries[pid] = 1.00000000e+05

   if branching:
      for pid in slha[1]:
         for k in range(len(slha[1][pid].decays)):
#            print 'deleting', slha[1][pid].decays[k]
            del slha[1][pid].decays[k]
      for pid in branching:
         if not slha[1].has_key(pid):
            slha[1][pid] = pyslha2.Particle(pid)
            print "[SLHATools.py] Created new decay object for pid %d" %pid
#         print 'number of decays:', len(slha[1][pid].decays)
         for decay in branching[pid]:
            ids = decay.split(',')
            for i in ids:
               i.replace(' ','')
            slha[1][pid].add_decay(branching[pid][decay], len(ids), ids)

   if totalwidth:
      for pid in totalwidth:
         if not slha[1].has_key(pid):
            slha[1][pid] = pyslha2.Particle(pid)
         slha[1][pid].totalwidth = totalwidth[pid]

   if filename:
      pyslha2.writeSLHAFile(filename, slha)
      return filename
   else:
      filename = tempfile.mkstemp()
      pyslha2.writeSLHAFile(filename[1], slha)
      return filename[1]
