#!/usr/bin/env python

"""
.. module:: SLHATools
    :synopsis: a collection of tools needed for use and manipulation of SLHA files
    
.. moduleauthor:: someone <email@example.com>
    
"""

""" A collection of tools needed for use and manipulation of SLHA files """

import pyslha2
import tempfile
import os

def createSLHAFile(topo, masses, filename = None):
   """ Creates an SLHA File for a certain Tx name and certain masses.

     :param topo: a string with the Tx name
     :param masses: a dictionary {pid: mass} where pid is an integer and mass is \
       an integer or float (by default the mass is set to 100000.). 
     :param filename: a string, by default a unique random filename will be generated. 
     :returns: the filename in string format
     """

   slha = pyslha2.readSLHAFile('../slha/%s.slha' %topo)

   for pid in masses:
      slha[0]['MASS'].entries[pid] = masses[pid]

   for pid in slha[0]['MASS'].entries:
      if not masses.has_key(pid):
         slha[0]['MASS'].entries[pid] = 1.00000000e+05

   if filename:
      pyslha2.writeSLHAFile(filename, slha)
      return filename
   else:
      filename = tempfile.mkstemp()
      pyslha2.writeSLHAFile(filename[1], slha)
      return filename[1]
