#!/usr/bin/env python

"""
.. module:: RCFile
    :synopsis: missing

.. moduleauthor:: missing <email@example.com>

"""


import os

def yesno ( B ):
  if B: return "yes"
  return "no"

rcfile=os.path.expanduser("~")+"/.smodelsrc"

print 'from rC=',rcfile
exists=os.path.exists ( rcfile )
#print "[RCFile.py] Check to see if %s exists: %s" % ( rcfile, yesno(exists))
if exists:
  execfile ( rcfile )
