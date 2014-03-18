#!/usr/bin/env python

"""
.. module:: rcFile
    :synopsis: When imported, ~/.smodelsrc is parsed.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def yesno ( B ):
  if B: return "yes"
  return "no"

def parseRCFile():
    import os
    rcfile=os.path.expanduser("~")+"/.smodelsrc"
    exists=os.path.exists ( rcfile )
    #print "[RCFile.py] Check to see if %s exists: %s" % ( rcfile, yesno(exists))
    if exists:
      execfile ( rcfile )

parseRCFile()
