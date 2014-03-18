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
      return True
    return False

parseRCFile()

if __name__ == "__main__":
    """ called as script, we check if there is a smodelsrc file """
    T=parseRCFile()
    print "Checking if a ~/.smodelsrc file exists:",
    if not T:
        print "no ~/.smodelsrc file found."
    else:
        print "found ~/.smodelsrc file."
