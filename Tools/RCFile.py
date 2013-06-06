
import os

def yesno ( B ):
  if B: return "yes"
  return "no"

rcfile=os.path.expanduser("~")+"/.smodelsrc"

exists=os.path.exists ( rcfile )
print "[RCFile.py] Check to see if %s exists: %s" % ( rcfile, yesno(exists))
if exists:
  execfile ( rcfile )
