#!/usr/bin/python

import set_path
import SModelS
from Tools import ROOTTools

print "ROOTTools.getRootVersion():",ROOTTools.getRootVersion()
print "ROOTTools.getRootVersion(true):",ROOTTools.getRootVersion(True)
print "ROOTTools.getRootPath():",ROOTTools.getRootPath()
print "ROOTTools.getRootLibraryPath():",ROOTTools.getRootLibraryPath()
print "ROOTTools.getRootPythonPath():",ROOTTools.getRootPythonPath()
print "SModelS.version()",SModelS.version()
print "SModelS.version(true)",SModelS.version(True)
