#!/usr/bin/env python

"""
.. module:: externalTool
    :synopsis: Wrapper code for external tools: base class

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

class ExternalTool:
  """ an external tool encapsulates a tool that we call via commands.getoutput.
      An ExternalTool defines how the tool is checked for proper installation,
      and how the tool is executed. """

  def installDirectory ( self ):
    T=self.executable_path
    P=T.rfind("/")
    if P==-1: return ""
    return self.executable_path[:P]

  def pathOfExecutable ( self ):
    return self.executable_path

  def basePath ( self ):
    """ get the base install path """
    import os, inspect
    base=os.path.dirname (  inspect.getabsfile( self.basePath) )
    return base
 

  def absPath ( self, path ):
    """ return the absolute path of <path>, replacing <install> with
        the install directory """
    if None: return None
    import os
    installdir=self.basePath()
    installdir=installdir.replace("tools","")
    path=path.replace("<install>",installdir)
    path=os.path.abspath( path )
    # path=path.replace("tools/","tools/external/")
    # path=path.replace("tools/","")
    return path

