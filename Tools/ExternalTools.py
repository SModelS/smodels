#!/usr/bin/env python

"""
.. module:: ExternalTools
    :synopsis: Wrapper code for external tools: pythia, nllfast, ...

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os

class ExternalTool:
  """ an external tool encapsulates a tool that we call via commands.getoutput.
      An ExternalTool defines how the tool is checked for proper installation,
      and how the tool is executed. """

  def __init__ ( self, name, path ):
    """ every tool has a name """
    self.name=name
    self.path=os.path.abspath(path)


class ToolBox:
  def __init__ ( self ):
    self.tools={}

  def add ( self, instance ):
    """ add a tool by passing an instance to this method """
    self.tools[instance.name]=instance

  def installationOk ( self, ok, colors ):
    import types
    green='\033[0;32m'
    red='\033[0;31m'
    reset='\033[;0m'
    if ok==True: 
      ret="installation ok!"
      if colors: ret="%s%s%s" % ( green, ret, reset )
      return ret
    ret="problem with installation"
    if type(ok)==types.StringType: ret+=" (%s)" % ok
    if colors: ret="%s%s%s" % ( red, ret, reset )
    return ret

  def checkInstallation ( self, colors=False ):
    ret="The following tools are found in the Toolbox:\n"
    for (name,instance) in self.tools.items():
      ok=instance.checkInstallation()
      ret+= "%-12s [%10s]:  %s\n" % ( name, instance.pathOfExecutable(), self.installationOk ( ok, colors ) )
    print ret

class ExternalPythia(ExternalTool):
  def __init__ ( self, name, executable_path, testparams, src_path=None ):
    """ 
      :param testparams: <executable_path> <testparams> will be run to check installation
    """ 
    self.name=name
    self.executable_path=os.path.abspath(executable_path)
    self.src_path=None
    if src_path: self.src_path=os.path.abspath(src_path)

  def pathOfExecutable ( self ):
    return self.executable_path

  def checkInstallation ( self ):
    """ checks if installation of tool looks ok by
        looking for executable and running it """
    if not os.path.exists( self.executable_path ): return "executable ``%s'' not found" % ( self.executable_path )
    if not os.access ( self.executable_path, os.X_OK ): return "%s is not executabloe" % self.executable
    return True

toolBox=ToolBox()
toolBox.add ( ExternalPythia( "pythia", "../pythia_lhe", "../pythia_src/" ) )
