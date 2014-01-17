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

  def pathOfExecutable ( self ):
    return self.executable_path


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

  def get ( self, tool ):
    """ get instance of tool from the toolbox """
    if not self.tools.has_key ( tool ): 
      print "[ToolBox] error: asking for non-existent tool ``%s''" % tool
      return None
    return self.tools[tool]

class ExternalPythia(ExternalTool):
  def __init__ ( self, executable_path="../pythia_lhe", test_params_path="../etc/external_lhe.test", 
                 src_path="../pythia_src", verbose=False ):
    """ 
      :param executable_path: location of executable, full path (pythia_lhe)
      :param test_params_path: location of the test config file, full path (external_lhe.test)
    """ 
    self.name="pythia"
    self.executable_path=os.path.abspath(executable_path)
    self.test_params_path=None
    if test_params_path: self.test_params_path=os.path.abspath ( test_params_path )
    self.src_path=None
    self.verbose=False
    if src_path: self.src_path=os.path.abspath(src_path)

  def pathOfExecutable ( self ):
    return self.executable_path

  def run ( self, cfg_file ):
    """ run pythia
      :params cfg_file: config file used
      :returns: stdout and stderr, or error message
    """
    cfg=os.path.abspath ( cfg_file )
    if not os.path.exists( cfg ): 
      return "config file ``%s'' not found" % ( cfg )
    import commands
    cmd="%s < %s" % ( self.executable_path, cfg_file )
    Out=commands.getoutput ( cmd )
    out=Out.split("\n")
    return out

  def checkInstallation ( self ):
    """ checks if installation of tool looks ok by
        looking for executable and running it """
    if not os.path.exists( self.executable_path ): return "executable ``%s'' not found" % ( self.executable_path )
    if not os.path.exists( self.test_params_path ): return "config file ``%s'' not found" % ( self.test_params_path )
    if not os.access ( self.executable_path, os.X_OK ): return "%s is not executabloe" % self.executable
    out=self.run ( self.test_params_path )
    out.pop()
    lines={ -1: " ********* Fraction of events that fail fragmentation cuts =  0.00000 *********" }
    for (nr, line) in lines.items():
      if out[nr].find(line)==-1:
        return "Something is wrong with the setup: "+str(out)
    return True

class ExternalNllFast7(ExternalTool):
  def __init__ ( self, executable_path="../nllfast/nllfast-1.2/nllfast_7TeV", 
                 cd_path="../nllfast/nllfast-1.2/",
                 test_params="gg cteq 500 600", src_path="../nllfast/nllfast-1.2/", verbose=False ):
    """ 
      :param executable_path: location of executable, full path (pythia_lhe)
      :param test_params_path: location of the test config file, full path (external_lhe.test)
    """ 
    self.name="nllfast7"
    self.executable_path=os.path.abspath(executable_path)
    self.cd_path=cd_path
    self.test_params=test_params
    self.verbose=verbose
    self.src_path=None
    if src_path: self.src_path=os.path.abspath(src_path)

  def compile ( self ):
    """ try to compile tool """
    cmd="cd %s; make" % self.src_path
    import commands
    out=commands.getoutput ( cmd )
    return out

  def run ( self, params ):
    """ run pythia
      :params cfg_file: config file used
      :returns: stdout and stderr, or error message
    """
    import commands
    cmd="cd %s; %s %s" % ( self.cd_path, self.executable_path, params )
    Out=commands.getoutput ( cmd )
    out=Out.split("\n")
    return out

  def checkInstallation ( self ):
    """ checks if installation of tool looks ok by
        looking for executable and running it """
    if not os.path.exists( self.executable_path ): return "executable ``%s'' not found" % ( self.executable_path )
    if not os.access ( self.executable_path, os.X_OK ): return "%s is not executabloe" % self.executable
    out=self.run ( self.test_params )
    lines={ -1: "500.     600.    0.193E+00  0.450E+00  0.497E+00" }
    for (nr, line) in lines.items():
      if out[nr].find(line)==-1:
        return "Something is wrong with the setup: "+str(out)
    return True


toolBox=ToolBox()
toolBox.add ( ExternalPythia( ) )
toolBox.add ( ExternalNllFast7( ) )
