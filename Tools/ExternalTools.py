#!/usr/bin/env python

"""
.. module:: ExternalTools
    :synopsis: Wrapper code for external tools: pythia, nllfast, ...

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

class ExternalTool:
  """ an external tool encapsulates a tool that we call via commands.getoutput.
      An ExternalTool defines how the tool is checked for proper installation,
      and how the tool is executed. """

  def pathOfExecutable ( self ):
    return self.executable_path

  def absPath ( self, path ):
    """ return the absolute path of <path>, replacing <install> with
        the install directory """
    if None: return None
    import os, inspect
    installdir=os.path.dirname (  inspect.getabsfile( self.absPath) )
    path=path.replace("<install>",installdir)
    path=os.path.abspath( path )
    path=path.replace("Tools/","")
    return path


class ExternalPythia(ExternalTool):
  def __init__ ( self, executable_path="<install>/pythia_lhe", 
                 test_params_path="<install>/etc/external_lhe.test", 
                 src_path="<install>/pythia_src", verbose=False ):
    """ 
      :param executable_path: location of executable, full path (pythia_lhe)
      :param test_params_path: location of the test config file, full path (external_lhe.test)
    """ 
    self.name="pythia"
    self.executable_path=self.absPath ( executable_path )
    self.test_params_path=self.absPath ( test_params_path )
    self.src_path=self.absPath ( src_path )
    self.verbose=False

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
    import os
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
    self.executable_path=self.absPath (executable_path)
    self.cd_path=cd_path
    self.test_params=test_params
    self.verbose=verbose
    self.src_path=None
    if src_path: self.src_path=self.absPath(src_path)

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
    import os
    if not os.path.exists( self.executable_path ): return "executable ``%s'' not found" % ( self.executable_path )
    if not os.access ( self.executable_path, os.X_OK ): return "%s is not executabloe" % self.executable
    out=self.run ( self.test_params )
    lines={ -1: "500.     600.    0.193E+00  0.450E+00  0.497E+00" }
    for (nr, line) in lines.items():
      if out[nr].find(line)==-1:
        return "Something is wrong with the setup: "+str(out)
    return True


