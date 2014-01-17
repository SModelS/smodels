#!/usr/bin/env python

"""
.. module:: ToolBox
    :synopsis: A singleton class that keeps track of all external tools

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

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

  def compile ( self ):
    for (name,instance) in self.tools.items():
      install_ok=instance.checkInstallation()
      if install_ok: continue
      print "[ToolBox] installation of",name,"not correct. Trying to compile."
      ok=instance.compile()

  def get ( self, tool ):
    """ get instance of tool from the toolbox """
    if not self.tools.has_key ( tool ):
      print "[ToolBox] error: asking for non-existent tool ``%s''" % tool
      return None
    return self.tools[tool]

from ExternalTools import ExternalPythia, ExternalNllFast7

toolBox=ToolBox()

toolBox.add ( ExternalPythia( ) )
toolBox.add ( ExternalNllFast7( ) )

if __name__ == "__main__":
  import argparse, types
  argparser = argparse.ArgumentParser(description='simple script to check if the tools are installed and compiled')
  argparser.add_argument('-c','--colors',help='colorful output', action='store_true')
  argparser.add_argument('-m','--make',help='compile packages if needed', action='store_true')
  args=argparser.parse_args()
  toolBox.compile()
  toolBox.checkInstallation(colors=args.colors )
