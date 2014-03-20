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


class ExternalPythia(ExternalTool):
  def __init__ ( self, executable_path="<install>/tools/external/pythia6/pythia_lhe", 
                 test_params_path="<install>/etc/pythia6_test.cfg", 
                 src_path="<install>/pythia6/", verbose=False ):
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
    import os, commands
    cfg=os.path.abspath ( cfg_file )
    if not os.path.exists( cfg ): 
      return "config file ``%s'' not found" % ( cfg )
    cmd="%s < %s" % ( self.executable_path, cfg_file )
    Out=commands.getoutput ( cmd )
    out=Out.split("\n")
    return out


  def clean ( self ):
    """ remote fort.61 and fort.68 """
    import os
    for File in [ "fort.61", "fort.68" ]:
      if os.path.exists ( File ): os.unlink ( File )

  def compile ( self ):
     """ compile pythia_lhe """
     print "[ExternalPythia] trying to compile pythia:"
     cmd="cd %s; make" % self.src_path
     import commands
     out=commands.getoutput ( cmd )
     print out

  def fetch ( self, verbose=True ):
    """ fetch and unpack tarball """
    import urllib, tarfile
    tempfile="/tmp/pythia.tar.gz"
    f=open( tempfile,"w")
    if verbose: print "[ExternalPythia] fetching tarball ...",
    R=urllib.urlopen("http://smodels.hephy.at/externaltools/pythia/pythia.tar.gz")
    l=R.readlines()
    for line in l:
      f.write ( line )
    R.close()
    f.close()
    if verbose:
      print "done."
      print "[ExternalPythia] untarring ...",
    tar=tarfile.open( tempfile )
    for item in tar:
      tar.extract ( item, self.src_path )
    if verbose: print "done."

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
    self.clean()
    return True

class ExternalNllFast7(ExternalTool):
  def __init__ ( self, executable_path="<install>/tools/external/nllfast/nllfast-1.2/nllfast_7TeV", 
                 cd_path="<install>/tools/external/nllfast/nllfast-1.2/",
                 test_params="gg cteq 500 600", src_path="<install>/tools/external/nllfast/nllfast-1.2/", verbose=False ):
    """ 
      :param executable_path: location of executable, full path (pythia_lhe)
      :param test_params_path: location of the test config file, full path (external_lhe.test)
    """ 
    self.name="nllfast7"
    self.executable_path=self.absPath (executable_path)
    self.cd_path=self.absPath(cd_path)
    self.test_params=test_params
    self.verbose=verbose
    self.src_path=None
    if src_path: self.src_path=self.absPath(src_path)

  def compile ( self ):
    """ try to compile tool """
    print "[ExternalNllfast7] trying to compile nllfast7:"
    cmd="cd %s; make" % self.src_path
    import commands
    out=commands.getoutput ( cmd )
    print out
    return True
  
  def fetch ( self, verbose=True ):
    """ fetch and unpack tarball """
    import urllib, tarfile
    tempfile="/tmp/nllfast7.tar.gz"
    f=open( tempfile,"w")
    if verbose: print "[ExternalNllfast7] fetching tarball ...",
    R=urllib.urlopen("http://smodels.hephy.at/externaltools/nllfast7/nllfast7.tar.gz")
    l=R.readlines()
    for line in l:
      f.write ( line )
    R.close()
    f.close()
    if verbose:
      print "done."
      print "[ExternalNllfast7] untarring ...",
    tar=tarfile.open( tempfile )
    for item in tar:
      # print item,self.src_path,
      tar.extract ( item, self.src_path+ "/" )
    if verbose: print "done."


  def unlink ( self, File ): 
    """ remove File.out """
    import os
    Fname="%s/%s.out" % ( self.cd_path, File )
    if os.path.exists ( Fname ):
      os.unlink ( Fname )

  def run_ ( self, params ):
    """ run nllfast7
      :params params: parameters used (e.g. gg cteq5 .... ) FIXME could have a fancier interface
      :returns: stdout and stderr, or error message
    """
    import commands
    cmd="cd %s; %s %s" % ( self.cd_path, self.executable_path, params )
    # print "cmd=",cmd
    Out=commands.getoutput ( cmd )
    out=Out.split("\n")
    return out

  def run ( self,  process, pdf, squarkmass, gluinomass ):
    """ run nllfast7
      :params process: which process: st, sb, gg, gdcpl, sdcpl, ss, sg, tot
      :params pdf: cteq=cteq6, mstw2008 
      :params squarkmass: squarkmass, None if squark decoupled
      :params gluinomass: gluinomass, None if gluino decoupled
      :returns: stdout and stderr, or error message
    """
    if not process in [ "st", "sb", "gg", "gdcpl", "sdcpl", "ss", "sg", "tot" ]:
      return None
    if not pdf in [ "cteq", "cteq6", "mstw", "mstw2008" ]: return None
    if not squarkmass: 
      return run_("%s %s %s") % ( process, pdf, gluinomass )
    if not gluinomass: 
      return run_("%s %s %s") % ( process, pdf, squarkmass )
    return run_("%s %s %s %s") % ( process, pdf, squarkmass, gluinomass )

  def checkInstallation ( self ):
    """ checks if installation of tool looks ok by
        looking for executable and running it """
    import os
    if not os.path.exists( self.executable_path ): return "executable ``%s'' not found" % ( self.executable_path )
    if not os.access ( self.executable_path, os.X_OK ): return "%s is not executabloe" % self.executable
    out=self.run_ ( self.test_params )
    lines={ -1: "500.     600.    0.193E+00  0.450E+00  0.497E+00" }
    for (nr, line) in lines.items():
      if out[nr].find(line)==-1:
        return "Something is wrong with the setup: "+str(out)
    self.unlink ( "gg" )
    return True


class ExternalNllFast8(ExternalTool):
  def __init__ ( self, executable_path="<install>/tools/external/nllfast/nllfast-2.1/nllfast_8TeV", 
                 cd_path="<install>/tools/external/nllfast/nllfast-2.1/",
                 test_params="gg cteq 500 600", src_path="<install>/tools/external/nllfast/nllfast-2.1/", verbose=False ):
    """ 
      :param executable_path: location of executable, full path (pythia_lhe)
      :param test_params_path: location of the test config file, full path (external_lhe.test)
    """ 
    self.name="nllfast8"
    self.executable_path=self.absPath (executable_path)
    self.cd_path=self.absPath(cd_path)
    self.test_params=test_params
    self.verbose=verbose
    self.src_path=None
    if src_path: self.src_path=self.absPath(src_path)

  def compile ( self ):
    """ try to compile tool """
    print "[ExternalNllfast8] trying to compile nllfast8:"
    cmd="cd %s; make" % self.src_path
    import commands
    out=commands.getoutput ( cmd )
    print out
    return True
  
  def fetch ( self, verbose=True ):
    """ fetch and unpack tarball """
    import urllib, tarfile
    tempfile="/tmp/nllfast8.tar.gz"
    f=open( tempfile,"w")
    if verbose: print "[ExternalNllfast8] fetching tarball ...",
    R=urllib.urlopen("http://smodels.hephy.at/externaltools/nllfast8/nllfast8.tar.gz")
    l=R.readlines()
    for line in l:
      f.write ( line )
    R.close()
    f.close()
    if verbose:
      print "done."
      print "[ExternalNllfast8] untarring ...",
    tar=tarfile.open( tempfile )
    for item in tar:
      # print item,self.src_path,
      tar.extract ( item, self.src_path+ "/" )
    if verbose: print "done."


  def unlink ( self, File ): 
    """ remove File.out """
    import os
    Fname="%s/%s.out" % ( self.cd_path, File )
    if os.path.exists ( Fname ):
      os.unlink ( Fname )

  def run_ ( self, params ):
    """ run nllfast8
      :params params: parameters used (e.g. gg cteq5 500 500 ) 
      :returns: stdout and stderr, or error message
    """
    import commands
    cmd="cd %s; %s %s" % ( self.cd_path, self.executable_path, params )
    # print "cmd=",cmd
    Out=commands.getoutput ( cmd )
    out=Out.split("\n")
    return out

  def run ( self,  process, pdf, squarkmass, gluinomass ):
    """ run nllfast8
      :params process: which process: st, sb, gg, gdcpl, sdcpl, ss, sg, tot
      :params pdf: cteq=cteq6, mstw2008 
      :params squarkmass: squarkmass, None if squark decoupled
      :params gluinomass: gluinomass, None if gluino decoupled
      :returns: stdout and stderr, or error message
    """
    if not process in [ "st", "sb", "gg", "gdcpl", "sdcpl", "ss", "sg", "tot" ]:
      return None
    if not pdf in [ "cteq", "cteq6", "mstw", "mstw2008" ]: return None
    if not squarkmass: 
      return run_("%s %s %s") % ( process, pdf, gluinomass )
    if not gluinomass: 
      return run_("%s %s %s") % ( process, pdf, squarkmass )
    return run_("%s %s %s %s") % ( process, pdf, squarkmass, gluinomass )

  def checkInstallation ( self ):
    """ checks if installation of tool looks ok by
        looking for executable and running it """
    import os
    if not os.path.exists( self.executable_path ): return "executable ``%s'' not found" % ( self.executable_path )
    if not os.access ( self.executable_path, os.X_OK ): return "%s is not executabloe" % self.executable
    out=self.run_ ( self.test_params )
    lines={ -1: "500.     600.    0.406E+00  0.873E+00  0.953E+00" }
    for (nr, line) in lines.items():
      if out[nr].find(line)==-1:
        return "Something is wrong with the setup: "+str(out)
    self.unlink ( "gg" )
    return True

class ExternalNllFast13(ExternalTool):
  def __init__ ( self, executable_path="<install>/tools/external/nllfast/nllfast-3.0/nllfast_13TeV", 
                 cd_path="<install>/tools/external/nllfast/nllfast-3.0/",
                 test_params="gg cteq 500 600", src_path="<install>/tools/external/nllfast/nllfast-3.0/", verbose=False ):
    """ 
      :param executable_path: location of executable, full path (pythia_lhe)
      :param test_params_path: location of the test config file, full path (external_lhe.test)
    """ 
    self.name="nllfast13"
    self.executable_path=self.absPath (executable_path)
    self.cd_path=self.absPath(cd_path)
    self.test_params=test_params
    self.verbose=verbose
    self.src_path=None
    if src_path: self.src_path=self.absPath(src_path)

  def compile ( self ):
    """ try to compile tool """
    print "[ExternalNllfast13] trying to compile nllfast13:"
    cmd="cd %s; make" % self.src_path
    import commands
    out=commands.getoutput ( cmd )
    print out
    return True
  
  def fetch ( self, verbose=True ):
    """ fetch and unpack tarball """
    import urllib, tarfile
    tempfile="/tmp/nllfast13.tar.gz"
    f=open( tempfile,"w")
    if verbose: print "[ExternalNllfast7] fetching tarball ...",
    R=urllib.urlopen("http://pauli.uni-muenster.de/~akule_01/nllfast/nllfast-3.0.tar.gz")
    l=R.readlines()
    for line in l:
      f.write ( line )
    R.close()
    f.close()
    if verbose:
      print "done."
      print "[ExternalNllfast13] untarring ...",
    tar=tarfile.open( tempfile )
    for item in tar:
      # print item,self.src_path,
      tar.extract ( item, self.src_path+ "/" )
    if verbose: print "done."


  def unlink ( self, File ): 
    """ remove File.out """
    import os
    Fname="%s/%s.out" % ( self.cd_path, File )
    if os.path.exists ( Fname ):
      os.unlink ( Fname )

  def run_ ( self, params ):
    """ run nllfast13
      :params params: parameters used (e.g. gg cteq5 .... ) FIXME could have a fancier interface
      :returns: stdout and stderr, or error message
    """
    import commands
    cmd="cd %s; %s %s" % ( self.cd_path, self.executable_path, params )
    # print "cmd=",cmd
    Out=commands.getoutput ( cmd )
    out=Out.split("\n")
    return out

  def run ( self,  process, pdf, squarkmass, gluinomass ):
    """ run nllfast13
      :params process: which process: st, sb, gg, gdcpl, sdcpl, ss, sg, tot
      :params pdf: cteq=cteq6, mstw2008 
      :params squarkmass: squarkmass, None if squark decoupled
      :params gluinomass: gluinomass, None if gluino decoupled
      :returns: stdout and stderr, or error message
    """
    if not process in [ "st", "sb", "gg", "gdcpl", "sdcpl", "ss", "sg", "tot" ]:
      return None
    if not pdf in [ "cteq", "cteq6", "mstw", "mstw2008" ]: return None
    if not squarkmass: 
      return run_("%s %s %s") % ( process, pdf, gluinomass )
    if not gluinomass: 
      return run_("%s %s %s") % ( process, pdf, squarkmass )
    return run_("%s %s %s %s") % ( process, pdf, squarkmass, gluinomass )

  def checkInstallation ( self ):
    """ checks if installation of tool looks ok by
        looking for executable and running it """
    import os
    if not os.path.exists( self.executable_path ): return "executable ``%s'' not found" % ( self.executable_path )
    if not os.access ( self.executable_path, os.X_OK ): return "%s is not executabloe" % self.executable
    out=self.run_ ( self.test_params )
    lines={ -1: "500.     600.    0.394E+01  0.690E+01  0.731E+01" }
    for (nr, line) in lines.items():
      if out[nr].find(line)==-1:
        return "Something is wrong with the setup: "+str(out)
    self.unlink ( "gg" )
    return True


class ExternalNllFast14(ExternalTool):
  def __init__ ( self, executable_path="<install>/tools/external/nllfast/nllfast-4.01dcpl/nllfast_14TeV", 
                 cd_path="<install>/tools/external/nllfast/nllfast-4.01dcpl/",
                 test_params="gdcpl cteq 500 600", src_path="<install>/tools/external/nllfast/nllfast-4.01dcpl/", verbose=False ):
    """ 
      :param executable_path: location of executable, full path (pythia_lhe)
      :param test_params_path: location of the test config file, full path (external_lhe.test)
    """ 
    self.name="nllfast14"
    self.executable_path=self.absPath (executable_path)
    self.cd_path=self.absPath(cd_path)
    self.test_params=test_params
    self.verbose=verbose
    self.src_path=None
    if src_path: self.src_path=self.absPath(src_path)

  def compile ( self ):
    """ try to compile tool """
    print "[ExternalNllfast14] trying to compile nllfast14:"
    cmd="cd %s; make" % self.src_path
    import commands
    out=commands.getoutput ( cmd )
    print out
    return True
  
  def fetch ( self, verbose=True ):
    """ fetch and unpack tarball """
    import urllib, tarfile
    tempfile="/tmp/nllfast14.tar.gz"
    f=open( tempfile,"w")
    if verbose: print "[ExternalNllfast7] fetching tarball ...",
    R=urllib.urlopen("http://pauli.uni-muenster.de/~akule_01/nllfast/nllfast-3.0.tar.gz")
    l=R.readlines()
    for line in l:
      f.write ( line )
    R.close()
    f.close()
    if verbose:
      print "done."
      print "[ExternalNllfast14] untarring ...",
    tar=tarfile.open( tempfile )
    for item in tar:
      # print item,self.src_path,
      tar.extract ( item, self.src_path+ "/" )
    if verbose: print "done."


  def unlink ( self, File ): 
    """ remove File.out """
    import os
    Fname="%s/%s.out" % ( self.cd_path, File )
    if os.path.exists ( Fname ):
      os.unlink ( Fname )

  def run_ ( self, params ):
    """ run nllfast14
      :params params: parameters used (e.g. gg cteq5 .... ) FIXME could have a fancier interface
      :returns: stdout and stderr, or error message
    """
    import commands
    cmd="cd %s; %s %s" % ( self.cd_path, self.executable_path, params )
    # print "cmd=",cmd
    Out=commands.getoutput ( cmd )
    out=Out.split("\n")
    return out

  def run ( self,  process, pdf, squarkmass, gluinomass ):
    """ run nllfast14
      :params process: which process: st, sb, gg, gdcpl, sdcpl, ss, sg, tot
      :params pdf: cteq=cteq6, mstw2008 
      :params squarkmass: squarkmass, None if squark decoupled
      :params gluinomass: gluinomass, None if gluino decoupled
      :returns: stdout and stderr, or error message
    """
    if not process in [ "st", "sb", "gg", "gdcpl", "sdcpl", "ss", "sg", "tot" ]:
      return None
    if not pdf in [ "cteq", "cteq6", "mstw", "mstw2008" ]: return None
    if not squarkmass: 
      return run_("%s %s %s") % ( process, pdf, gluinomass )
    if not gluinomass: 
      return run_("%s %s %s") % ( process, pdf, squarkmass )
    return run_("%s %s %s %s") % ( process, pdf, squarkmass, gluinomass )

  def checkInstallation ( self ):
    """ checks if installation of tool looks ok by
        looking for executable and running it """
    import os
    if not os.path.exists( self.executable_path ): return "executable ``%s'' not found" % ( self.executable_path )
    if not os.access ( self.executable_path, os.X_OK ): return "%s is not executabloe" % self.executable
    out=self.run_ ( self.test_params )
    lines={ -1: "500.    0.235E+02  0.346E+02  0.362E+02" }
    for (nr, line) in lines.items():
      if out[nr].find(line)==-1:
        return "Something is wrong with the setup: "+str(out)
    self.unlink ( "gg" )
    return True

class ExternalNllFast33(ExternalTool):
  def __init__ ( self, executable_path="<install>/tools/external/nllfast/nllfast-5.01dcpl/nllfast_33TeV", 
                 cd_path="<install>/tools/external/nllfast/nllfast-5.01dcpl/",
                 test_params="gdcpl cteq 500 600", src_path="<install>/tools/external/nllfast/nllfast-3.0/", verbose=False ):
    """ 
      :param executable_path: location of executable, full path (pythia_lhe)
      :param test_params_path: location of the test config file, full path (external_lhe.test)
    """ 
    self.name="nllfast33"
    self.executable_path=self.absPath (executable_path)
    self.cd_path=self.absPath(cd_path)
    self.test_params=test_params
    self.verbose=verbose
    self.src_path=None
    if src_path: self.src_path=self.absPath(src_path)

  def compile ( self ):
    """ try to compile tool """
    print "[ExternalNllfast33] trying to compile nllfast33:"
    cmd="cd %s; make" % self.src_path
    import commands
    out=commands.getoutput ( cmd )
    print out
    return True
  
  def fetch ( self, verbose=True ):
    """ fetch and unpack tarball """
    import urllib, tarfile
    tempfile="/tmp/nllfast33.tar.gz"
    f=open( tempfile,"w")
    if verbose: print "[ExternalNllfast7] fetching tarball ...",
    R=urllib.urlopen("http://pauli.uni-muenster.de/~akule_01/nllfast/nllfast-3.0.tar.gz")
    l=R.readlines()
    for line in l:
      f.write ( line )
    R.close()
    f.close()
    if verbose:
      print "done."
      print "[ExternalNllfast33] untarring ...",
    tar=tarfile.open( tempfile )
    for item in tar:
      # print item,self.src_path,
      tar.extract ( item, self.src_path+ "/" )
    if verbose: print "done."


  def unlink ( self, File ): 
    """ remove File.out """
    import os
    Fname="%s/%s.out" % ( self.cd_path, File )
    if os.path.exists ( Fname ):
      os.unlink ( Fname )

  def run_ ( self, params ):
    """ run nllfast33
      :params params: parameters used (e.g. gg cteq5 .... ) FIXME could have a fancier interface
      :returns: stdout and stderr, or error message
    """
    import commands
    cmd="cd %s; %s %s" % ( self.cd_path, self.executable_path, params )
    # print "cmd=",cmd
    Out=commands.getoutput ( cmd )
    out=Out.split("\n")
    return out

  def run ( self,  process, pdf, squarkmass, gluinomass ):
    """ run nllfast33
      :params process: which process: st, sb, gg, gdcpl, sdcpl, ss, sg, tot
      :params pdf: cteq=cteq6, mstw2008 
      :params squarkmass: squarkmass, None if squark decoupled
      :params gluinomass: gluinomass, None if gluino decoupled
      :returns: stdout and stderr, or error message
    """
    if not process in [ "st", "sb", "gg", "gdcpl", "sdcpl", "ss", "sg", "tot" ]:
      return None
    if not pdf in [ "cteq", "cteq6", "mstw", "mstw2008" ]: return None
    if not squarkmass: 
      return run_("%s %s %s") % ( process, pdf, gluinomass )
    if not gluinomass: 
      return run_("%s %s %s") % ( process, pdf, squarkmass )
    return run_("%s %s %s %s") % ( process, pdf, squarkmass, gluinomass )

  def checkInstallation ( self ):
    """ checks if installation of tool looks ok by
        looking for executable and running it """
    import os
    if not os.path.exists( self.executable_path ): return "executable ``%s'' not found" % ( self.executable_path )
    if not os.access ( self.executable_path, os.X_OK ): return "%s is not executabloe" % self.executable
    out=self.run_ ( self.test_params )
    lines={ -1: "500.    0.257E+03  0.383E+03  0.393E+03" }
    for (nr, line) in lines.items():
      if out[nr].find(line)==-1:
        return "Something is wrong with the setup: "+str(out)
    self.unlink ( "gg" )
    return True

