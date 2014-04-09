#!/usr/bin/env python

"""
.. module:: externalPythia6
        :synopsis: Wrapper code for pythia6

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import setPath
from externalTool import ExternalTool
from tools import logger
import tempfile
import os, shutil

class ExternalPythia6(ExternalTool):
    def __init__(self,
                 config_file="<install>/etc/pythia.card", 
                 executable_path="<install>/tools/external/pythia6/pythia_lhe",
                 src_path="<install>/pythia6/", verbose=False):
        """ 
        :param config_file: location of the config file, full path. 
             We will make a copy of this config file, and provide tools
             to meddle with the content, so a template file can be provided 
             here. For how to meddle with the content, 
             see .replaceInCfgFile and .setParameter.
        :param executable_path: location of executable, full path (pythia_lhe)
        
        """ 
        self.name="pythia6"
        self.executable_path=self.absPath ( executable_path )
        self.src_path=self.absPath ( src_path )
        self.verbose=False
        self.tempdir=tempfile.mkdtemp()
        cfgfile=self.checkFileExists ( config_file )
        shutil.copy ( cfgfile, self.tempdir+"/temp.cfg" )

    def checkFileExists ( self, File ):
        """ check if file exists, raise an IOError if it doesnt.
            return absolute file name if it does. """
        nFile=self.absPath ( File )
        if not os.path.exists ( nFile ):
            raise IOError ( "config file %s does not exist" % nFile )
        return nFile

    def tempDirectory( self ): 
        """ simply returns the temporary directory name """
        return self.tempdir

    def unlink( self, unlinkdir=True ):
        """ remove temporary files.
            :param unlinkdir: remove temp directory completely
        """
        import os
        logger.debug ( "unlinking %s " % self.tempdir )
        for File in [ "fort.61", "fort.68" ]:
            if os.path.exists ( self.tempdir + "/" + File ):
                os.unlink ( self.tempdir + "/" + File )
        if unlinkdir:
            for File in [ "temp.cfg" ]:
                os.unlink ( self.tempdir + "/" + File )
            if os.path.exists ( self.tempdir ):
                os.rmdir ( self.tempdir )

    def replaceInCfgFile ( self, replacements = { "NEVENTS": 10000, "SQRTS":8000 } ):
        """ replace certain strings in the cfg file by other strings,
            similar to setParameter.
            this is introduced as a simple mechanism to make certain
            changes to the parameter file.
            :param replacements: a dictionary of strings and values. 
            The strings will be replaced with the values.
            The dictionary keys must be strings present in the cfg file.
        """
        f=open ( self.tempdir+"/temp.cfg" )
        lines=f.readlines()
        f.close()
        f=open ( self.tempdir+"/temp.cfg", "w" )
        for line in lines:
            for (key,value) in replacements.items():
                line=line.replace(key,str(value))
            f.write ( line )
        f.close()

    def setParameter ( self, param="MSTP(163)", value=6 ):
        """ modifies the cfg file, similar to .replaceInCfgFile.
            It will set param to value, overwriting possible old values.
        """
        f=open ( self.tempdir+"/temp.cfg" )
        lines=f.readlines()
        f.close()
        f=open ( self.tempdir+"/temp.cfg", "w" )
        for line in lines: ## copy all but lines with "param"
            if not param in line: f.write ( line )
        f.write ( "%s=%s\n" % ( param, str(value) ) )
        f.close()

    def run(self, slhafile, cfgfile=None ):
        """
        Run Pythia.
        
        :param slhafile: slhafile used
        :param cfgfile: optionally supply a new cfg file. If not supplied,
                        use the one supplied at construction time.
                        this config file will not be touched or copied; 
                        it will be taken as is.
        :returns: stdout and stderr, or error message
        
        """
        import os, commands, shutil
        slha=self.checkFileExists ( slhafile )
        cfg=self.tempdir+"/temp.cfg"
        if cfgfile!=None:
            cfg=self.absPath ( cfgfile )
            logger.info ( "running with %s" % cfg )
        shutil.copy ( slha, self.tempdir+"/fort.61" )
        cmd="cd %s ; %s < %s" % \
             ( self.tempdir, self.executable_path, cfg )
        logger.debug ( "Now running %s " % cmd )
        # cmd="%s < %s" % ( self.executable_path, cfg_file )
        Out=commands.getoutput ( cmd )
        # logger.debug ( "output: %s" % str(Out) ) 
        ## out=Out.split("\n")
        self.unlink ( unlinkdir=False )
        return Out

    def compile(self):
        """ compile pythia_lhe """
        logger.info("Trying to compile pythia:")
        cmd="cd %s; make" % self.src_path
        import commands
        out=commands.getoutput ( cmd )
        logger.info(out)
        

    def fetch(self, verbose=True):
        """ fetch and unpack tarball """
        import urllib, tarfile
        tempfile="/tmp/pythia.tar.gz"
        f=open( tempfile,"w")
        if verbose:
            logger.info("Fetching tarball ...")
        R=urllib.urlopen("http://smodels.hephy.at/externaltools/pythia/pythia.tar.gz")
        l=R.readlines()
        for line in l:
            f.write ( line )
        R.close()
        f.close()
        if verbose:
            logger.info("done.")
            logger.info("Untarring ...")
        tar=tarfile.open( tempfile )
        for item in tar:
            tar.extract ( item, self.src_path )
        if verbose:
            logger.info("done.")
        

    def checkInstallation(self):
        """ checks if installation of tool looks ok by
                looking for executable and running it """
        import os
        if not os.path.exists( self.executable_path ): 
            logger.error("executable ``%s'' not found" % (self.executable_path))
            return False
        if not os.access(self.executable_path, os.X_OK ): 
            logger.error("%s is not executabloe" % self.executable)
            return False
        import SModelS
        slhafile=SModelS.installDirectory()+"/inputFiles/slha/andrePT4.slha"
        try:
            Out=self.run ( slhafile, "<install>/etc/pythia_test.card" )
            out=Out.split("\n")

            out.pop()
            lines={ -1: " ********* Fraction of events that fail fragmentation cuts =  0.00000 *********" }
            for (nr, line) in lines.items():
                if out[nr].find(line)==-1:
                    logger.error("Something is wrong with the setup: " )
                    logger.error("expected >>>%s<<< found >>>%s<<<" % ( line, out[nr] ) )
                    return False
        except Exception,e:
            logger.error ( "Something is wrong with the setup: exception %s" % e )
        return True
    
if __name__ == "__main__":
    tool=ExternalPythia6()
    print("installed:" + str(tool.installDirectory()))
    import os
    def ok ( B ):
        if B: return "ok"
        return "error"
    td_exists=ok ( os.path.exists ( tool.tempDirectory() ) )
    print("temporary directory: %s: %s"  % ( str(tool.tempDirectory() ), td_exists ) )
    print("check:" + ok(tool.checkInstallation()))
    tool.replaceInCfgFile ({ "NEVENTS": 1, "SQRTS":8000 })
    tool.setParameter ("MSTP(163)","6" )
    import SModelS
    print("run:")
    slhafile=SModelS.installDirectory()+"/inputFiles/slha/andrePT4.slha" 
    out=tool.run ( slhafile )
    #print out
    tool.unlink()
