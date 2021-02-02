#!/usr/bin/env python3

from smodels.tools.databaseClient import DatabaseClient
from smodels.experiment.databaseObj import Database
import socket, os, subprocess, copy, time #, tempfile

class ProxyDBCreater:
    def __init__ ( self, inputfile, rundir, verbose="info" ):
        self.inputfile = inputfile
        self.rundir = rundir
        self.nr = 0
        self.verbstring = verbose
        if type(verbose)==int:
            self.verbose = verbose
        else:
            verbose = verbose.lower()
            verbs = { "err": 10, "warn": 20, "info": 30, "debug": 40 }
            self.verbose = 50
            for k,v in verbs.items():
                if k in verbose:
                    self.verbose = v
        self.database = Database ( self.inputfile )

    def create ( self, servername, serverport ):
        if servername == None:
            servername = socket.gethostname()
            self.pprint ( "determined servername as '%s'" % servername )
        if serverport == None:
            serverport = 31770
        self.servername = servername
        self.serverport = serverport
        self.database.client = DatabaseClient ( servername, serverport,
                verbose=self.verbstring, rundir = self.rundir, clientid = self.nr )
        for e,expRes in enumerate(self.database.expResultList):
            for d,dataset in enumerate(expRes.datasets):
                for t,txn in enumerate(dataset.txnameList):
                    self.database.expResultList[e].datasets[d].txnameList[t].dbClient = copy.copy ( self.database.client )
                    del self.database.expResultList[e].datasets[d].txnameList[t].txnameData.tri
                    if txn.txnameDataExp != None:
                        del self.database.expResultList[e].datasets[d].txnameList[t].txnameDataExp.tri

    def pprint ( self, *args ):
        if self.verbose > 25:
            print ( "[proxyDBCreater-%s] %s" % \
                    ( time.strftime("%H:%M:%S"), " ".join(map(str,args)) ) )

    def store ( self, outputfile ):
        """ store the outputfile """
        self.outputfile = outputfile
        self.pprint ( "writing to %s" % outputfile )
        if os.path.exists ( outputfile ):
            os.unlink ( outputfile )
        ## first create it as temporary file, then move
        tempf = outputfile + ".tmp" # tempfile.mktemp ( suffix=".pcl" )
        self.database.createBinaryFile ( tempf )
        #cmd = f"mv {tempf} {outputfile}"
        #subprocess.getoutput ( cmd )
        os.rename ( tempf, outputfile ) ## would only work on same device

    def symlink ( self ):
        """ set a symlink from self.outputfile to default.pcl """
        dirname = os.path.dirname ( self.outputfile )
        symfile = f"{dirname}/default.pcl"
        self.pprint ( "setting a symlink from %s to %s" % \
                      ( self.outputfile, symfile ) )
        if os.path.exists ( symfile ):
            os.unlink ( symfile )
        cmd = f"ln -s {self.outputfile} {symfile}"
        subprocess.getoutput ( cmd )

    def run ( self, really ):
        """ now run the server
        :param really: if False, then only write out command
        """
        dirname = os.path.dirname ( __file__ )
        inputfile = self.inputfile
        #if not "/" in inputfile:
        #    inputfile = os.getcwd() + "/" + inputfile
        servercmd = "%s/databaseServer.py -R %s -p %d -d %s -v %s" % \
                      ( dirname, self.rundir, self.serverport, inputfile, self.verbstring )
        if really:
            self.pprint ( "starting a server on %s: %s" % \
                          ( self.servername, servercmd ) )
            import subprocess
            a = subprocess.getoutput ( servercmd )
            self.pprint ( "output %s" % a )
        else:
            print ( "not started a server. you can start one yourself:" )
            self.pprint ( servercmd )


def main ( args ): ## needed for smodelsTools
    creater = ProxyDBCreater ( args.inputfile, args.rundir, args.verbose )
    creater.create( args.servername, args.serverport )
    creater.store ( args.outputfile )
    if args.symlink:
        creater.symlink()
    creater.run ( args.run )

if __name__ == "__main__":
    import sys
    args = ""
    for i in sys.argv[1:]:
        if " " in i or "," in i:
            i = '"%s"' % i
        args += i + " "
    cmd = "./smodelsTools.py proxydb %s" % args 
    print ( "call %s" % cmd )
    #a = subprocess.getoutput ( cmd )
    #print ( a )
