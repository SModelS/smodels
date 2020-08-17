#!/usr/bin/env python3

from smodels.tools.databaseClient import DatabaseClient
from smodels.experiment.databaseObj import Database
import socket, os, subprocess, copy

class ProxyDBCreater:
    def __init__ ( self, inputfile, verbose=True ):
        self.inputfile = inputfile
        self.verbose = verbose
        self.database = Database ( self.inputfile )

    def create ( self, servername, serverport ):
        if servername == None:
            servername = socket.gethostname()
            self.pprint ( "determined servername as '%s'" % servername )
        if serverport == None:
            serverport = 31770
        self.servername = servername
        self.serverport = serverport
        self.database.client = DatabaseClient ( servername, serverport, verbose="warn",
                                                logfile = "dbclient.log" )
        for e,expRes in enumerate(self.database.expResultList):
            for d,dataset in enumerate(expRes.datasets):
                for t,txn in enumerate(dataset.txnameList):
                    self.database.expResultList[e].datasets[d].txnameList[t].dbClient = copy.copy ( client )
                    del self.database.expResultList[e].datasets[d].txnameList[t].txnameData.tri
                    if txn.txnameDataExp != None:
                        del self.database.expResultList[e].datasets[d].txnameList[t].txnameDataExp.tri

    def pprint ( self, *args ):
        if self.verbose:
            print ( "[ProxyDBCreater]", " ".join(map(str,args)) )

    def store ( self, outputfile ):
        """ store the outputfile """
        self.outputfile = outputfile 
        self.pprint ( "writing to %s" % outputfile )
        self.database.createBinaryFile ( outputfile )

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
        if not "/" in inputfile:
            inputfile = os.getcwd() + "/" + inputfile
        servercmd = "%s/databaseServer.py -p %d -d %s -v error" % \
                      ( dirname, self.serverport, inputfile ) 
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
    creater = ProxyDBCreater ( args.inputfile )
    creater.create( args.servername, args.serverport )
    creater.store ( args.outputfile )
    creater.run ( args.run )
    if args.symlink:
        creater.symlink()
