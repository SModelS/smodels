#!/usr/bin/env python3

from smodels.tools.databaseClient import DatabaseClient
from smodels.experiment.databaseObj import Database

class ProxyDBCreater:
    def __init__ ( self, inputfile ):
        self.inputfile = inputfile
        self.database = Database ( self.inputfile )

    def create ( self, servername, serverport ):
        self.servername = servername
        self.serverport = serverport
        client = DatabaseClient ( servername, serverport )
        for e,expRes in enumerate(self.database.expResultList):
            for d,dataset in enumerate(expRes.datasets):
                for t,txn in enumerate(dataset.txnameList):
                    self.database.expResultList[e].datasets[d].txnameList[t].dbClient = client
                    del self.database.expResultList[e].datasets[d].txnameList[t].txnameData.tri
                    if txn.txnameDataExp != None:
                        del self.database.expResultList[e].datasets[d].txnameList[t].txnameDataExp.tri

    def pprint ( self, *args ):
        print ( "[ProxyDBCreater]", " ".join(map(str,args)) )

    def store ( self, outputfile ):
        """ store the outputfile """
        self.pprint ( "writing to %s" % outputfile )
        self.database.createBinaryFile ( outputfile )

    def run ( self, really ):
        """ now run the server
        :param really: if False, then only write out command
        """
        servercmd = "../smodels/tools/databaseServer.py -p %d -d %s" % \
                      ( self.serverport, self.inputfile ) 
        if really:
            self.pprint ( "starting a server on %s as:" % self.servername )
            import subprocess
            subprocess.getoutput ( servercmd )
        else:
            print ( "not started a server. you can start one yourself:" )
            self.pprint ( servercmd )
                                  

def main ( args ): ## needed for smodelsTools
    creater = ProxyDBCreater ( args.inputfile )
    creater.create( args.servername, args.serverport )
    creater.store ( args.outputfile )
    creater.run ( args.run )
