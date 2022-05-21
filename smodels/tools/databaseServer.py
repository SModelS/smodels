#!/usr/bin/env python3

from smodels.experiment.databaseObj import Database
from smodels.tools.physicsUnits import GeV
from smodels.tools.caching import Cache
import socket, atexit, time, os, sys, copy
import unum

unum.Unum.VALUE_FORMAT = "%0.16E"

servers = [ [] ]

def shutdownAll ():
    if len(servers[0])==0:
        return
    print ( "[databaseServer] shutting down %d servers" % len(servers[0]) )
    for i in servers[0]:
        i.shutdown( fromwhere = "atexit" )
    servers[0] = []

class DatabaseServer:
    def __init__ ( self, dbpath, servername = None, port = None, verbose = "info",
                   rundir = "./", logfile = "@@rundir@@/dbserver.log" ):
        self.rundir = rundir
        logfile = logfile.replace("@@rundir@@",rundir)
        verbose = verbose.lower()
        verbs = { "err": 10, "warn": 20, "info": 30, "debug": 40 }
        self.verbose = 50
        self.logfile = logfile
        for k,v in verbs.items():
            if k in verbose:
                self.verbose = v
        self.setStatus ( "ramping" )
        self.pprint ( "starting a server at %s" % time.asctime() )
        if port == None:
            port = 31770
            while self.is_port_in_use ( port ):
                port += 1
            self.pprint ( "using first free port %d" % port )
        if servername == None:
            servername = socket.gethostname()
            self.pprint ( "determined servername as '%s'" % servername )
        self.servername = servername
        self.dbpath = dbpath
        self.t0 = time.time()
        self.port = port
        self.packetlength = 256
        self.nlookups = 0
        servers[0].append ( self )
        if "SLURM_JOB_ID" in os.environ:
            self.pprint ( "slurm job id is %s" % os.environ["SLURM_JOB_ID"] )
        self.setStatus ( "ramped" )

    def setStatus ( self, status ):
        """ servers have a status file that tells us if they are running """
        statusfile = self.rundir + "/serverstatus.log"
        if status == "down" and os.path.exists ( statusfile ):
            os.unlink ( statusfile )
            return
        with open ( statusfile, "wt" ) as f:
            f.write ( status + "\n" )
            f.close()

    def is_port_in_use(self, port):
        """ check if port is in use """
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            return s.connect_ex(('localhost', port)) == 0

    def run ( self, nonblocking = False ):
        """ run server
        :param nonblock: run in nonblocking mode (not yet implemented)
        """ 
        if nonblocking:
            pid = os.fork()
            if pid == 0:
                return
        self.initialize()

    def logServerStats ( self ):
        """ log our stats upon exit """
        self.pprint ( "server stats" )
        self.pprint ( "============" )
        self.pprint ( "total number of lookups: %d" % self.nlookups )
        
    def shutdown ( self, fromwhere = "unknown" ):
        self.pprint ( f"Received shutdown request from {fromwhere}" )
        self.logServerStats ()
        self.setStatus ( "down" )
        self.finish()
        if hasattr ( self, "connection" ):
            self.connection.close()
            del self.connection
        if hasattr ( self, "socket" ):
            self.socket.close()
            del self.socket
        if fromwhere != "atexit":
            try:  ## remove from list of servers
                servers[0].remove ( self )
            except:
                print ( "[databaseServer] couldnt remove %s" % str(self ) )
        sys.exit()

    def parseData ( self, data ):
        """ parse the data packet """
        data=data[2:-1]
        self.log ( 'received "%s"' % data )
        if data.startswith ( "shutdown" ):
            self.shutdown( fromwhere = "client" )
            return
        if not data.startswith ( "query " ):
            self.pprint ( "I dont understand the data packet %s" % data )
            return
        data=data[6:] ## remove the query statement
        ret = self.lookUpResult ( data )
        self.log ( 'sending result of "%s" back to the client' % ret )
        ret = (str(ret)+" "*32)[:32]
        self.connection.sendall ( bytes(ret,"utf-8") )

    def lookUpResult ( self, data ):
        self.nlookups += 1
        if self.nlookups % 20000 == 0:
            self.pprint ( "looked up %dth result" % self.nlookups )
        tokens = data.split(":")
        anaId = tokens[1]
        dType = ":". join ( tokens[2:-2] )

        txname = tokens[-2]
        massv = eval(tokens[-1])
        massvunits = copy.deepcopy ( massv )
        for ibr,br in enumerate(massv):
            for iel,el in enumerate(br):
                    massvunits[ibr][iel]=el*GeV
        expected = False 
        if tokens[0] == "exp":
            expected = True
        self.log ( 'looking up for %s,%s,%s,%s' % \
                      ( anaId, dType, txname, massv ) )
        for exp in self.expResults:
            if not exp.globalInfo.id == anaId:
                continue
            for ds in exp.datasets:
                if dType == "ul" and ds.getType() != "upperLimit":
                    continue
                if dType != "ul" and dType != ds.getID():
                    continue
                if dType != "ul" and ds.getType() != "efficiencyMap":
                    continue
                for txn in ds.txnameList:
                    if txn.txName != txname:
                        continue
                    coords = txn.txnameData.dataToCoordinates ( massv, txn.txnameData._V,
                             txn.txnameData.delta_x ) 
                    res = None
                    if expected:
                        if txn.txnameDataExp != None:
                            res = txn.txnameDataExp.getValueForPoint ( coords )
                    else:
                        res = txn.txnameData.getValueForPoint ( coords )
                    # print ( "now query", massv, anaId, ds.getType(), txname, ":", res )
                    return str(res)
        return "None"

    def finish ( self ):
        if hasattr ( self, "connection" ):
            self.connection.close()

    def listen ( self ):
        try:
            self.log ( 'connection from', self.client_address )

            # Receive the data in small chunks and retransmit it
            while True:
                data = self.connection.recv( self.packetlength )
                if data:
                    self.parseData ( str(data) )
                else:
                    self.log ( 'no more data from %s:%s' % \
                               ( self.client_address[0], self.client_address[1] ) )
                    break
        finally:
            # Clean up the connection
            self.finish()

    def log ( self, *args ):
        if self.verbose > 35:
            print ( "[databaseServer]", " ".join(map(str,args)) )
            with open ( self.logfile, "at" ) as f:
                f.write ( "[databaseServer-%s] %s\n" % \
                          ( time.strftime("%H:%M:%S"), " ".join(map(str,args)) ) )
                f.close()

    def pprint ( self, *args ):
        if self.verbose > 25:
            print ( "[databaseServer]", " ".join(map(str,args)) )
            with open ( self.logfile, "at" ) as f:
                f.write ( "[databaseServer-%s] %s\n" % \
                          ( time.strftime("%H:%M:%S"), " ".join(map(str,args)) ) )
                f.close()

    def initialize( self ):
        self.setStatus ( "initialized" )
        Cache.n_stored = 20000 ## crank up the caching
        self.db = Database ( self.dbpath )
        self.expResults = self.db.expResultList
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.server_address = ( self.servername, self.port )
        self.pprint ( 'starting up on %s port %s' % self.server_address )
        self.pprint ( 'I will be serving database %s at %s' % \
                      (self.db.databaseVersion, self.dbpath ) )
        try:
            self.sock.bind( self.server_address )
        except OSError as e:
            self.pprint ( "exception %s. is host ''%s'' reachable?" % \
                          ( e, self.server_address ) )
            sys.exit(-1)
        # Listen for incoming connections
        self.sock.listen(1)

        atexit.register ( shutdownAll )

        while True:
            # Wait for a connection
            self.setStatus ( "waiting" )
            self.log ( 'waiting for a connection' )
            self.connection, self.client_address = self.sock.accept()
            self.listen()

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='an instance of a database server' )
    argparser.add_argument ( '-d', '--dbpath',
            help='The database path [./original.pcl]',
            type=str, default="./original.pcl" )
    argparser.add_argument ( '-R', '--rundir',
            help='The rundir [./]',
            type=str, default="./" )
    argparser.add_argument ( '-p', '--port',
            help='port to listen to [31770]',
            type=int, default=None )
    argparser.add_argument ( '-v', '--verbose',
            help='verbosity [info]',
            type=str, default="info" )
    argparser.add_argument ( '-s', '--servername',
            help='server name, if not specified then determined from socket [None]',
            type=str, default=None )
    args = argparser.parse_args()

    server = DatabaseServer ( args.dbpath, args.servername, args.port, args.verbose,
                              args.rundir )
    server.run()
