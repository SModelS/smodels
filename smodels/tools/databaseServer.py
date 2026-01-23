#!/usr/bin/env python3

from smodels.experiment.databaseObj import Database
from smodels.base.physicsUnits import GeV
import socket, atexit, time, os, sys, copy
from smodels.statistics.basicStats import observed, apriori, aposteriori, NllEvalType
import unum
from typing import Union

unum.Unum.VALUE_FORMAT = "%0.16E"

servers = [ [] ]

def shutdownAll ():
    if len(servers[0])==0:
        return
    print ( f"[databaseServer] shutting down {len(servers[0])} servers" )
    for i in servers[0]:
        i.shutdown( fromwhere = "atexit" )
    servers[0] = []

class DatabaseServer:
    def __init__ ( self, dbpath : str, servername : Union[str,None] = None,
            port : Union[int,None] = None, verbose : str = "info",
            rundir : str = "./", logfile : str = "@@rundir@@/dbserver.log" ):
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
        self.pprint ( f"starting a server at {time.asctime()}" )
        if port == None:
            port = 31770
            while self.is_port_in_use ( port ):
                port += 1
            self.pprint ( f"using first free port {port}" )
        if servername == None:
            servername = socket.gethostname()
            self.pprint ( f"determined servername as '{servername}'" )
        self.servername = servername
        self.dbpath = dbpath
        self.t0 = time.time()
        self.port = port
        self.packetlength = 256
        self.nlookups = 0
        servers[0].append ( self )
        if "SLURM_JOB_ID" in os.environ:
            self.pprint ( f"slurm job id is {os.environ['SLURM_JOB_ID']}" )
        self.setStatus ( "ramped" )

    def setStatus ( self, status : str ):
        """ servers have a status file that tells us if they are running """
        statusfile = self.rundir + "/serverstatus.log"
        if status == "down" and os.path.exists ( statusfile ):
            os.unlink ( statusfile )
            return
        with open ( statusfile, "wt" ) as f:
            f.write ( status + "\n" )
            f.close()

    def is_port_in_use(self, port : int ) -> bool:
        """ check if port is in use """
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            return s.connect_ex(('localhost', port)) == 0

    def run ( self, nonblocking : bool = False ):
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
        self.pprint ( f"total number of lookups: {self.nlookups}" )

    def shutdown ( self, fromwhere : str = "unknown" ):
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
                print ( f"[databaseServer] couldnt remove {str(self)}" )
        sys.exit()

    def parseData ( self, data : str ):
        """ parse the data packet """
        data=data[2:-1]
        self.log ( f'received "{data}"' )
        if data.startswith ( "shutdown" ):
            self.shutdown( fromwhere = "client" )
            return
        if not data.startswith ( "query " ):
            self.pprint ( f"I dont understand the data packet {data}" )
            return
        data=data[6:] ## remove the query statement
        ret = self.lookUpResult ( data )
        self.log ( f'sending result of "{ret}" back to the client' )
        ret = (str(ret)+" "*32)[:32]
        self.connection.sendall ( bytes(ret,"utf-8") )

    def lookUpResult ( self, data : str ) -> str:
        self.nlookups += 1
        if self.nlookups % 20000 == 0:
            self.pprint ( f"looked up {self.nlookups}th result" )
        tokens = data.split(":")
        anaId = tokens[1]
        dType = ":". join ( tokens[2:-2] )

        txname = tokens[-2]
        massv = eval(tokens[-1])
        massvunits = copy.deepcopy ( massv )
        for ien,mass in enumerate(massvunits):
            massvunits[ien]=mass*GeV
        evaluationType = observed
        if tokens[0] == "exp":
            evaluationType = apriori
        self.log ( f'looking up for {anaId},{dType},{txname},{massv}' )
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
                    res = None
                    if evaluationType != observed:
                        if txn.txnameDataExp != None:
                            res = txn.txnameDataExp.getValueFor ( massv )
                    else:
                        res = txn.txnameData.getValueFor ( massv )
                    # print ( "now query", massv, anaId, ds.getType(), txname, ":", res )
                    return str(res)
        return "None"

    def finish ( self ):
        if hasattr ( self, "connection" ):
            self.connection.close()

    def pprintClientAddr ( self ):
        """ super simple convenience function """
        return f"{self.client_address[0]}:{self.client_address[1]}"

    def listen ( self ):
        try:
            self.log ( f'connection from {self.pprintClientAddr()}' )

            # Receive the data in small chunks and retransmit it
            while True:
                data = self.connection.recv( self.packetlength )
                if data:
                    self.parseData ( str(data) )
                else:
                    self.log ( f'no more data from {self.pprintClientAddr()}' )
                    break
        finally:
            # Clean up the connection
            self.finish()

    def log ( self, *args ):
        if self.verbose > 35:
            print ( "[databaseServer]", " ".join(map(str,args)) )
            with open ( self.logfile, "at" ) as f:
                asctime = time.strftime("%H:%M:%S")
                f.write(f"[databaseServer-{asctime}] {' '.join(map(str,args))}\n")
                f.close()

    def pprint ( self, *args ):
        if self.verbose > 25:
            print ( "[databaseServer]", " ".join(map(str,args)) )
            with open ( self.logfile, "at" ) as f:
                asctime = time.strftime("%H:%M:%S")
                f.write(f"[databaseServer-{asctime}] {' '.join(map(str,args))}\n")
                f.close()

    def initialize( self ):
        self.setStatus ( "initialized" )
        self.db = Database ( self.dbpath )
        self.expResults = self.db.expResultList
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.server_address = ( self.servername, self.port )
        self.pprint ( f'starting up on {self.servername} port {self.port}' )
        self.pprint ( f'I will be serving database {self.db.databaseVersion} at {self.dbpath}' )
        try:
            self.sock.bind( self.server_address )
        except OSError as e:
            self.pprint ( f"exception {e}: is host '{self.servername}:{self.port}' reachable?" )
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
