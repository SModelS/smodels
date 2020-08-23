#!/usr/bin/env python3

import socket, sys, time, random, multiprocessing, random, os
from smodels.tools.physicsUnits import fb, pb
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError

class DatabaseClient:
    def __init__ ( self, servername = None, port = None, verbose = "info", 
                   rundir = "./", logfile = "@@rundir@@/dbclient.log",
                   clientid = -1 ):
        verbose = verbose.lower()
        verbs = { "err": 10, "warn": 20, "info": 30, "debug": 40 }
        self.cache = {}
        self.clientid = clientid
        self.nqueries = 0
        self.verbose = 50
        self.rundir = rundir
        logfile = logfile.replace("@@rundir@@",rundir )
        self.logfile = logfile
        for k,v in verbs.items():
            if k in verbose:
                self.verbose = v
        self.logfile = logfile
        self.pprint ( "starting a database client at %s" % time.asctime() )
        self.packetlength = 256
        if port == None:
            port = 31770
            self.pprint ( "using default port %d" % port )
        self.port = port
        if servername == None:
            servername = socket.gethostname()
            self.pprint ( "determined servername as '%s'" % servername )
        self.servername = servername
        self.maxtries = 25 ## max numbers of trying to connect
        sstatus = self.findServerStatus()
        self.pprint ( "connecting to %s port %s, rundir is %s. server status is '%s'" % \
                      ( self.servername, self.port, self.rundir, sstatus ) )
        if "SLURM_JOB_ID" in os.environ:
            self.pprint ( "slurm job id is %s" % os.environ["SLURM_JOB_ID"] )

    def send_shutdown ( self ):
        """ send shutdown request to server """
        self.initialize()
        self.send ( "shutdown", amount_expected = 0 )

    def saveStats ( self ):
        self.pprint ( "client stats after %d queries" % self.nqueries )
        self.pprint ( "=================================" )
        self.pprint ( "number of results in cache: %d" % len(self.cache) )
        maxhits=[]
        for i in self.cache.values():
            maxhits.append ( i[1] )
        self.pprint ( "highest number of hits: %d" % max(maxhits))
        self.pprint ( "lowest number of hits: %d" % min(maxhits))

    def getWaitingTime ( self ):
        """ compute a waiting time between attempts, from self.ntries """
        return random.uniform ( 1, 5 ) + 6*self.ntries**2

    def query ( self, msg ):
        """ query a certain result, msg is eg.
            obs:ATLAS-SUSY-2016-07:ul:T1:[[5.5000E+02,4.5000E+02],[5.5000E+02,4.5000E+02]]
        """
        if msg.startswith ( "query " ):
            msg = msg.replace( "query ", "" )
        if not hasattr ( self, "cache" ):
            self.cache={}
        if not hasattr ( self, "nqueries" ):
            self.nqueries = 0
        self.nqueries+= 1
        if self.nqueries % 100000 == 0:
            self.saveStats()
        if msg in self.cache:
            self.cache[msg][1]+=1
            return self.cache[msg][0]
        self.initialize()
        ret = self.send ( "query " + msg )
        self.cache[msg] = [ ret, 1 ]
        return ret

    def clearCache ( self ):
        maxsize = 100000
        if len(self.cache)<maxsize:
            return
        minclass = 1 ## remove entries that got called that many times
        while len(self.cache)>maxsize/2:
            newcache = {}
            for msg,ret in self.cache.items():
                if ret[1]>minclass:
                    newcache[msg]=ret
            self.cache = newcache
            minclass += 1

    def send ( self, message, amount_expected=32 ):
        """ send the message.
        :param amount_expected: how many return bytes do you expect
        """
        try:
            message = bytes ( message, "UTF-8" ) 
            # Send data
            # msg = b'query obs:ATLAS-SUSY-2017-01:SRHad-Low:TChiWH:[[500,100],[500,100]]'
            self.log ( 'sending "%s"' % message )
            self.ntries = 0
            while self.ntries < self.maxtries:
                try:
                    self.sock.sendall(message)

                    # Look for the response
                    amount_received = 0
                    
                    self.log ( 'sent message' )
                    if amount_expected <= 0:
                        return
                    
                    while amount_received < amount_expected:
                        data = self.sock.recv( self.packetlength )
                        amount_received += len(data)
                    data = str(data)[2:-1]
                    data = data.replace(" [fb]","*fb")
                    data = data.replace(" [pb]","*pb")
                    data = eval(data)
                    self.log ( 'received "%s"' % ( data ) )
                    return data

                except (ConnectionRefusedError,ConnectionResetError,BrokenPipeError,ConnectionAbortedError) as e:
                    dt = self.getWaitingTime()
                    self.ntries += 1
                    self.log ( 'could not connect to %s. trying again in %d seconds' % \
                               ( self.nameAndPort(), dt ) )
                    time.sleep ( dt )
            self.pprint ( f"could not connect in send, after trying {self.ntries} times. aborting" )
            raise SModelSError ( f"Could not connect to database in send, tried {self.ntries} times" )
            
        finally:
            self.log ( 'closing socket' )
            self.sock.close()
            del self.sock

    def setDefaults ( self ):
        """ put in some defaults if data members dont exist """
        if not hasattr ( self, "clientid" ):
            self.clientid = 0
        if not hasattr ( self, "maxtries" ):
            self.maxtries = 25
        if not hasattr ( self, "rundir" ):
            self.rundir = os.getcwd()
        if not hasattr ( self, "logfile" ):
            self.logfile = self.rundir + "/dbclient.log"

    def findServerStatus ( self ):
        serverfile = self.rundir + "/serverstatus.log"
        if not os.path.exists ( serverfile ):
            return "not found"
        with open ( serverfile, "rt" ) as f:
            status = f.read()
            f.close()
        return status.strip()

    def log ( self, *args ):
        if type(self.verbose)==str or self.verbose > 35:
            self.setDefaults()
            print ( "[databaseClient%d] %s" % ( self.clientid, " ".join(map(str,args)) ) )
            with open ( self.logfile, "at" ) as f:
                f.write ( "[databaseClient%d-%s] %s\n" % \
                     ( self.clientid, time.strftime("%H:%M:%S"), " ".join(map(str,args)) ) )
                f.close()

    def pprint ( self, *args ):
        if type(self.verbose)==str or self.verbose > 25:
            self.setDefaults()
            print ( "[databaseClient%d] %s" % ( self.clientid, " ".join(map(str,args)) ) )
            with open ( self.logfile, "at" ) as f:
                f.write ( "[databaseClient%d-%s] %s\n" % \
                      ( self.clientid, time.strftime("%H:%M:%S"), " ".join(map(str,args)) ) )
                f.close()

    def nameAndPort ( self ):
        return "%s:%d" % ( self.servername, self.port )

    def initialize( self ):
        if hasattr ( self, "sock" ):
            return ## already initialized
        # Create a TCP/IP socket
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout ( 240 )

        # Connect the socket to the port where the server is listening
        self.server_address = ( self.servername, self.port )
        self.ntries = 0
        if not hasattr ( self, "maxtries" ):
            self.maxtries = 25
        while self.ntries < self.maxtries:
            try:
                self.sock.connect( self.server_address )
                return
            except (socket.timeout,OSError,ConnectionRefusedError,ConnectionResetError,BrokenPipeError,ConnectionAbortedError) as e:
                dt = self.getWaitingTime()
                self.ntries += 1
                self.log ( 'could not connect to %s after %d times. trying again in %d seconds' % \
                           ( self.nameAndPort(), self.ntries, dt ) )
                time.sleep ( dt )
        self.pprint ( f'could not connect to database in initialize, after trying {self.ntries} times. aborting' )
        raise SModelSError ( "Could not connect to database in initialize, tried %d times" % self.ntries )


def stresstest( args ):
    """ this is one process in the stress test """
    verbosity, servername, port = "error", args[0], args[1]
    nr = args[2]
    client = DatabaseClient ( servername, port, verbose = verbosity,
                              logfile="@@rundir@@/dbclient%d.log" % nr,
                              clientid = nr )
    mmother = random.uniform ( 200, 900 )
    mlsp = random.uniform ( 0, mmother )
    for i in range(1000):
        if random.uniform(0,1)>.9:
            mmother = random.uniform ( 200, 900 )
            mlsp = random.uniform ( 0, mmother )
        msg = "obs:ATLAS-SUSY-2017-01:SRHad-Low:TChiWH:[[%.2f,%.2f],[%.2f,%.2f]]" % \
               ( mmother, mlsp, mmother, mlsp )
        # print ( "client #%d" % pid )
        client.query ( msg )
    print ( "finished %d" % nr )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='an instance of a database client' )
    argparser.add_argument ( '-s', '--servername',
            help='The server name [None]',
            type=str, default=None )
    argparser.add_argument ( '-v', '--verbosity',
            help='Verbosity [info]',
            type=str, default="info" )
    argparser.add_argument ( '-R', '--rundir',
            help='run directory [./]',
            type=str, default="./" )
    argparser.add_argument ( '-x', '--shutdown',
            help='Send shutdown command to server',
            action = "store_true" )
    argparser.add_argument ( '--stresstest',
            help='stress test the server, send a lot of requests',
            action = "store_true" )
    argparser.add_argument ( '-q', '--query',
            help='query message [obs:ATLAS-SUSY-2017-01:SRHad-Low:TChiWH:[[500,100],[500,100]]]',
            type=str, default="obs:ATLAS-SUSY-2017-01:SRHad-Low:TChiWH:[[500,100],[500,100]]" )
    argparser.add_argument ( '-p', '--port',
            help='port to listen to [31770]',
            type=int, default=None )
    args = argparser.parse_args()
    if args.stresstest:
        t0 = time.time()
        nproc = 20
        margs = []
        for i in range(nproc):
            margs.append ((args.servername, args.port, i ) )
        p = multiprocessing.Pool( nproc )
        p.map ( stresstest, margs )
        dt = time.time() - t0 
        print ( "[databaseClient] stress test took %.2f seconds" % dt )
        sys.exit()
    client = DatabaseClient ( args.servername, args.port, args.verbosity, 
                              rundir=args.rundir, clientid = 0 )
    client.initialize()
    if args.shutdown:
        client.send_shutdown()
        sys.exit()
    client.query ( args.query )
    client.saveStats()
