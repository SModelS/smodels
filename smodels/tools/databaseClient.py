#!/usr/bin/env python3

import socket, sys, time, random, os
from smodels.tools.physicsUnits import fb, pb
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError

class DatabaseClient:
    def __init__ ( self, servername = None, port = None, verbose = "info", 
                   logfile = "dbclient.log" ):
        verbose = verbose.lower()
        verbs = { "err": 10, "warn": 20, "info": 30, "debug": 40 }
        self.verbose = 50
        self.logfile = logfile
        for k,v in verbs.items():
            if k in verbose:
                self.verbose = v
        self.logfile = logfile
        self.packetlength = 256
        if port == None:
            port = 31770
            self.pprint ( "using default port %d" % port )
        self.port = port
        if servername == None:
            servername = socket.gethostname()
            self.pprint ( "determined servername as '%s'" % servername )
        self.servername = servername

    def send_shutdown ( self ):
        """ send shutdown request to server """
        self.initialize()
        self.send ( "shutdown", amount_expected = 0 )

    def query ( self, msg ):
        """ query a certain result, msg is eg.
            obs:ATLAS-SUSY-2016-07:ul:T1:[[5.5000E+02,4.5000E+02],[5.5000E+02,4.5000E+02]]
        """
        self.initialize()
        if not msg.startswith ( "query " ):
            msg = "query " + msg
        ret = self.send ( msg )
        return ret

    def send ( self, message, amount_expected=32 ):
        """ send the message.
        :param amount_expected: how many return bytes do you expect
        """
        try:
            message = bytes ( message, "UTF-8" ) 
            # Send data
            # msg = b'query obs:ATLAS-SUSY-2017-01:SRHad-Low:TChiWH:[[500,100],[500,100]]'
            self.pprint ( 'sending "%s"' % message )
            self.ntries = 0
            while self.ntries < 10:
                try:
                    self.sock.sendall(message)

                    # Look for the response
                    amount_received = 0
                    
                    self.pprint ( 'sent message' )
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
                    self.ntries += 1
                    dt = random.uniform ( 1, 5 ) + 5*self.ntries
                    self.pprint ( 'could not connect to %s. trying again in %d seconds' % \
                                  ( self.nameAndPort(), dt ) )
                    time.sleep ( dt )
            self.pprint ( f"could not connect after trying {self.ntries} times. aborting" )
            raise SModelSError ( "Could not connect to database" )
            
        finally:
            self.pprint ( 'closing socket' )
            self.sock.close()
            del self.sock

    def log ( self, *args ):
        if type(self.verbose)==str or self.verbose > 35:
            print ( "[databaseClient]", " ".join(map(str,args)) )
        with open ( self.logfile, "at" ) as f:
            f.write ( "[databaseClient] %s\n" % " ".join(map(str,args)) )
            f.close()

    def pprint ( self, *args ):
        if type(self.verbose)==str or self.verbose > 25:
            print ( "[databaseClient]", " ".join(map(str,args)) )
        with open ( self.logfile, "at" ) as f:
            f.write ( "[databaseServer] %s\n" % " ".join(map(str,args)) )
            f.close()


    def nameAndPort ( self ):
        return "%s:%d" % ( self.servername, self.port )

    def initialize( self ):
        if hasattr ( self, "sock" ):
            return ## already initialized
        # Create a TCP/IP socket
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout ( 120 )

        # Connect the socket to the port where the server is listening
        self.server_address = ( self.servername, self.port )
        self.pprint ( 'connecting to %s port %s' % self.server_address )
        self.ntries = 0
        while self.ntries < 20:
            try:
                self.sock.connect( self.server_address )
                return
            except (ConnectionRefusedError,ConnectionResetError,BrokenPipeError,ConnectionAbortedError) as e:
                self.ntries += 1
                dt = random.uniform ( 1, 5 ) + 10*self.ntries
                self.pprint ( 'could not connect to %s. trying again in %d seconds' % \
                              ( self.nameAndPort(), dt ) )
                time.sleep ( dt )
                self.pprint ( f'could not connect after trying {self.ntries} times. trying again' )
        self.pprint ( f'could not connect after trying {self.ntries} times. aborting' )
        raise SModelSError ( "Could not connect to database" )

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
        for i in range(50):
            pid = os.fork()
            if pid != 0:
                client = DatabaseClient ( args.servername, args.port, args.verbosity )
                for i in range(10000):
                    msg = "obs:ATLAS-SUSY-2017-01:SRHad-Low:TChiWH:[[500,100],[500,100]]"
                    print ( "client #%d" % pid )
                    client.query ( msg )
        sys.exit()
    client = DatabaseClient ( args.servername, args.port, args.verbosity, "./dbclient.log" )
    client.initialize()
    if args.shutdown:
        client.send_shutdown()
        sys.exit()
    client.query ( args.query )
