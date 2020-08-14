#!/usr/bin/env python3

import socket, sys
from smodels.tools.physicsUnits import fb, pb

class DatabaseClient:
    def __init__ ( self, servername = None, port = None, verbose = "info" ):
        self.verbose = verbose
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
        return self.send( msg )

    def send ( self, message, amount_expected=32 ):
        """ send the message.
        :param amount_expected: how many return bytes do you expect
        """
        try:
            message = bytes ( message, "UTF-8" ) 
            # Send data
            # msg = b'query obs:ATLAS-SUSY-2017-01:SRHad-Low:TChiWH:[[500,100],[500,100]]'
            self.pprint ( 'sending "%s"' % message )
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
            self.pprint ( 'received "%s", %s' % ( data, type(data) ) )
            return data

        finally:
            self.pprint ( 'closing socket' )
            self.sock.close()
            del self.sock


    def pprint ( self, *args ):
        if self.verbose == "info":
            print ( "[databaseClient]", " ".join(map(str,args)) )

    def initialize( self ):
        if hasattr ( self, "sock" ):
            return ## already initialized
        # Create a TCP/IP socket
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        # Connect the socket to the port where the server is listening
        self.server_address = ( self.servername, self.port )
        self.pprint ( 'connecting to %s port %s' % self.server_address )
        self.sock.connect( self.server_address )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='an instance of a database client' )
    argparser.add_argument ( '-s', '--servername',
            help='The server name [None]',
            type=str, default=None )
    argparser.add_argument ( '-x', '--shutdown',
            help='Send shutdown command to server',
            action = "store_true" )
    argparser.add_argument ( '-q', '--query',
            help='query message [obs:ATLAS-SUSY-2017-01:SRHad-Low:TChiWH:[[500,100],[500,100]]]',
            type=str, default="obs:ATLAS-SUSY-2017-01:SRHad-Low:TChiWH:[[500,100],[500,100]]" )
    argparser.add_argument ( '-p', '--port',
            help='port to listen to [31770]',
            type=int, default=None )
    args = argparser.parse_args()
    client = DatabaseClient ( args.servername, args.port )
    client.initialize()
    if args.shutdown:
        client.send_shutdown()
    else:
        client.query ( args.query )
