"""
.. module:: debug
   :synopsis: just a simple method to help with debugging

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def printTo (  *args  ):
    import time
    with open ( "debug.log", "at" ) as f:
        # f.write ( text + "\n" )
        line = " ".join ( map ( str,  args ) )
        print ( line )
        f.write ( f"{time.asctime()}: {line}\n" )
        f.close()
