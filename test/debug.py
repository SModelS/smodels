"""
.. module:: debug
   :synopsis: just a simple method to help with debugging

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def printTo (  *args  ):
    with open ( "check.txt", "at" ) as f:
        # f.write ( text + "\n" )
        line = " ".join ( map ( str,  args ) )
        print ( line )
        f.write ( line + "\n" )
        f.close()
