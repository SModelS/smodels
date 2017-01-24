"""
.. module:: colors
   :synopsis: Simple facility to handle console colors

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

class Colors:
    def __init__ ( self ):
        self.on = True
        setattr ( self, "error", self.red )
        setattr ( self, "warn", self.yellow )
        setattr ( self, "info", self.green )
        setattr ( self, "debug", "" )

    @property
    def magenta( self ):
        if self.on == False: return ""
        return '\033[0;35m'

    @property
    def green( self ):
        if self.on == False: return ""
        return '\033[0;32m'

    @property
    def red ( self ):
        if self.on == False: return ""
        return '\033[0;31m'

    @property
    def yellow ( self ):
        if self.on == False: return ""
        return '\033[0;33m'

    @property
    def cyan ( self ):
        if self.on == False: return ""
        return '\033[0;36m'

    @property
    def blue ( self ):
        if self.on == False: return ""
        return '\033[0;34m'

    @property
    def reset ( self ):
        if self.on == False: return ""
        return '\033[;0m'

colors = Colors()

if __name__ == "__main__":
    for i in dir ( colors ):
        if "__" in i: continue
        if i in [ "on", "reset" ]: continue
        print ( getattr ( colors, i ) + i + colors.reset )
