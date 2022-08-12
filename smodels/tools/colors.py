"""
.. module:: colors
   :synopsis: Simple facility to handle console colors

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

if __name__ == "__main__":
    for i in dir ( colors ):
        if "__" in i: continue
        if i in [ "on", "reset" ]: continue
        print ( getattr ( colors, i ) + i + colors.reset )
