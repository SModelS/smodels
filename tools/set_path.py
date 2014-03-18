#!/usr/bin/python

"""                                                                                   
.. module:: set_path
    :synopsis: sets the path such that e.g. "from tools import rootTools"
               works correctly. Called as a script, the path is printed.
                                                                                      
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>             
                                                                                      
"""

def configure():
    """ get the path name of this file, remove set_path.py, 
        remove the last subdir, the remaining string should be the
        base path name """
    import sys, inspect, os
    base=os.path.dirname ( os.path.realpath ( inspect.getabsfile(configure) ) )
    pos=base.rfind("/")
    base=base[:pos+1]
    sys.path.append ( base )
    return base

configure()

if __name__ == "__main__":
    """ called as a script, we simply print out the path """
    print "The following string is appended to the path variable:",configure()
