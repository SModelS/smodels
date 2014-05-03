#!/usr/bin/python

"""                                                                                   
.. module:: setPath
   :synopsis: Sets the path such that e.g. "from tools import rootTools" works
   correctly. Called as a script, the path is printed.
                                                                                      
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>             
                                                                                      
"""

from __future__ import print_function
import sys
import inspect
import os

def configure():
    """
    Gets the path name of this file, remove set_path.py, remove the last
    subdir, the remaining string should be the base path name.
    
    """
    base = os.path.dirname(os.path.realpath(inspect.getabsfile(configure)))
    pos = base.rfind("/")
    base = base[:pos+1]
    sys.path.append(base)
    return base


configure()


if __name__ == "__main__":
    """
    Called as a script, we simply print out the path.
    
    """
    print("The following string is appended to the path variable:",
          configure())
