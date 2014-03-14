""" this sets the path so we can write e.g. "from theory import blahblah """

def configure():
    """ get the path name of this file, remove set_path.py, 
        remove the last subdir, the remaining string should be the
        base path name """
    import sys, inspect, os
    base=os.path.dirname ( os.path.realpath ( inspect.getabsfile(configure) ) )
    pos=base.rfind("/")
    base=base[:pos+1]
    sys.path.append ( base )

configure()
