## simple code snippet to set the path
import sys, os
try:
    import smodels ## already works; nothing needs to be done
except:
    sys.path.append ( os.path.join ( *([".."]*4) ) )
