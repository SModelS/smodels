""" totally trivial stuff that is nice for testing """

green='\033[0;32m'
red='\033[0;31m'
reset='\033[;0m'

def ok ( A,B,verbose=True ):
  if A==B: return "%sok.%s" % ( green, reset )
  if verbose:
    return "%sfailed. [%s]%s" % ( red, B, reset )
  return "%sfailed. %s" % ( red, reset )
