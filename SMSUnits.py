from unum import Unum

""" this introduces physical units (GeV,fb) to the framework """

# description see
# http://home.scarlet.be/be052320/Unum.html
# http://home.scarlet.be/be052320/docs.html

# install it with:
# sudo pip install ez_setup
# sudo easy_install Unum
# or sth similar

# can be switched off with a single switch
useUnits=True

Unum.VALUE_FORMAT = "%0.2E"
Unum.UNIT_HIDE_EMPTY = True

fb=Unum.unit('fb')
pb=Unum.unit('pb', 1000 * fb)

eV=Unum.unit('eV')
keV=Unum.unit('keV',10**3*eV)
MeV=Unum.unit('MeV',10**6*eV)
GeV=Unum.unit('GeV',10**9*eV)
TeV=Unum.unit('TeV',10**12*eV)

def addunit ( value, unitstring ):
  """ a function that can add units to values, but also
      makes it easy to turn this functionality off, in case "units" isnt installed
      """
  if value==None: return value
  if not useUnits: return value
  if useUnits: 
    import types as t
    ## for convenience, we add units also to tuples, lists, and dictionaries
    if type(value)==t.ListType:
      return [ addunit(x,unitstring) for x in value ]
    if type(value)==t.TupleType:
      return tuple ( [ addunit(x,unitstring) for x in value ] )
    if type(value)==t.DictType:
      ret={}
      for (k,v) in value.items():
        ret[k]=addunit(v,unitstring)
      return ret
    if type(value)!=t.FloatType and type(value)!=t.IntType:
      return value
    if unitstring=="GeV":
      return value * GeV
    if unitstring=="TeV":
      return value * TeV
    if unitstring=="fb":
      return value * fb
    if unitstring=="pb":
      return value * pb
    if unitstring=="fb-1":
      return value / fb
    print "[SMSUnits.py] Warning: dont know what to do with unit",unitstring
  return value

def rmvunit ( value, unitstring ):
  """ a function that can remove units from values, but also
      makes it easy to turn this functionality off, in case "units" isnt installed
      """
  if not useUnits: return value
  if useUnits: 
    import types as t
    if type(value) != type(1.*GeV):
      return value
    if unitstring=="GeV":
      return value.asNumber(GeV)
    if unitstring=="TeV":
      return value.asNumber(TeV)
    if unitstring=="fb":
      return value.asNumber(fb)
    if unitstring=="pb":
      return value.asNumber(pb)
    if unitstring=="fb-1":
      return value.asNumber(1/fb)

    print "[SMSUnits.py] Warning: dont know what to do with unit",unitstring
    return value

