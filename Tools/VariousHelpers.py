""" various helper classes """

def getMaxLum(List):
  """ Goes through all analyses in the list and returns the maximum luminosity.
  If one of the luminosities has not been define, return None. """
  from PhysicsUnits import addunit
  maxlum = addunit(0.,'fb-1')
  for Ana in List:
    if type(Ana.lum) == type(1) or type(Ana.lum) == type(1.): Ana.lum = addunit(Ana.lum,'fb-1')
    if Ana.lum.asNumber():
      maxlum = max(maxlum,Ana.lum)
    else:
      return None
  return maxlum

def getInstallationBase():
  """ get the directory name of the installation of SMSDecomposition """
  import os
  basedir=os.getcwd()
  for i in [ "bin", "test", "regression", "Tools", "Theory", "Experiment" ]:
    n=len(i)
    if basedir[-n:]==i: basedir=basedir[:-n]
  return basedir
