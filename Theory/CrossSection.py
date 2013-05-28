class CrossSection:
  """ basically a wrapper around this complicated result dictionary, to 
      make it easier to use this dictionary """

  def __init__ ( self, data ):
    self.data=data

  # make it behave much like a dictionary
  def __len__ ( self ): return len(self.data)
  def __getitem__ ( self, i ): return self.data[i]
  def items ( self ): return self.data.items()
  def __str__ ( self ): return str(self.data)

  def lhefile ( self, sqrts ):
    if sqrts==7: return self.data["lhe7file"]
    if sqrts==8: return self.data["lhe8file"]
    return None
