"""
.. module:: inclusiveObjects
   :synopsis: Module holding multipurpose inclusive objects
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

import unum


class InclusiveValue(int):
    """
    An inclusive number class. It will return True when compared to any other integer, float or Unum object.
    """
    
    def __init__(self):
        int.__init__(self)
        
    def __str__(self):
        return '*'   
    
    def __repr__(self):
        return self.__str__()

    def __cmp__(self,other):
        if isinstance(other,(int,float,unum.Unum)):
            return 0
        else:
            return -1

    def __eq__(self,other):
        return self.__cmp__(other) == 0 
    
    def __ne__(self,other):
        return self.__cmp__(other) != 0
             
    
class InclusiveList(list):    
    """
    An inclusive list class. It will return True when compared to any other list object.
    """
    
    def __init__(self):
        list.__init__(self)
        
    def __str__(self):
        return '[*]'    

    def __repr__(self):
        return self.__str__()

    def __cmp__(self,other):
        if isinstance(other,list):
            return 0
        else:
            return -1

    def __eq__(self,other):
        return self.__cmp__(other) == 0  
    
    def __ne__(self,other):
        return self.__cmp__(other) != 0
