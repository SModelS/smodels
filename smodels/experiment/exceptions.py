"""
.. module:: experiment.exceptions
   :synopsis: Contains exceptions for SModelS's experiment package.

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

class DatabaseNotFoundException(Exception):
    """
    This exception is used when the database cannot be found.
    
    """
    def __init__(self, value):
        self.value = value
        Exception.__init__(self, value)
        

    def __str__(self):
        return repr(self.value)
