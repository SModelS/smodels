"""
.. module:: experimentExceptions
   :synopsis: Contains exceptions for SModelS's experiment package.

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

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


class SModelSExperimentError(Exception):
    """
    Class to define SModelS specific errors
    """

    def __init__(self, value=None):
        self.value = value
        Exception.__init__(self, value)

    def __str__(self):
        return repr(self.value)
