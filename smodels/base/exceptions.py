"""
.. module:: baseExceptions
   :synopsis: Contains exceptions for SModelS's base packages.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


class SModelSBaseError(Exception):
    """
    Class to define SModelS specific errors
    """

    def __init__(self, value=None):
        self.value = value
        Exception.__init__(self, value)

    def __str__(self):
        return repr(self.value)
