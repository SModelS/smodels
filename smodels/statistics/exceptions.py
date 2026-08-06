"""
.. module:: statsExceptions
   :synopsis: Contains exceptions for SModelS's statistics package.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


class SModelSStatisticsError(Exception):
    """
    Class to define SModelS specific errors
    """

    def __init__(self, value: object = None):
        self.value = value
        Exception.__init__(self, value)

    def __str__(self) -> str:
        return repr(self.value)
