"""
.. module:: decompositionExceptions
   :synopsis: Contains exceptions for SModelS's theory package.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from typing import Optional


class SModelSDecompositionError(Exception):
    """
    Class to define SModelS specific errors
    """

    def __init__(self, value: Optional[str] = None):
        self.value: Optional[str] = value
        Exception.__init__(self, value)

    def __str__(self) -> str:
        return repr(self.value)
