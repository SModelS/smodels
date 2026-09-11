"""
.. module:: types
   :synopsis: Common type aliases used throughout the smodels package.

.. moduleauthor:: SModelS developers
"""

from typing import TypeAlias
from os import PathLike

# very simple alias for anything that looks like a path
PathType: TypeAlias = str | PathLike[str]
