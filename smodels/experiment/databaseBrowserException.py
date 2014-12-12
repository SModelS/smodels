"""
.. module:: databaseException
   :synopsis: Exception for databaseBrowser.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""
class DatabaseNotFoundException(Exception):
    pass

class InvalidExperimentException(Exception):
    pass

class InvalidRunRestrictionException(Exception):
    pass

class InvalidRunException(Exception):
    pass

class InvalidExperimentIDException(Exception):
    pass

class InvalidTxNameException(Exception):
    pass

class ROOTNotFoundException(Exception):
    pass

class MassParametrizationException(Exception):
    pass

class InvalidInfotxtFileException(Exception):
    pass

class InvalidSmspytFileException(Exception):
    pass

class InvalidInfoFieldException(Exception):
    pass

class InvalidFieldValueException(Exception):
    pass