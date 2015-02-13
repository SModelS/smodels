"""
.. module:: stringTools
   :synopsis: Holds all code snippets that meddle with strings

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def concatenateLines ( oldcontent ):
    """ of all lines in the list "oldcontent", concatenate the ones
        that end with \  or , """
    content=[] ## concatenate lines that end with "," or "\"
    tmp=""
    for line in oldcontent:
        tmp+=line.strip()
        if tmp != "" and tmp[-1] not in [ ",", '\\' ]:
            content.append ( tmp )
            tmp=""
        if tmp != "" and tmp[-1] == '\\':
            tmp=tmp[:-1] # remove trailing \ (but keep trailing ,)
    return content
