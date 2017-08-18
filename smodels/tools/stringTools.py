"""
.. module:: stringTools
   :synopsis: Holds all code snippets that meddle with strings

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from smodels.tools.smodelsLogging import logger

def cleanWalk ( topdir ):
    """ perform os.walk, but ignore all hidden files and directories """
    import os
    ret = []
    for root, d_, f_ in os.walk ( topdir ):
        isHidden=False
        tokens = root.split("/")
        for token in tokens:
            if len(token)>0 and token[0]==".":
                isHidden=True
                break
        if isHidden:
            continue
        dirs,files = [],[]
        for d in d_:
            if not d[0]==".":
                dirs.append ( d )
        for f in f_:
            if not f[0]==".":
                files.append ( f )
        ret.append ( [ root, dirs, files ] )
    return ret

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
