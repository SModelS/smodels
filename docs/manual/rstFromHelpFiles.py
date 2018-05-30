#!/usr/bin/env python3

from __future__ import print_function

import subprocess

def isEmpty ( string ):
    return string.strip() == ""

def write ( c, f ):
    lines = c.split("\n")
    isIn = ""
    lastIn = isIn
    for line in lines:
        new=line.replace ( "optional arguments", "*arguments*" )
        new=new.replace ( "usage: ", "   " )
        if line[:3] == "  -":
            isIn = "newargument"
            while "[" in new:
                p1 = new.find ( "[" )
                p2 = new.find ( "]" )
                new = new[:p1]+new[p2+1:]
            while new[-1] == " ":
                new=new[:-1]
            while " ," in new:
                new = new.replace ( " ,", "," )
        if isIn == "afterusage" and not isEmpty(line):
            isIn = "afterdescription"
            continue
        if isIn == "usage" and isEmpty ( line ):
            isIn = "afterusage"
        if line[:5] == "usage":
            isIn = "usage" 
        if not isIn in [ "usage", "argument" ]:
            new += "\n"
        if isIn == "newargument" and lastIn == "argument":
            new = "\n" + new
        hasNewLine = new[-1:] == "\n"
        if False:
            pnew = new.replace ( "\n", "" )
            print ( "isIn='%s' lastIn='%s' now write: >>%s ...<< hasNL=%d" % \
                       ( isIn, lastIn, pnew[:30], hasNewLine ) )
        f.write ( new )
        lastIn = isIn
        if line[:3] == "  -":
            isIn = "argument"

def run ( tool, rstfile ):
    with open("source/%s.rst" % rstfile, "w" ) as f:
        c = subprocess.getoutput ( "../../smodelsTools.py %s -h" % tool )
        write ( c, f )

def runSModelS ():
    with open("source/RunSModelS.rst", "w" ) as f:
        c = subprocess.getoutput ( "../../runSModelS.py -h" )
        write ( c, f )


run ( "xseccomputer", "XSecComputer" )
run ( "lhechecker", "LheChecker" )
run ( "slhachecker", "SlhaChecker" )
run ( "database-browser", "DatabaseBrowser" )
run ( "toolbox", "ToolBox" )
run ( "fixpermissions", "FixPermissions" )
runSModelS()
