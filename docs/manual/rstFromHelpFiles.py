#!/usr/bin/python

from __future__ import print_function

import commands

def isEmpty ( string ):
    return string.strip() == ""

def run ( tool, rstfile ):
    with open("source/%s.rst" % rstfile, "w" ) as f:
        c = commands.getoutput ( "../../smodelsTools.py %s -h" % tool )
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
            pnew = new.replace ( "\n", "" )
            if False:
                print ( "isIn='%s' lastIn='%s' now write: >>%s<< hasNL=%d" % \
                           ( isIn, lastIn, pnew[:30], hasNewLine ) )
            f.write ( new )
            lastIn = isIn
            if line[:3] == "  -":
                isIn = "argument"


run ( "xseccomputer", "XSecComputer" )
