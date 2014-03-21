#!/usr/bin/env python

"""
.. module:: asciiGraphs
   :synopsis: This unit contains a simple routine to draw ASCII-art Feynman
   graphs.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
import logging

logger = logging.getLogger(__name__)

def printParticle_(label):
    """
    Renames a few particles for the asciidraw routine. Do not call directly.
    
    """
    if label == "jet":
        label = "q"
    label = label + "     "
    return label[:2]

def drawBranch_(branch, upwards, labels, htmlFormat, border, l):
    """
    Draws a single branch. Should only be used via .asciidraw, not directly.
    
    """
    length = 0
    lines=["   ","----"] 
    labels="   "
    if border and upwards:
        lines=[" |    "," | ----"]
        labels=" |    "
    if border and not upwards:
        lines=[" |    "," | ----"]
        labels=" |    "

    for insertions in branch.particles:
        if len(insertions) == 0:
            lines[0] += " "
            lines[1] += "*"
            continue
        lines[1] += "*----"
        if len(insertions) == 1:
            labels += " " + printParticle_(insertions[0]) + "  "
            lines[0] += " |   "
        if len(insertions) == 2:
            labels += printParticle_(insertions[0]) + " " + \
                    printParticle_(insertions[1])
            if upwards:
                lines[0]+="\\ /  "
            else:
                lines[0]+="/ \\  "
        if len(insertions) > 2:
            logger.error("case for n-body decay, n>3 not yet. implemented. \
                          Please implement.")
            sys.exit(0)

    order = [0, 1]
    if not upwards: order = [1, 0]
    html = "<br>"
    lengthdiff = l - len(lines[0])/5
    if border:
        if l == 2:
            lines[0] += " "
            lines[1] += " "
            labels += " "
        labels += " " + " "*(5*lengthdiff) + " |"
        lines[0] += " "*(5*lengthdiff + 0) + "  |"
        lines[1] += " "*(5*lengthdiff + 0) + " |"
    if border and upwards:
        print " /" + "-"*(4*l + 4) + "\\"
    if htmlFormat:
        print html
    if upwards and labels:
        print labels
    if htmlFormat:
        print html
    for i in order:
        print lines[i]
    if htmlFormat:
        print html
    if not upwards and labels:
        print labels
    if htmlFormat:
        print html
    if border and not upwards:
        print " \\" + "-"*(4*l + 4) + "/"

def asciidraw(element, labels=True, html=False, border=False):
    """
    Draws a simple ASCII graph on the screen.
    
    """
    l = []
    for (ct, branch) in enumerate(element.branches):
        l.append(int( str(branch).count("[")))
    for (ct, branch) in enumerate(element.branches):
        drawBranch_(branch, upwards=(ct == 0), labels=labels, htmlFormat=html,
                    border=border, l=max(l))

if __name__ == "__main__":
    import set_path, argparse, types
    import SModelS
    from theory import LHEReader, lheDecomposer, crossSection

    argparser = argparse.ArgumentParser(description = 'simple tool that is \
            meant to draw lessagraphs, as an ascii plot') 
    argparser.add_argument('-T', nargs='?', help = 'Tx name, will look up lhe \
            file in ../regression/Tx_1.lhe. Will be overriden by the "--lhe" \
            argument', type=types.StringType, default='T1')
    argparser.add_argument('-l', '--lhe', nargs='?', help = 'lhe file name, \
            supplied directly. Takes precedence over "-T" argument.', 
            type=types.StringType, default='')
    argparser.add_argument('-b', '--border', help='draw a border around the \
            graph', action='store_true')
    args = argparser.parse_args()

    filename = "%s/lhe/%s_1.lhe" % (SModelS.installDirectory(), args.T)
    if args.lhe != "":
        filename = args.lhe

    reader = LHEReader.LHEReader(filename)
    Event = reader.next()
    element = lheDecomposer.elementFromEvent(Event,
                                             crossSection.XSectionList())
    asciidraw(element, border=args.border)                    
