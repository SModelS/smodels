#!/usr/bin/env python3

"""
.. module:: asciiGraph
   :synopsis: Contains a simple routine to draw ASCII-art Feynman-like graphs.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
import sys
import argparse
from smodels.theory import decomposer
from smodels.tools.smodelsLogging import logger
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model



def _printParticle(label):
    """
    Rename particles for the asciidraw routine.
    
    """
    if label == "jet":
        label = "q"
    label = label + "     "
    return label[:2]


def _drawBranch(branch, upwards, htmlFormat, border, l):
    """
    Draw a single branch.
    
    """
    ret=""
    lines = ["   ", "----"]
    labels = "   "
    if border and upwards:
        lines = [" |    ", " | ----"]
        labels = " |    "
    if border and not upwards:
        lines = [" |    ", " | ----"]
        labels = " |    "

    for insertions in branch.evenParticles:
        if len(insertions) == 0:
            lines[0] += " "
            lines[1] += "*"
            continue
        lines[1] += "*----"
        if len(insertions) == 1:
            labels += " " + _printParticle(insertions[0].label) + "  "
            lines[0] += " |   "
        if len(insertions) == 2:
            labels += _printParticle(insertions[0].label) + " " + \
                    _printParticle(insertions[1].label)
            if upwards:
                lines[0] += "\\ /  "
            else:
                lines[0] += "/ \\  "
        if len(insertions) > 2:
            logger.error("n > 3 for n-body decay not yet implemented.")
            sys.exit(0)

    order = [0, 1]
    if not upwards:
        order = [1, 0]
    html = "<br>"
    lengthdiff = int ( l - len(lines[0]) / 5 )
    if border:
        if l == 2:
            lines[0] += " "
            lines[1] += " "
            labels += " "
        labels += " " + " "*(5 * lengthdiff) + " |"
        lines[0] += " "*(5 * lengthdiff + 0) + "  |"
        lines[1] += " "*(5 * lengthdiff + 0) + " |"
    if border and upwards:
        ret+=" /" + "-"*(4 * l + 4) + "\\\n"
    if htmlFormat:
        ret+=html+"\n"
    if upwards and labels:
        ret+=labels+"\n"
    if htmlFormat:
        ret+=html+"\n"
    for i in order:
        ret+=lines[i]+"\n"
    if htmlFormat:
        ret+=html+"\n"
    if not upwards and labels:
        ret+=labels+"\n"
    if htmlFormat:
        ret+=html+"\n"
    if border and not upwards:
        ret+=" \\" + "-"*(4 * l + 4) + "/\n"
    return ret


def asciidraw(element, labels=True, html=False, border=False):
    """
    Draw a simple ASCII graph on the screen.
    
    """
    ret=""
    l = []
    for (ct, branch) in enumerate(element.branches):
        l.append(int(str(branch).count("[")))
    for (ct, branch) in enumerate(element.branches):
        ret+=_drawBranch(branch, upwards=(ct == 0), htmlFormat=html,
                    border=border, l=max(l))
    return ret


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="simple tool that is "
                                        "meant to draw lessagraphs, as an "
                                        "ascii plot")
    argparser.add_argument('-l', '--lhe', help="LHE file name",
                           type=str, required=True )
    argparser.add_argument('-b', '--border', action='store_true',
                           help="draw a border around the graph")
    args = argparser.parse_args()

    filename = args.lhe

    model = Model(BSMparticles=BSMList, SMparticles=SMList)
    model.updateParticles(inputFile=filename)
    topList = decomposer.decompose(model)
    element = topList.getElements()[0]
    print(asciidraw(element, border=args.border))
