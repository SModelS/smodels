#!/usr/bin/env python

"""
.. module:: simpleExample
   :synopsis: Basic use case for the SModelS framework.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import elementtree.ElementTree as ET

smstoplist = [1,2,3,4,55,6]
toplist = ET.Element("toplist")
for itop,topi in enumerate(smstoplist):
  top = ET.SubElement(toplist,"top")
  top.attrib["id"] = "top"+str(itop)
  nels = ET.SubElement(top,"graph")
  nels.text = str(topi)
  
tree = ET.ElementTree(toplist)
tree.write("text.xml")
