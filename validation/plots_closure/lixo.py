#!/usr/bin/env python

import sys


f1 = open('T2-ATLAS.dat','r')
f2 = open('bak/T2_.dat','r')
data1 = f1.read()
data2 = f2.read()
pts1 = data1[:data1.find('#END')-1].split('\n')
pts2 = data2[:data2.find('#END')-1].split('\n')


for ipt,pt1 in enumerate(pts1):
  x,y,res,lim,cond,tot = pt1.split()
  pt2 = pts2[ipt]
  x2,y2,res2,lim2 = pt2.split()
  if x != x2: print x,y,'x differ',x,x2
  if y != y2: print x,y,'y differ',y,y2
  if lim != lim2: print x,y,'lim differ',lim,lim2
  if res != res2: print x,y,'res differ',res,res2
  
sys.exit()  




