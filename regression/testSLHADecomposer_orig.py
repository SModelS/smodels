#!/usr/bin/python

"""regression test for SLHADecomposer """

import pickle, math
import set_path
import Theory.SLHADecomposer as DEC
import Tools.PhysicsUnits as UNIT

i = 0
j = 0
k = 0
red = '\033[0;31m'
green = '\033[0;32m'
reset = '\033[0;0m'

dic = pickle.load(file("dictT2.pkl"))
sigmacut = UNIT.addunit(0.1,'fb')   
slhafile = "../slha/T2_600_375_xsec.slha"

cs = dic.crossSections()
gtopo = pickle.load(file("TopologyList.pkl"))
gtopo_new = DEC.decompose( slhafile ) 

for g in gtopo_new:
#  print 'gtopology:',g, g.describe()
  if g.isEqual( gtopo[i] ):   #compares GTop
    print 'GTop',g, '%sOK!%s' %( green, reset )
  else:
    print 'GTop',g,'%sfailed!%s' %( red, reset ) 
  j = 0
  for el in g.ElList:
#    print 'element:', el, el.describe()
#    print 'cross section:', el.weight
    if el.isEqual( gtopo[i].ElList[j] ):    #compares EElement
      print 'EElement',el, '%sOK!%s' %( green, reset )
    else:
      print 'EElement',el, '%sfailed!%s' %( red, reset )
    for key in el.weight:
      if math.fabs(UNIT.rmvunit(el.weight[key], 'fb')-UNIT.rmvunit(gtopo[i].ElList[j].weight[key],'fb'))> 0.00001:
        print 'different weight for EElement', el, '.%sfailed!%s' %( red, reset )
        print 'weight_orig:', gtopo[i].ElList[j].weight
        print 'weight_new:', el.weight
      else:
        print 'weight for EElement',el,'%sOK!%s' %( green, reset ) 
    k = 0
    for b in el.B:
#      print 'branch:', b, b.describe()
#      print 'branch_old:', gtopo[i].ElList[j].B[k], gtopo[i].ElList[j].B[k].describe()
#      print 'particles branch:', b.particles
#      print 'particles branch_old:', gtopo[i].ElList[j].B[k].particles
      if b.isEqual( gtopo[i].ElList[j].B[k] ):     #compares BElement
        print 'BElement',b, '%sOK!%s' %( green, reset )
      else:
        print 'BElement',b, '%sfailed!%s' %( red, reset ) 
      k +=1
    j += 1
  i += 1

