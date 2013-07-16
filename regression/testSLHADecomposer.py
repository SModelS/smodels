#!/usr/bin/env python

"""
.. module:: testSLHADecomposer
   :synopsis: Compares output of SLHA decomposition with given value.
    
.. moduleauthor:: Wolfgang magerl <wolfgang.magerl@gmail.com>
    
"""
import unittest
import setPath
import pickle, math
import Theory.SLHADecomposer as DEC
import Tools.PhysicsUnits as UNIT

class Test(unittest.TestCase):


    def testSLHADecomposition(self):
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
            if g.isEqual( gtopo[i] ):                     #compares GTop
                print 'GTop',g, '%sOK!%s' %( green, reset )
            else:
                print 'GTop',g,'%sfailed!%s' %( red, reset ) 
            j = 0
            for el in g.ElList:
                if el.isEqual( gtopo[i].ElList[j] ):        #compares EElement
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
                    if b.isEqual( gtopo[i].ElList[j].B[k] ):         #compares BElement
                        print 'BElement',b, '%sOK!%s' %( green, reset )
                    else:
                        print 'BElement',b, '%sfailed!%s' %( red, reset ) 
                    k +=1
                j += 1
            i += 1
            
        with open('4.log', 'r') as logFile:
            data = logFile.read().replace('\n', '')
        
        self.assertEqual(resultsMerge, data)


if __name__ == "__main__":
    unittest.main()