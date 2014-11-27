#!/usr/bin/env python

"""
.. module:: runCompleteTestSuite
   :synopsis: Runs all test suites.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com> 
    
"""
import sys
sys.path.append("../")

import unittest

def main():
    unittest.TextTestRunner().run( unittest.TestLoader().discover("./") )

if __name__ == "__main__":
    main()
