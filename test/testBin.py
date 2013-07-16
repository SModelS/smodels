#!/usr/bin/env python

"""
.. module:: testBin
   :synopsis: Runs all scripts in /bin. Asserts that no exceptions are thrown.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
"""
#import nose
import os


class Test(object):
    
    def __init__(self):
        self.directory = "../bin/"

    def testScripts(self):        
        for script in os.listdir(self.directory):
            if os.path.splitext(script)[1] != ".py":
                continue
            if script in ["run_all.py", "set_path.py", "RefXSecCalculator.py"]:
                continue
            yield self.runScript, script            
    
    def runScript(self, script):
        cmd = "python " + self.directory + script + " > /dev/null"
        ret = os.system(cmd)
        #self.assertEqual(ret, 0)
        assert ret == 0
  

if __name__ == "__main__":
    unittest.main()