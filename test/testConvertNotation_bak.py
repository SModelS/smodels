#!/usr/bin/env python3

"""
.. module:: convertNotationTest

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
from testbook import testbook
import unittest

notebook = '../notebooks-Examples/convertNotation.ipynb'

class convertNotationTest(unittest.TestCase):

    def testConvertNotation(self):

        importCells = [0]
        testCells = [[2],[3],[4],[4,5,6],[7]]
        expectedOutputs = ['(PV > anyOdd(1),anyOdd(2)), (anyOdd(1) > anyOdd(3),e-,nu), (anyOdd(3) > MET,jet,jet), (anyOdd(2) > MET,L,nu)',
        '(PV > squark(1),squark(2)), (squark(1) > gluino(3),e-,nu), (gluino(3) > MET,jet,jet), (squark(2) > HSCP,L,nu)',
        '(PV > gluino(1),squark(2)), (gluino(1) > MET,jet,jet), (squark(2) > HSCP,L,nu)',
        '[(PV, gluino), (PV, squark), (gluino, MET), (gluino, jet), (gluino, jet), (squark, HSCP), (squark, L), (squark, nu)]',
        '[(PV, gluino), (PV, squark), (gluino, X), (gluino, X), (gluino, jet), (squark, HSCP), (squark, L), (squark, nu), (X, nu), (X, nu), (X, MET), (X, e-)]']

        for ic,cell in enumerate(testCells):
            with testbook(notebook, execute=importCells + cell) as tb:
                output = tb.ref("output")
                # print(output)
                self.assertTrue(output == expectedOutputs[ic])




if __name__ == "__main__":
    unittest.main()
