import unittest
import Overlaps

class testOverlap(unittest.TestCase):

    test1 = {
        'start':10,
        'end':20
        }

    test2 = {
        'start':5,
        'end':15
        }

    test3 = {
        'start':15,
        'end':25
        }

    test4 = {
        'start':5,
        'end':25
        }

    test5 = {
        'start':12,
        'end':17
        }

    test6 = {
        'start':1,
        'end':5
        }

    test7 = {
        'start':25,
        'end':30
        }

    test8 = {
        'start':1,
        'end':10
        }

    test9 = {
        'start':20,
        'end':25
        }

    def testPass(self):
        self.assertTrue(True)

    def testFail(self):
        self.assertTrue(False)

    def testAllCombinations(self):
        assertTrue(Overlap.overlapBases()

unittest.main()
