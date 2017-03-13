import unittest
import Overlaps

class testOverlaps(unittest.TestCase):

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

    test10 = {
        'start':1,
        'end':11
        }

    test11 = {
        'start':19,
        'end':25
        }

    def testPass(self):
        self.assertTrue(True)

    def testFail(self):
        self.assertTrue(False)

    def testBooleanOverlap(self):
        self.assertTrue(Overlaps.overlap(self.test1, self.test2))
        self.assertTrue(Overlaps.overlap(self.test1, self.test3))
        self.assertTrue(Overlaps.overlap(self.test1, self.test4))
        self.assertTrue(Overlaps.overlap(self.test1, self.test5))
        self.assertFalse(Overlaps.overlap(self.test1, self.test6))
        self.assertFalse(Overlaps.overlap(self.test1, self.test7))
        self.assertFalse(Overlaps.overlap(self.test1, self.test8))
        self.assertFalse(Overlaps.overlap(self.test1, self.test9))
        self.assertTrue(Overlaps.overlap(self.test1, self.test10))
        self.assertTrue(Overlaps.overlap(self.test1, self.test11))

    def testOverlapBases(self):
        self.assertEqual(Overlaps.overlapBases(self.test1, self.test2),5)
        self.assertEqual(Overlaps.overlapBases(self.test1, self.test3),5)
        self.assertEqual(Overlaps.overlapBases(self.test1, self.test4),10)
        self.assertEqual(Overlaps.overlapBases(self.test1, self.test5),5)
        self.assertEqual(Overlaps.overlapBases(self.test1, self.test6),0)
        self.assertEqual(Overlaps.overlapBases(self.test1, self.test7),0)
        self.assertEqual(Overlaps.overlapBases(self.test1, self.test8),0)
        self.assertEqual(Overlaps.overlapBases(self.test1, self.test9),0)
        self.assertEqual(Overlaps.overlapBases(self.test1, self.test10),1)
        self.assertEqual(Overlaps.overlapBases(self.test1, self.test11),1)


unittest.main()
