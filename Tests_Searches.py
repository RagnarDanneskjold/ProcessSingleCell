import unittest
import Searches

class testRangeBsearch(unittest.TestCase):

    search_data = [
            {'start':10,'end':50},    #0
            {'start':20,'end':80},    #1
            {'start':20,'end':100},   #2
            {'start':150, 'end':170}, #3
            {'start':170, 'end':200}, #4
            {'start':210,'end':220},  #5
            {'start':250, 'end':400}  #6
    ]

    case1 = [2,8]
    case2 = [5,15]
    case3 = [12,15]
    case4 = [15,25]
    case5 = [25,35]
    case6 = [80,90]
    case7 = [120,130]
    case8 = [100,150]
    case9 = [100,155]
    case10 = [160,165]
    case11 = [190,215]
    case12 = [300,400]
    case13 = [500,600]

    def testPass(self):
        self.assertTrue(True)

    def testFail(self):
        self.assertTrue(False)

    def testUpper(self):
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case1[0],self.case1[1],self.search_data),0)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case2[0],self.case2[1],self.search_data),1)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case3[0],self.case3[1],self.search_data),1)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case4[0],self.case4[1],self.search_data),3)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case5[0],self.case5[1],self.search_data),3)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case6[0],self.case6[1],self.search_data),3)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case7[0],self.case7[1],self.search_data),3)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case8[0],self.case8[1],self.search_data),3)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case9[0],self.case9[1],self.search_data),4)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case10[0],self.case10[1],self.search_data),4)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case11[0],self.case11[1],self.search_data),6)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case12[0],self.case12[1],self.search_data),7)
        self.assertEqual(Searches.rangeBsearchUpper(
            self.case13[0],self.case13[1],self.search_data),7)

    def testLower(self):
        self.assertEqual(Searches.rangeBsearchLower(
            self.case1[0],self.case1[1],self.search_data),0)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case2[0],self.case2[1],self.search_data),0)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case3[0],self.case3[1],self.search_data),0)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case4[0],self.case4[1],self.search_data),0)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case5[0],self.case5[1],self.search_data),0)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case6[0],self.case6[1],self.search_data),2)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case7[0],self.case7[1],self.search_data),3)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case8[0],self.case8[1],self.search_data),3)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case9[0],self.case9[1],self.search_data),3)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case10[0],self.case10[1],self.search_data),3)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case11[0],self.case11[1],self.search_data),4)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case12[0],self.case12[1],self.search_data),6)
        self.assertEqual(Searches.rangeBsearchLower(
            self.case13[0],self.case13[1],self.search_data),7)

    def testRangeAll(self):
        empty_data = []

        self.assertEqual(Searches.rangeBsearch(
            self.case1[0],self.case1[1],empty_data),[])

        self.assertEqual(Searches.rangeBsearch(
            self.case1[0],self.case1[1],self.search_data),[])
        self.assertEqual(Searches.rangeBsearch(
            self.case2[0],self.case2[1],self.search_data),[0,1])
        self.assertEqual(Searches.rangeBsearch(
            self.case3[0],self.case3[1],self.search_data),[0,1])
        self.assertEqual(Searches.rangeBsearch(
            self.case4[0],self.case4[1],self.search_data),[0,3])
        self.assertEqual(Searches.rangeBsearch(
            self.case5[0],self.case5[1],self.search_data),[0,3])
        self.assertEqual(Searches.rangeBsearch(
            self.case6[0],self.case6[1],self.search_data),[2,3])
        self.assertEqual(Searches.rangeBsearch(
            self.case7[0],self.case7[1],self.search_data),[])
        self.assertEqual(Searches.rangeBsearch(
            self.case8[0],self.case8[1],self.search_data),[])
        self.assertEqual(Searches.rangeBsearch(
            self.case9[0],self.case9[1],self.search_data),[3,4])
        self.assertEqual(Searches.rangeBsearch(
            self.case10[0],self.case10[1],self.search_data),[3,4])
        self.assertEqual(Searches.rangeBsearch(
            self.case11[0],self.case11[1],self.search_data),[4,6])
        self.assertEqual(Searches.rangeBsearch(
            self.case12[0],self.case12[1],self.search_data),[6,7])
        self.assertEqual(Searches.rangeBsearch(
            self.case13[0],self.case13[1],self.search_data),[])
        

unittest.main()
