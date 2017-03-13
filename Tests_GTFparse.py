import unittest
import GTFparse
from copy import deepcopy

class testParseValidLine(unittest.TestCase):
    def isUnempty(self, data_test):
        if not data_test:
            return False

        return True

    def keysExist(self, data_test):
        if (
                'chrom' not in data_test or
                'feature' not in data_test or 
                'start' not in data_test or 
                'end' not in data_test or 
                'strand' not in data_test or
                'gene_id' not in data_test
                ):
            return False

        return True


    def compareData(self, data_true, data_test):
        if (    
                data_true['chrom'] == data_test['chrom'] and
                data_true['feature'] == data_test['feature'] and
                data_true['start'] == data_test['start'] and
                data_true['end'] == data_test['end'] and
                data_true['strand'] == data_test['strand'] and
                data_true['gene_id'] == data_test['gene_id']
                ):
            return True

        return False
  
    def testPass(self):
        self.assertTrue(True)

    def testFail(self):
        self.assertTrue(False)

    def testGeneString(self):
        gene_string = "chr1\tENSEMBL\tgene\t3102016\t3102125\t.\t+\t.\tgene_id \"ENSMUSG00000064842.1\"; gene_type \"snRNA\"; gene_status \"KNOWN\"; gene_name \"Gm26206\"; level 3;"

        true_data = {
            'chrom':'chr1',
            'feature':'gene',
            'start':3102016,
            'end':3102125,
            'strand':'+',
            'gene_id':'ENSMUSG00000064842.1'
            }

        test_data = GTFparse.parseValidLine(gene_string)
        self.assertTrue(self.isUnempty(test_data))
        self.assertTrue(self.keysExist(test_data))
        self.assertTrue(self.compareData(true_data, test_data))

    def testTranscriptString(self):
        transcript_string = "chr1\tENSEMBL\ttranscript\t3102016\t3102125\t.\t+\t.\tgene_id \"ENSMUSG00000064842.1\"; transcript_id \"ENSMUST00000082908.1\"; gene_type \"snRNA\"; gene_status \"KNOWN\"; gene_name \"Gm26206\"; transcript_type \"snRNA\"; transcript_status \"KNOWN\"; transcript_name \"Gm26206-201\"; level 3; transcript_support_level \"NA\"; tag \"basic\";"

        true_data = {
            'chrom':'chr1',
            'feature':'transcript',
            'start':3102016,
            'end':3102125,
            'strand':'+',
            'gene_id':'ENSMUSG00000064842.1'
            }

        test_data = GTFparse.parseValidLine(transcript_string)
        self.assertFalse(self.isUnempty(test_data))

class testSortGTFStructure(unittest.TestCase):
    unsort_data = {
        'chr1': {
            '+': [
                {
                    'gene_id':'GENE101',
                    'start':30,
                    'end':50,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE001',
                    'start':25,
                    'end':60,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE201',
                    'start':40,
                    'end':70,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE701',
                    'start':10,
                    'end':20,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE301',
                    'start':750,
                    'end':80,
                    'reads':list()
                    }
                ],
            '-': [
                {
                    'gene_id':'GENE101',
                    'start':30,
                    'end':50,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE001',
                    'start':25,
                    'end':60,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE201',
                    'start':40,
                    'end':70,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE701',
                    'start':10,
                    'end':20,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE301',
                    'start':750,
                    'end':80,
                    'reads':list()
                    }
                ]
            },
        'chr5': {
            '+': [
                {
                    'gene_id':'GENE101',
                    'start':30,
                    'end':50,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE001',
                    'start':25,
                    'end':60,
                    'reads':list()
                    }
                ],
            '-': [
                {
                    'gene_id':'GENE101',
                    'start':30,
                    'end':50,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE001',
                    'start':25,
                    'end':60,
                    'reads':list()
                    }
                ]
            },
        'chr3': {
            '+': [
                {
                    'gene_id':'GENE101',
                    'start':30,
                    'end':50,
                    'reads':list()
                    }
                ],
            '-': [
                ]
            }

        }

    sort_data = {
        'chr1': {
            '+': [
                {
                    'gene_id':'GENE701',
                    'start':10,
                    'end':20,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE001',
                    'start':25,
                    'end':60,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE101',
                    'start':30,
                    'end':50,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE201',
                    'start':40,
                    'end':70,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE301',
                    'start':75,
                    'end':80,
                    'reads':list()
                    }
                ],
            '-': [
                {
                    'gene_id':'GENE701',
                    'start':10,
                    'end':20,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE001',
                    'start':25,
                    'end':60,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE101',
                    'start':30,
                    'end':50,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE201',
                    'start':40,
                    'end':70,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE301',
                    'start':75,
                    'end':80,
                    'reads':list()
                    }
                ]
            },
        'chr3': {
            '+': [
                {
                    'gene_id':'GENE101',
                    'start':30,
                    'end':50,
                    'reads':list()
                    }
                ],
            '-': [
                ]
            },
        'chr5': {
            '+': [
                {
                    'gene_id':'GENE001',
                    'start':25,
                    'end':60,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE101',
                    'start':30,
                    'end':50,
                    'reads':list()
                    }
                ],
            '-': [
                {
                    'gene_id':'GENE001',
                    'start':25,
                    'end':60,
                    'reads':list()
                    },
                {
                    'gene_id':'GENE101',
                    'start':30,
                    'end':50,
                    'reads':list()
                    }
                ]
            }

        }

    def compare_sorts(self, data_true, data_test):
        for chrom in data_test:
            for strand in data_test[chrom]:
                for i in range(len(data_test[chrom][strand])):
                    if data_test[chrom][strand][i]['gene_id'] != data_true[chrom][strand][i]['gene_id']:
                        return False

        return True

    def testPass(self):
        self.assertTrue(True)

    def testFail(self):
        self.assertTrue(False)


    def testBase(self):
        self.assertFalse(self.compare_sorts(self.unsort_data, self.sort_data))
        self.assertTrue(self.compare_sorts(self.sort_data, self.sort_data))

        testset = deepcopy(self.sort_data)
        self.assertTrue(self.compare_sorts(testset, self.sort_data))

    def testUnsort(self):
        testset = deepcopy(self.unsort_data)
        GTFparse.sortGTFStructure(testset)
        self.assertTrue(self.compare_sorts(testset, self.sort_data))

    def testSort(self):
        testset = deepcopy(self.sort_data)
        GTFparse.sortGTFStructure(testset)
        self.assertTrue(self.compare_sorts(testset, self.sort_data))

class testParseGTFFile(unittest.TestCase):
    filename = "TEST.gtf"

    data_true = {
        'chr4': {
            '+': [
                {
                    'gene_id':'NAME1',
                    'start':10,
                    'end':20,
                    'reads':list()
                    }

                ],
            '-': [
                {
                    'gene_id':'NAME4',
                    'start':1,
                    'end':6,
                    'reads':list()
                    },
                {
                    'gene_id':'NAME3',
                    'start':5,
                    'end':15,
                    'reads':list()
                    }
                ]
            },
        'chr2': {
            '+':[
                ],
            '-':[
                {
                    'gene_id':'NAME2',
                    'start':35,
                    'end':70,
                    'reads':list()
                    }
                ]

            }

        }

    def make_testfile(self):
        fp = open(self.filename, "w")
        header = "##HEADERLINE1\n##HEADERLINE2\n"
        line1 = "chr4\tHAVANA\tgene\t10\t20\t.\t+\t.\tgene_id \"NAME1\";\n"
        line2 = "chr2\tHAVANA\tgene\t35\t70\t.\t-\t.\tgene_id \"NAME2\";\n"
        line3 = "chr2\tHAVANA\ttranscript\t35\t70\t.\t+\t.\tgene_id \"NAME2.1\"; transcript_id \"TRANSCRIPT2\";\n"
        line4 = "chr4\tENSEMBL\tgene\t5\t15\t.\t-\t.\tgene_id \"NAME3\";\n"
        line5 = "chr4\tENSEMBL\tgene\t1\t6\t.\t-\t.\tgene_id \"NAME4\";\n"
        
        fp.write(header + line1 + line2 + line3 + line4 + line5)
        fp.close()

    def compare_data(self,data_true, data_test):
        for chrom in data_true:
            if chrom not in data_test:
                return False
            for strand in data_true[chrom]:
                if strand not in data_test[chrom]:
                    return False
                if len(data_true[chrom][strand]) != len(data_test[chrom][strand]):
                    return False
                for i in range(len(data_true[chrom][strand])):
                    if data_true[chrom][strand][i]['gene_id'] != data_test[chrom][strand][i]['gene_id']:
                        return False
                    if data_true[chrom][strand][i]['start'] != data_test[chrom][strand][i]['start']:
                        return False
                    if data_true[chrom][strand][i]['end'] != data_test[chrom][strand][i]['end']:
                        return False
                    if data_true[chrom][strand][i]['reads'] != data_test[chrom][strand][i]['reads']:
                        return False

        return True

    def testPass(self):
        self.assertTrue(True)

    def testFail(self):
        self.assertTrue(False)

    def testFile(self):
        self.make_testfile()
        fp = open(self.filename, "r")
        data_test = GTFparse.parseGTFFile(fp)
        fp.close()

        self.assertTrue(self.compare_data(self.data_true, data_test))

unittest.main()



