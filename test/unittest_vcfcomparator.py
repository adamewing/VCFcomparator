#!/usr/bin/env python

import unittest
import vcfcomparator as vc

class testVCFcomparator(unittest.TestCase):
    def setUp(self):
        self.vcf_list = ['test_data/test_A.vcf.gz','test_data/test_B.vcf.gz']
        self.vcf_handles = vc.openVCFs(self.vcf_list)
        self.comparison = vc.compareVCFs(self.vcf_handles[0], self.vcf_handles[1])
        self.matchedSNV = self.comparison.vartype['SNV'][0]
        self.matchedSV  = self.comparison.vartype['SV'][0]
        self.overlapSV = self.comparison.vartype['SV'][1]

    def testComparisonMatchedSNV(self):
        self.assertEqual(self.comparison.matched('SNV'), 2)

    def testComparisonMatchedSV(self):
        self.assertEqual(self.comparison.matched('SV'),2)

    def testComparisonAltMatchedSV(self):
        self.assertEqual(self.comparison.altmatched('SV'),1)

    def testMatchedSNVMatched(self):
        self.assertTrue(self.matchedSNV.matched())

    def testMatchedSNVScore(self):
        self.assertEqual(self.matchedSNV.score(), 1.0)

    def testMatchedSVMatched(self):
        self.assertTrue(self.matchedSV.matched())

    def testMatchedSVScore(self):
        self.assertEqual(self.matchedSV.score(), 1.0)

    def testOverlapSVScore(self):
        self.assertLess(self.overlapSV.score(),1.0)
        self.assertGreater(self.overlapSV.score(),0.0)

if __name__ == '__main__':
    unittest.main()
