#!/usr/bin/env python

import unittest
import vcfcomparator as vc

class testVCFcomparator(unittest.TestCase):
    def setUp(self):
        self.vcf_list = ['test_data/test_A.vcf.gz','test_data/test_B.vcf.gz']
        self.vcf_handles = vc.openVCFs(self.vcf_list)
        self.comparison = vc.compareVCFs(self.vcf_handles[0], self.vcf_handles[1])

        self.matchedSNV = self.comparison.vartype['SNV'][0]
        self.matchedPFSNV = self.comparison.vartype['SNV'][1] # disagreement on filter and somatic
        self.unmatchedSNV = self.comparison.vartype['SNV'][2]

#        self.matchedINDEL = self.comparison.vartype['INDEL'][0]
#        self.unmatchedINDEL = self.comparison.vartype['INDEL'][1]

        self.matchedSV  = self.comparison.vartype['SV'][0]
        self.overlapSV = self.comparison.vartype['SV'][1]

#        self.overlapCNV = self.comparison.vartype['CNV'][0]
#        self.unmatchedCNV = self.comparison.vartype['CNV'][1]

    ## SNV tests ##

    def testComparisonMatchedSNV(self):
        self.assertEqual(self.comparison.matched('SNV'), 3)

    def testMatchedSNVMatched(self):
        self.assertTrue(self.matchedSNV.matched())

    def testMatchedSNVUnmatched(self):
        self.assertFalse(self.unmatchedSNV.matched())

    def testMatchedSNVScore(self):
        self.assertEqual(self.matchedSNV.score(), 1.0)

    def testMatchedPPSNV(self): # both SNVs pass filter
        self.assertTrue(self.matchedSNV.both_pass())

    def testMatchedPFSNV(self): # one pass SNV filter
        self.assertTrue(self.matchedPFSNV.has_pass())
        self.assertFalse(self.matchedPFSNV.both_pass())

    def testCountPassAgree(self):
        p = self.comparison.count_agree_pass('SNV')
        self.assertGreater(p,0.0)

    def testCountFailAgree(self):
        p = self.comparison.count_agree_fail('SNV')
        self.assertEqual(p,0.0)

    def testCountPassDisagree(self):
        p = self.comparison.count_disagree_pass('SNV')
        self.assertGreater(p,0.0)

    def testCountSomaticAgree(self):
        p = self.comparison.count_agree_somatic('SNV')
        self.assertGreater(p,0.0)

    def testCountSomaticDisagree(self):
        p = self.comparison.count_disagree_somatic('SNV')
        self.assertGreater(p,0.0)

    def testSNVSummarySumScores(self):
        s = self.comparison.sum_scores('SNV')
        self.assertGreater(s,0.0)

    ## INDEL tests ##

    ## SV tests ##

    def testComparisonMatchedSV(self):
        self.assertEqual(self.comparison.matched('SV'),2)

    def testComparisonAltMatchedSV(self):
        self.assertEqual(self.comparison.altmatched('SV'),1)

    def testMatchedSVMatched(self):
        self.assertTrue(self.matchedSV.matched())

    def testMatchedSVScore(self):
        self.assertEqual(self.matchedSV.score(), 1.0)
        self.assertEqual(self.matchedSV.score(density=True), 1.0)

    def testOverlapSVScore(self):
        self.assertLess(self.overlapSV.score(),1.0)
        self.assertLess(self.overlapSV.score(density=True),1.0)
        self.assertGreater(self.overlapSV.score(),0.0)
        self.assertGreater(self.overlapSV.score(density=True),0.0)

        # turning density on should decrease the score for a partial overlap
        self.assertLess(self.overlapSV.score(density=True),self.overlapSV.score(density=False))

    ## CNV tests ##

    def testSVSummarySumScores(self):
        s_w = self.comparison.sum_scores('SV', weight=True)
        s_uw = self.comparison.sum_scores('SV', weight=False)
        self.assertGreater(s_w,0.0)
        self.assertGreater(s_uw,0.0)
        self.assertFalse(s_uw == s_w)

if __name__ == '__main__':
    unittest.main()
