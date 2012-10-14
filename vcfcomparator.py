#!/usr/bin/env python

'''
VCF comparator: compares mutation calls in VCF formatted files
(c) 2012 Adam Ewing (ewingad at soe.ucsc.edu)
see LICENSE.txt for licensing inforomation
'''

import vcf
import argparse
import gzip
import sys
from re import search
from os.path import exists
from scipy.stats import norm
from numpy import interp

## classes ##

class Comparison:
    ''' stores the result of a one-way comparison vcfA --> vcfB
        imprements functions to report things about the comparison '''
    def __init__(self):
        self.vartype = {}
        self.vartype['SNV']   = []
        self.vartype['INDEL'] = []
        self.vartype['CNV']   = []
        self.vartype['SV']    = []

    def summary(self):
        for type in self.vartype.keys():
            print type, self.matched(type), "of", len(self.vartype[type]), "alt", self.altmatched(type)

    def matched(self, type):
        ''' count number of matched variants of <type> '''
        return reduce(lambda x, y: x+y.matched(), self.vartype[type], 0)

    def altmatched(self,type):
        ''' count number of alternate matches '''
        return reduce(lambda x, y: x+len(y.altmatch), self.vartype[type], 0)

    def count_agree_somatic(self,type):
        ''' count number of matches that agree on somatic status given variant type '''
        return reduce(lambda x, y: x+y.both_somatic(), self.vartype[type], 0)

    def count_disagree_somatic(self,type):
        ''' count number of matches that disagree on somatic status given variant type '''
        return reduce(lambda x, y: x+(y.matched() and y.has_somatic() and not y.both_somatic()), self.vartype[type], 0)

    def count_agree_pass(self,type):
        ''' count number of matches that agree on passing filters given variant type'''
        return reduce(lambda x, y: x+y.both_pass(), self.vartype[type], 0)

    def count_agree_fail(self,type):
        ''' count number of matches that agree on not passing filters given variant type '''
        return reduce(lambda x, y: x+(y.has_pass()==False), self.vartype[type], 0)

    def count_disagree_pass(self,type):
       ''' count number of matches where one call passed and the other did not given variant type '''
       return reduce(lambda x, y: x+(y.matched() and y.has_pass() and not y.both_pass()), self.vartype[type], 0)

    def sum_scores(self,type):
        ''' return sum of all scores for variant type '''
        return reduce(lambda x, y: x+y.score(), self.vartype[type], 0.0)

class Variant:
    ''' base class for variant types 
        vcf_recA and vcf_recB are vcf._Record objects or None '''
    def __init__(self, vcf_recA, vcf_recB, somatic=False):
        self.recA = vcf_recA
        self.recB = vcf_recB

        # if there is more than one match, the rest are stored here
        self.altmatch = []

    def __str__(self):
        return str(self.recA) + "\t" + str(self.recB)

    def interval_score(self, density=False):
        ''' scoring function for intervals, based on amount of overlap '''
        # get pos,end +/- confidence intervals if present
        iv_a = self.get_conf_interval(self.recA)
        iv_b = self.get_conf_interval(self.recB)

        ol_coords = get_overlap_coords(iv_a, iv_b)
        ol_width = ol_coords[1] - ol_coords[0]
        assert ol_width > 0

        len_a = iv_a[1] - iv_a[0]
        len_b = iv_b[1] - iv_b[0]
        assert len_a > 0
        assert len_b > 0


        s = float(2*ol_width)/float(len_a+len_b)
        if density: # weight overlap based on cumulative density of normal distribution
            den_a = 1.0
            den_b = 1.0
            if ol_width < len_a:
                den_a = interval_density(iv_a, ol_coords)
            if ol_width < len_b:
                den_b = interval_density(iv_b, ol_coords)
            s = s * den_a * den_b

        return s

    def get_conf_interval(self,rec):
        ''' return confidence interval as (start-ci, end+ci)'''
        cipos = ciend = 0
        
        if 'CIPOS' in rec.INFO:
            cipos = abs(min(rec.INFO.get('CIPOS')))
        if 'CIEND' in rec.INFO:
            ciend = abs(max(rec.INFO.get('CIEND')))

        end = rec.INFO.get('END')   
        if not end:
            end = rec.POS

        return rec.POS-cipos, end+ciend

    def matched(self):
        if self.recA and self.recB:
            return True
        return False

    def has_somatic(self):
        ''' return True if either call is somatic '''
        if not self.matched():
            if self.recA.INFO.get('SS') == 'Somatic':
                return True
            return False

        ss = [self.recA.INFO.get('SS'), self.recB.INFO.get('SS')]
        if 'Somatic' in ss:
            return True
        return False

    def both_somatic(self):
        ''' return True if both calls are somatic '''
        if not self.matched():
            return False

        ss = [self.recA.INFO.get('SS'), self.recB.INFO.get('SS')]
        if 'Germline' in ss and 'Somatic' in ss:
            return False
        return True

    def has_pass(self):
        ''' return True if either filter is PASS '''
        if not self.matched():
            if not self.recA.FILTER:
                return True
            return False

        if not self.recA.FILTER or not self.recB.FILTER:
            return True
        return False

    def both_pass(self):
        ''' return True if both filters are PASS '''
        if not self.matched():
            return False

        if not self.recA.FILTER and not self.recB.FILTER:
            return True
        return False


class SNV (Variant):
    ''' single nucleotide variant subclass '''
    def type(self):
        if self.is_transition:
            return 'transition'
        else:
            return 'transversion'

    def score(self):
        if self.matched():
            return 1.0
        return 0.0


class INDEL (Variant):
    ''' short insertion/deletion subclass '''
    def type(self):
        pass
    def score(self, density=False):
        return self.interval_score(density)

class SV (Variant):
    ''' structural variant subclass '''
    def score(self, density=False):
        return self.interval_score(density)

class CNV (Variant):
    ''' copy number variant subclass '''
    def score(self, density=False):
        return self.interval_score(density)

## functions ##

def get_overlap_coords(iv_a, iv_b):
    ''' return start and end coordinates of overlap between iv_a and iv_b
        return 0,0 if no overlap
        iv_a or iv_b [0] is start and [1] is end '''
    if min(iv_a[1], iv_b[1]) - max(iv_a[0], iv_b[0]) > 0: # is there overlap?
        return max(iv_a[0], iv_b[0]), min(iv_a[1], iv_b[1])
    return 0,0

def interval_density(iv, window):
    ''' returns cumulative density of a window over a rescaled interval (iv)'''

    # set left = 0 and rescale to (-3,3) (adjusting will change how quickly the score function drops off at extremes)
    norm_scale = [-3.0, 3.0]
    iv_scale = [0.0, float(iv[1]-iv[0])]

    rel_w_start = window[0] - iv[0]
    rel_w_end = window[1] - iv[0]

    w_s  = interp(rel_w_start, iv_scale, norm_scale)
    w_e  = interp(rel_w_end, iv_scale, norm_scale)

    return norm_density(w_s, w_e)

def norm_density(start,end):
    ''' return cumulative denisty of normal distribution (from start to end)'''
    start = float(start)
    end = float(end)

    rv = norm()
    d_left = rv.cdf(start)
    d_right = 1.0 - rv.cdf(end)

    return 1.0 - d_left - d_right

def vcfVariantMatch(recA, recB):
    ''' return True if SNV/INDEL/SV/CNV intervals match given critera for each variant type '''

    # SNVs have to have the same position, ref allele, and alt allele
    if recA.is_snp and recB.is_snp:
        if recA.POS == recB.POS and recA.REF == recB.REF and recA.ALT == recB.ALT:
            return True

    # indels have to be within w_indel of each other, have same ref and alt alleles
    if recA.is_indel and recB.is_indel:
        if recA.REF == recB.REF and recA.ALT == recB.ALT:
            return True

    # SVs have to be within w_sv of each other, pass vcfIntervalMatch
    if recA.is_sv and recB.is_sv and recA.INFO.get('SVTYPE') == recB.INFO.get('SVTYPE') == 'BND': 
        orientA = orientSV(str(recA.ALT[0]))
        orientB = orientSV(str(recB.ALT[0]))
        if orientA == orientB:
            return True

    return False 

def compareVCFs(h_vcfA, h_vcfB, w_indel=50, w_sv=1000, mask=None):
    ''' does most of the work - unidirectional comparison vcfA --> vcfB
        h_vcfA and h_vcfB are pyvcf handles (vcf.Reader) '''

    cmp = Comparison()

    for recA in h_vcfA:
        if mask:
            pass # FIXME not implemented

        match = False
        type = None
        variant = None
        w = 0

        if recA.is_snp:
            type = 'SNV'
            variant = SNV(recA, None)

        elif recA.is_indel:
            type = 'INDEL'
            variant = INDEL(recA, None)
            w = w_indel

        elif recA.is_sv and recA.INFO.get('SVTYPE') == 'BND':
            type = 'SV'
            variant = SV(recA, None)
            w = w_sv

        if type: # only compare known variant types
            for recB in h_vcfB.fetch(recA.CHROM, recA.start-w, recA.end+w):
                if vcfVariantMatch(recA, recB):
                    if match: # handle one-to-many matches
                        variant.altmatch.append(recB)
                    else:
                        variant.recB = recB
                        match = True
            cmp.vartype[type].append(variant)

    return cmp

def orientSV(alt):
    '''
    REF   ALT    Meaning
    s     t[p[   piece extending to the right of p is joined after t
    s     t]p]   reverse comp piece extending left of p is joined after t
    s     ]p]t   piece extending to the left of p is joined before t
    s     [p[t   reverse comp piece extending right of p is joined before t
    '''
    orient = alt # return info line by default

    if search('^[A-Z]\[',alt):
        orient = 'right_of_p_after_t'

    elif search('^[A-Z]\]',alt):
        orient = 'left_of_p_after_t'

    elif search('^\]',alt):
        orient = 'left_of_p_before_t'

    elif search('^\[',alt):
        orient = 'right_of_p_before_t'

    return orient

def openVCFs(vcf_list):
    ''' return list of vcf file handles '''
    vcf_handles = []

    for vcf_file in vcf_list:
        try:
            vcf_handles.append(vcf.Reader(filename=vcf_file,compressed=True))
        except IOError as e:
            sys.stderr.write(str(e) + ' -- is this an indexed tabix file?\n')
            sys.exit()

    return vcf_handles

def parseVCFs(vcf_list, maskfile=None):
    ''' handle the list of vcf files
        handle errors '''

    vcf_handles = openVCFs(vcf_list) 

    # give some warnings if the comparison is between different samples
    try:
        if len(vcf_handles[0].metadata['SAMPLE']) != len(vcf_handles[1].metadata['SAMPLE']):
            sys.write.stderr("warning: vcfs have different numbers of samples")
        else:
            for i in range(len(vcf_handles[0].metadata['SAMPLE'])):
                if vcf_handles[0].metadata['SAMPLE'][i]['Individual'] != vcf_handles[1].metadata['SAMPLE'][i]['Individual']:
                    sys.stderr.write("warning: comparing two different individuals: " +
                                     vcf_handles[0].metadata['SAMPLE'][i]['Individual'] + " " +
                                     vcf_handles[1].metadata['SAMPLE'][i]['Individual'] + "\n")
                    break
    except KeyError as e:
        sys.stderr.write("Invalid header, make sure vcf has for following metadata: " + str(e) + "\n")
        sys.exit()

    # compare VCFs
    try:
        sys.stderr.write(vcf_list[0] + " --> " + vcf_list[1] + "\n")
        result01 = compareVCFs(vcf_handles[0], vcf_handles[1])
        result01.summary()

        vcf_handles = openVCFs(vcf_list) 

        sys.stderr.write(vcf_list[1] + " --> " + vcf_list[0] + "\n")
        result10 = compareVCFs(vcf_handles[1], vcf_handles[0])
        result10.summary()

    except ValueError as e:
        sys.stderr.write('error while comparing vcfs: ' + str(e) + '\n')

def main(args):
    parseVCFs(args.vcf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compares two sorted VCF files and (optionally) masks regions.')
    parser.add_argument(metavar='<vcf_file>', dest='vcf', nargs=2, help='tabix-indexed files in VCF format')
    parser.add_argument('-m', '--mask', dest='maskfile', default=None, help='BED file of masked intervals') 
    args = parser.parse_args()
    main(args)
