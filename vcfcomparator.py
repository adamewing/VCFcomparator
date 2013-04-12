#!/usr/bin/env python

'''
VCF comparator: compares mutation calls in VCF formatted files
Distributed under MIT license, see LICENSE.txt
Contact: Adam Ewing (ewingad@soe.ucsc.edu)
'''

import vcf
import pysam
import argparse
import itertools
import sys
import time
import re
import os
import pp

## classes ##

class Comparison:
    ''' stores the result of a one-way comparison vcfA --> vcfB
        imprements functions to report things about the comparison '''
    def __init__(self):
        self.vartype = {}
        self.vartype['SNV']   = []
        self.vartype['INDEL'] = []
        #self.vartype['CNV']   = []
        self.vartype['SV']    = []

    def matched(self, vtype):
        ''' count number of matched variants of <vtype> that PASS '''
        return reduce(lambda x, y: x+(y.matched() and y.both_pass()), self.vartype[vtype], 0)

    def matched_any(self, vtype):
        ''' count number of matched variants of <vtype> that PASS '''
        return reduce(lambda x, y: x+y.matched(), self.vartype[vtype], 0)

    def unmatched(self, vtype):
        ''' count total number of unmatched variants of <vtype> that PASS '''
        return reduce(lambda x, y: x+(not y.matched() and y.has_pass()), self.vartype[vtype], 0)

    def unmatched_any(self, vtype):
        ''' count total number of unmatched variants of <vtype> that PASS '''
        return reduce(lambda x, y: x+(not y.matched()), self.vartype[vtype], 0)

    def unmatched_germline(self, vtype):
        ''' count number of unmatched germline variants of <vtype> that PASS'''
        return reduce(lambda x, y: x+(not y.matched() and y.has_pass() and y.has_germline()), self.vartype[vtype], 0)

    def unmatched_somatic(self,vtype):
        ''' count number of unmatched somatic variants of <vtype> that PASS'''
        return reduce(lambda x, y: x+(not y.matched() and y.has_pass() and y.has_somatic()), self.vartype[vtype], 0)

    def altmatched(self,vtype):
        ''' count number of alternate matches '''
        return reduce(lambda x, y: x+len(y.altmatch), self.vartype[vtype], 0)

    def count_agree_somatic(self,vtype):
        ''' count number of matches that pass and agree on somatic status given variant type '''
        return reduce(lambda x, y: x+(y.matched() and y.both_somatic() and y.both_pass()), self.vartype[vtype], 0)

    def count_agree_germline(self,vtype):
        ''' count number of matches that pass and agree on germline status given variant type '''
        return reduce(lambda x, y: x+(y.matched() and y.both_germline() and y.both_pass()), self.vartype[vtype], 0)

    def count_disagree_somatic(self,vtype):
        ''' count number of matches that pass and disagree on somatic status given variant type '''
        return reduce(lambda x, y: x+(y.matched() and y.both_pass() and y.has_somatic() and not y.both_somatic()), self.vartype[vtype], 0)

    def count_agree_pass(self,vtype): # not useful?
        ''' count number of matches that agree on passing filters given variant type'''
        return reduce(lambda x, y: x+y.both_pass(), self.vartype[vtype], 0)

    def count_agree_fail(self,vtype): # not useful?
        ''' count number of matches that agree on not passing filters given variant type '''
        return reduce(lambda x, y: x+(y.has_pass()==False), self.vartype[vtype], 0)

    def count_disagree_pass(self,vtype): 
       ''' count number of matches where one call passed and the other did not given variant type '''
       return reduce(lambda x, y: x+(y.matched() and y.has_pass() and not y.both_pass()), self.vartype[vtype], 0)

    def sum_scores(self,vtype):
        ''' return sum of all scores for variant type '''
        return reduce(lambda x, y: x+y.score(), self.vartype[vtype], 0.0)

class Variant:
    ''' base class for variant types 
        vcf_recA and vcf_recB are vcf._Record objects or None '''
    def __init__(self, vcf_recA, vcf_recB):
        self.recA = vcf_recA
        self.recB = vcf_recB

        # if there is more than one match, the rest are stored here
        self.altmatch = []

        self.truth = None # will be a VCF record if not None

    def __str__(self):
        return str(self.recA) + "\t" + str(self.recB)

    def set_left(self, vcf_recB):
        ''' sets record B only if it is not already set '''
        if not self.recB:
            self.recB = vcf_recB
            return True
        return False

    def interval_score(self):
        ''' scoring function for intervals, based on amount of overlap '''
        # get pos,end +/- confidence intervals if present
        iv_a = get_conf_interval(self.recA)
        iv_b = get_conf_interval(self.recB)

        ol_coords = get_overlap_coords(iv_a, iv_b)
        ol_width = ol_coords[1] - ol_coords[0]
        assert ol_width > 0

        len_a = iv_a[1] - iv_a[0]
        len_b = iv_b[1] - iv_b[0]
        assert len_a > 0
        assert len_b > 0

        s = float(2*ol_width)/float(len_a+len_b)

        return s

    def somatic_in_format(self, rec):
        SS = []
        for sample in rec.samples:
            calldata = sample.data
            if 'SS' in calldata._fields:
                SS.append(calldata.SS)

        if '2' in SS or 2 in SS:
            return True
        return False

    def matched(self):
        if self.recA and self.recB:
            return True
        return False

    def truthed(self):
        if self.recT:
            return True
        return False

    def recA_somatic(self):
        ''' return True if recA is somatic '''
        if str(self.recA.INFO.get('SS')).upper() in ['SOMATIC', '2']:
            return True

        if self.recA.INFO.get('SOMATIC'):
            return True

        if self.somatic_in_format(self.recA):
            return True

        return False

    def recB_somatic(self):
        ''' return True if recB is somatic '''
        if not self.matched():
            return False

        if str(self.recB.INFO.get('SS')).upper() in ['SOMATIC', '2']:
            return True

        if self.recB.INFO.get('SOMATIC'):
            return True

        if self.somatic_in_format(self.recB):
            return True

        return False

    def has_somatic(self):
        ''' return True if either call is somatic '''
        return self.recA_somatic() or self.recB_somatic()

    def both_somatic(self):
        ''' return True if both calls are somatic '''
        if not self.matched():
            return False

        return self.recA_somatic() and self.recB_somatic()

    def has_germline(self):
        ''' return True if either call is germline '''
        if not self.matched():
            if not self.has_somatic():
                return True
            return False

        if not self.both_somatic():
            return True

        return False

    def both_germline(self):
        ''' return True if both calls are germline '''
        if not self.matched():
            return False

        if self.has_somatic():
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
    def vtype(self):
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
    def vtype(self):
        pass
    def score(self):
        if self.matched():
            return 1.0
        return 0.0

class SV (Variant):
    ''' structural variant subclass '''
    def score(self):
        if self.matched():
            return self.interval_score()
        return 0.0

class CNV (Variant):
    ''' copy number variant subclass '''
    def score(self):
        if self.matched():
            return self.interval_score()
        return 0.0

class Segment:
    ''' used for segmenting the genome into chunks for threading '''
    def __init__(self, line=None):
        self.chrom  = None
        self.start  = None
        self.end    = None
        self.length = None

        if line is not None:
            self.line = line.strip()
            c = line.strip().split()

            self.chrom  = c[0]
            self.length = int(c[1])
            self.start  = 0
            self.end    = int(c[1])

    def __gt__(self, other):
        return self.length > other.length

    def __str__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end) + " (" + str(self.length) + " bp)"

    def left(self):
        lseg = Segment()
        lseg.chrom = self.chrom
        lseg.start = self.start
        lseg.end = self.start + int(self.length/2)
        lseg.length = lseg.end - lseg.start
        return lseg

    def right(self):
        rseg = Segment()
        rseg.chrom = self.chrom
        rseg.start = self.end - int(self.length/2)
        rseg.end = self.end
        rseg.length = rseg.end - rseg.start
        return rseg

## functions ##

def get_conf_interval(rec, w_indel=0):
    ''' return confidence interval as (start-ci, end+ci), if rec is an indel, w_indel is added to interval'''
    cipos_start = cipos_end = 0
    ciend_start = ciend_end = 0
    
    if 'CIPOS' in rec.INFO:
        if len(rec.INFO.get('CIPOS')) == 2:
            cipos_start, cipos_end = map(abs, rec.INFO.get('CIPOS'))
        else:
            cipos_start = cipos_end = abs(rec.INFO.get('CIPOS')[0])

    if 'CIEND' in rec.INFO:
        if len(rec.INFO.get('CIEND')) == 2:
            ciend_start, ciend_end = map(abs, rec.INFO.get('CIEND'))
        else:
            ciend_start = ciend_end = abs(rec.INFO.get('CIEND')[0])

    # default to using POS+1 as end
    end = rec.POS+1

    # if 'END' is specified, use as end
    if 'END' in rec.INFO:
        try:
            end = rec.INFO.get('END')[0] + ciend_end
        except TypeError:
            end = rec.INFO.get('END') + ciend_end
    else:
        end += cipos_end

    if rec.is_indel:
        cipos_start += w_indel 
        end += w_indel

    return rec.POS-cipos_start, end

def get_overlap_coords(iv_a, iv_b):
    ''' return start and end coordinates of overlap between iv_a and iv_b
        return 0,0 if no overlap
        iv_a or iv_b [0] is start and [1] is end '''
    if min(iv_a[1], iv_b[1]) - max(iv_a[0], iv_b[0]) > 0: # is there overlap?
        return max(iv_a[0], iv_b[0]), min(iv_a[1], iv_b[1])
    return 0,0

def vcfVariantMatch(recA, recB):
    ''' return True if SNV/INDEL/SV/CNV intervals match given critera for each variant type '''

    # SNVs have to have the same position, ref allele, and alt allele
    if recA.is_snp and recB.is_snp:
        if recA.POS == recB.POS and recA.REF == recB.REF and recA.ALT == recB.ALT:
            return True

    # indels must have same ref and alt alleles
    if recA.is_indel and recB.is_indel:
        if recA.REF == recB.REF and recA.ALT == recB.ALT:
            return True

    # SVs have to be within w_sv of each other, pass vcfIntervalMatch
    if recA.is_sv and recB.is_sv and recA.INFO.get('SVTYPE') == recB.INFO.get('SVTYPE') == 'BND': 
        orientA = orientSV(str(recA.ALT[0]))
        orientB = orientSV(str(recB.ALT[0]))
        if orientA == orientB and vcfIntervalMatch(recA, recB):
            return True
    return False 

def vcfIntervalMatch(recA, recB):
    ''' match SV/CNV intervals using POS/END/CIPOS/CIEND '''
    assert recA.INFO.get('SVTYPE') == recB.INFO.get('SVTYPE')

    iv_A = get_conf_interval(recA)
    iv_B = get_conf_interval(recB)

    if sum(get_overlap_coords(iv_A, iv_B)) > 0.0:
        return True
    return False

# copy of file handle for snv iteration and interval fetch
def compareVCFs(h_vcfA, h_interval_vcfB, verbose=False, w_indel=0, w_sv=1000, mask=None, truth=None, chrom=None, fetch_start=0, fetch_end=int(1e9)): 
    ''' does most of the work - unidirectional comparison vcfA --> vcfB
        h_vcfA and h_vcfB are pyvcf handles (vcf.Reader) '''

    h_snv_vcfB = vcf.Reader(filename=h_interval_vcfB.filename, compressed=h_interval_vcfB.filename.endswith('.gz'))

    cmp = Comparison()

    # keep match symmetric by adding B records already seen to altmatch (intervals only)
    used_B_interval = {}

    recnum = 0
    nskip = 0

    # "int(1e9)" is just a value larger than any hg19 chromosome, fetch(chrom,start) not supported
    for recA in h_vcfA.fetch(chrom,fetch_start,fetch_end):
        recnum += 1
        if mask:
            # skip variants on chromosomes not in mask
            if recA.CHROM not in mask.contigs:
                nskip += 1
                continue

            if len(list(mask.fetch(recA.CHROM, recA.POS, recA.POS+1))) > 0:
                nskip += 1
                continue

        if verbose:
            if recnum % 10000 == 0:
                localtime = time.asctime(time.localtime(time.time()))
                sys.stderr.write(str(localtime) + ": " + os.path.basename(h_vcfA.filename) + " vs " 
                                 + os.path.basename(h_interval_vcfB.filename) + ": " + str(recnum) 
                                 + " records compared, pos: " + str(recA.CHROM) + ":"
                                 + str(recA.POS) + " masked: " + str(nskip) + "\n")

        match = False
        vtype = None
        variant = None
        w = 0

        if recA.is_snp:
            vtype = 'SNV'
            variant = SNV(recA, None)

        elif recA.is_indel:
            vtype = 'INDEL'
            variant = INDEL(recA, None)
            w = w_indel

        elif recA.is_sv and recA.INFO.get('SVTYPE') == 'BND':
            vtype = 'SV'
            variant = SV(recA, None)
            w = w_sv

        elif recA.ALT == 'CNV':
            vtype = 'CNV'
            variant = CNV(recA, None)

        if vtype in ('SNV', 'SV', 'INDEL'): # only compare intervals for known variant types
            w_start = recA.start-w
            w_end = recA.end+w
            if w_start < 1:
                w_start = 1

            # try to find a match in the other VCF
            try:
                for recB in h_interval_vcfB.fetch(recA.CHROM, w_start, w_end):
                    if vcfVariantMatch(recA, recB):
                        if match: # handle one-to-many matches
                            variant.altmatch.append(recB)
                        else:
                            assert variant.recB is None

                            # special case for intervals
                            if vtype in ('INDEL','SV','CNV') and sv_uid(recB) in used_B_interval:
                                variant.altmatch.append(recB)

                            elif variant.set_left(recB):
                                used_B_interval[sv_uid(recB)] = recA
                                match = True
            except:
                sys.stderr.write(' '.join(("warning: couldn't fetch from region:", str(recA.CHROM), str(w_start), str(w_end), "\n")))

            # compare to truth if present
            if truth is not None:
                #try:
                for recT in truth.fetch(recA.CHROM, w_start, w_end):
                    if vcfVariantMatch(recA, recT):
                        recA.truth = recT
                #except:
                #    sys.stderr.write(' '.join(("warning: truth VCF missing region:", str(recA.CHROM), str(w_start), str(w_end), "\n")))

            cmp.vartype[vtype].append(variant)

    return cmp

def sv_uid(rec):
    ''' makes a (hopefully) unique id for an SV record '''
    fields = (rec.CHROM,rec.POS,rec.ID,rec.REF,rec.ALT,rec.QUAL,rec.FILTER,rec.INFO)
    return ','.join(map(str,fields))

def orientSV(alt):
    '''
    REF   ALT    Meaning
    s     t[p[   piece extending to the right of p is joined after t
    s     t]p]   reverse comp piece extending left of p is joined after t
    s     ]p]t   piece extending to the left of p is joined before t
    s     [p[t   reverse comp piece extending right of p is joined before t
    '''
    orient = alt # return info line by default

    if re.search('^[A-Z]\[',alt):
        orient = 'right_of_p_after_t'

    elif re.search('^[A-Z]\]',alt):
        orient = 'left_of_p_after_t'

    elif re.search('^\]',alt):
        orient = 'left_of_p_before_t'

    elif re.search('^\[',alt):
        orient = 'right_of_p_before_t'

    return orient

def summary(compAB_list, compBA_list, outfile=None, chrom=None, start=None, end=None):
    ''' given A --> B comparison and B --> A comparison, output summary stats to outfile, or to stdout if outfile=None '''

    out = []
    out.append(' '.join(('chrom','start','end','vtype','A_only_pass_somatic','A_only_pass_germline', 'B_only_pass_somatic',
                         'B_only_pass_germline','match_pass_somatic','match_pass_germline', 'A_pass_disagree_somatic',
                         'B_pass_disagree_somatic')))

    for vtype in compAB_list[0].vartype.keys():
        assert compBA_list[0].vartype.has_key(vtype)

        n_only_A_somatic  = 0 
        n_only_B_somatic  = 0
        n_only_A_germline = 0
        n_only_B_germline = 0
        n_agree_som       = 0
        n_agree_germ      = 0
        n_disagree_som_A  = 0
        n_disagree_som_B  = 0
        n_shared_AB       = 0
        n_shared_BA       = 0
    
        for compAB, compBA in itertools.izip(compAB_list, compBA_list):
            n_shared_AB       += compAB.matched(vtype)
            n_shared_BA       += compBA.matched(vtype)
            n_only_A_somatic  += compAB.unmatched_somatic(vtype)
            n_only_B_somatic  += compBA.unmatched_somatic(vtype)
            n_only_A_germline += compAB.unmatched_germline(vtype)
            n_only_B_germline += compBA.unmatched_germline(vtype)
            n_agree_som       += compAB.count_agree_somatic(vtype)
            n_agree_germ      += compAB.count_agree_germline(vtype)
            n_disagree_som_A  += compAB.count_disagree_somatic(vtype)
            n_disagree_som_B  += compBA.count_disagree_somatic(vtype)

        outstr = map(str, (chrom, start, end, vtype, n_only_A_somatic, n_only_A_germline, n_only_B_somatic, n_only_B_germline, 
                           n_agree_som, n_agree_germ, n_disagree_som_A, n_disagree_som_B))
        out.append(' '.join(outstr))

    if n_shared_AB != n_shared_BA:
        if chrom is not None:
            sys.stderr.write("region " + chrom + ":" + str(start) + "-" + str(end) + ": ")
            sys.stderr.write("warning: overlap was not symmetric (A-->B: " + str(n_shared_AB) + "), (B-->A: " + str(n_shared_BA) + 
                             ") using A-->B\n")

    f = sys.stdout
    if outfile:
        f = open(outfile, 'w')
    f.write('\n'.join(out) + "\n")
    if outfile:
        f.close()

def outputVCF(comparison_list, inVCFhandle, outdir):
    ''' write VCF files for matched and unmatched records, for matched variants, output the record from sample A'''
    ifname = os.path.basename(inVCFhandle.filename)
    assert ifname.endswith('.vcf.gz')

    if outdir is not None:
        if not os.path.exists(outdir):
            sys.stderr.write("creating output directory: " + outdir + "\n")
            os.makedirs(outdir)
        ifname = outdir + '/' + ifname

    ofname_match = re.sub('vcf.gz$', 'matched.vcf', ifname)
    vcfout_match   = vcf.Writer(file(ofname_match, 'w'), inVCFhandle)

    ofname_unmatch = re.sub('vcf.gz$', 'unmatched.vcf', ifname)
    vcfout_unmatch = vcf.Writer(file(ofname_unmatch, 'w'), inVCFhandle)

    match = 0
    unmatch = 0

    for comparison in comparison_list:
        for vtype in comparison.vartype.keys():
            for var in comparison.vartype[vtype]:
                if var.matched():
                    vcfout_match.write_record(var.recA)
                    match += 1
                else:
                    vcfout_unmatch.write_record(var.recA)
                    unmatch +=1

    print ofname_match + ":", match
    print ofname_unmatch + ":", unmatch

    vcfout_match.close()
    vcfout_unmatch.close()

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

def parseVCFs(vcf_list, maskfile=None, truthvcf=None, chrom=None, start=None, end=None, verbose=False):
    ''' handle the list of vcf files and handle errors '''
    assert len(vcf_list) == 2
    vcf_handles = openVCFs(vcf_list) 
    assert len(vcf_handles) == 2

    tabix_mask = None
    if maskfile is not None:
        try:
            tabix_mask = pysam.Tabixfile(maskfile)
        except:
            sys.stderr.write("could not read mask: " + maskfile + "  is it a tabix-indexed bgzipped BED?\n")
            sys.exit()

    tabix_truth = None
    if truthvcf is not None:
        try:
            tabix_truth = openVCFs([truthvcf])[0]
        except:
            sys.stderr.write("could not read mask: " + truthvcf + "  is it a tabix-indexed bgzipped VCF?\n")
            sys.exit()
            
    # compare VCFs
    try:
        sys.stderr.write(vcf_list[0] + " --> " + vcf_list[1] + "\n")
        resultAB = compareVCFs(vcf_handles[0], vcf_handles[1], verbose=verbose, mask=tabix_mask, truth=tabix_truth, chrom=chrom, fetch_start=start, fetch_end=end)

        # reload vcfs to reset iteration
        vcf_handles = openVCFs(vcf_list) 

        sys.stderr.write(vcf_list[1] + " --> " + vcf_list[0] + "\n")
        resultBA = compareVCFs(vcf_handles[1], vcf_handles[0], verbose=verbose, mask=tabix_mask, truth=tabix_truth, chrom=chrom, fetch_start=start, fetch_end=end)

        return resultAB, resultBA, vcf_handles

    except ValueError as e:
        sys.stderr.write('error while comparing vcfs: ' + str(e) + '\n')

def split_genome(chroms, n, minlen=1e6, verbose=False):
    ''' used externally, repeatedly split the largest segment, repect chromosome boundaries '''
    assert n > 0
    segs = []
    with open(chroms, 'r') as chromfile:
        for line in chromfile:
            seg = Segment(line=line)
            if seg.length > minlen:
                segs.append(seg)

    segs.sort()

    # break the longest chunk in half until we have enough
    while n > len(segs):
        bigseg = segs.pop()
        segs.append(bigseg.left())
        segs.append(bigseg.right())
        segs.sort()

    assert n <= len(segs)

    jobs = []
    for i in range(n):
        jobs.append([])

    while len(segs) > 0:
        for j in range(n):
            if len(segs) == 0:
                break
            jobs[j].append(segs.pop())

    if verbose:
        print "-"*60
        print "segmented into",n,"jobs:"
        for i in range(len(jobs)):
            print "job",str(i) + ":",','.join(map(str,jobs[i]))
        print "-"*60

    return jobs

def runList(args, seg_list):
    ''' used by external script to parallelize jobs '''    
    resultsAB = []
    resultsBA = []
    vcf_handles = None

    for seg in seg_list: # Segment
        resultAB, resultBA, vcf_handles = parseVCFs(args.vcf, maskfile=args.maskfile, truthvcf=args.truth, chrom=seg.chrom, start=seg.start, end=seg.end, verbose=args.verbose)
        resultsAB.append(resultAB)
        resultsBA.append(resultBA)
        if args.verbose:
            summary([resultAB], [resultBA], chrom=seg.chrom, start=seg.start, end=seg.end)

    return resultsAB, resultsBA, vcf_handles

def main(args):
    resultAB, resultBA, vcf_handles = parseVCFs(args.vcf, maskfile=args.maskfile, truthvcf=args.truth, chrom=args.chrom, start=int(args.start), end=int(args.end), verbose=args.verbose)
    summary([resultAB], [resultBA], chrom=args.chrom, start=args.start, end=args.end)
    outputVCF([resultAB], vcf_handles[0], args.outdir)
    outputVCF([resultBA], vcf_handles[1], args.outdir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compares two sorted VCF files and (optionally) masks regions.')
    parser.add_argument(metavar='<vcf_file>', dest='vcf', nargs=2, help='tabix-indexed files in VCF format')
    parser.add_argument('-m', '--mask', dest='maskfile', default=None, help='tabix-indexed BED file of masked intervals') 
    parser.add_argument('-o', '--outdir', dest='outdir', default=None, help='directory for output')
    parser.add_argument('-t', '--truth', dest='truth', default=None, help='also compare results to a "truth" VCF (should be sorted and tabix-indexed)')
    parser.add_argument('-c', '--chrom', dest='chrom', default=None, help='limit to one chromosome')
    parser.add_argument('-s', '--start', dest='start', default=0, help='start position')
    parser.add_argument('-e', '--end', dest='end', default=int(1e9), help='end position') 
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='verbose mode for debugging')
    args = parser.parse_args()
    main(args)


