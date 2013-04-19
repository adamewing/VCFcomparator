#!/usr/bin/env python

''' provides a means to query a tabix-indexed VCF file using another VCF and/or various parameters '''

import argparse
import sys
import vcf
import os
import pysam



def is_somatic(rec):
    if str(rec.INFO.get('SS')).upper() in ['SOMATIC', '2']:
        return True

    if rec.INFO.get('SOMATIC'):
        if str(rec.INFO.get('SS')).upper() == 'LOH':
            return False
        return True

    if somatic_in_format(rec):
        return True

    return False

def somatic_in_format(rec):
    SS = []
    for sample in rec.samples:
        calldata = sample.data
        if 'SS' in calldata._fields:
            SS.append(calldata.SS)

    if '2' in SS or 2 in SS:
        return True
    return False

def vcfVariantMatch(recA, recB):
    ''' return True if SNVs/INDELs match given critera for each variant type '''

    if recA.is_snp and recB.is_snp:
        if recA.POS == recB.POS and recA.REF == recB.REF and recA.ALT == recB.ALT:
            return True

    if recA.is_indel and recB.is_indel:
        if recA.POS == recB.POS and recA.REF == recB.REF and recA.ALT == recB.ALT:
            return True

    return False

def main(args):
    if None == args.queryvcf == args.querybed == args.excludevcf and not args.passonly and not args.failonly and not args.maskbed:
        sys.exit("nothing to do!")

    if True == args.passonly == args.failonly:
        sys.exit("called with both -p/--passonly and -f/--failonly, null result")

    assert args.vcf[0].endswith('.vcf') or args.vcf[0].endswith('.vcf.gz') 
    vcfin = vcf.Reader(filename=args.vcf[0])
    vcfout = vcf.Writer(sys.stdout, vcfin)

    # initalize query BED if specified
    qbed = None
    if args.querybed is not None:
        assert args.querybed.endswith('.gz') and os.path.exists(args.querybed + '.tbi')
        qbed = pysam.Tabixfile(args.querybed) 

    # initalize mask BED if specified
    mbed = None
    if args.maskbed is not None:
        assert args.maskbed.endswith('.gz') and os.path.exists(args.maskbed + '.tbi')
        mbed = pysam.Tabixfile(args.maskbed) 

    # initalize query VCF if specified
    qvcf = None
    if args.queryvcf is not None:
        assert args.queryvcf.endswith('vcf.gz') and os.path.exists(args.queryvcf + '.tbi')
        qvcf = vcf.Reader(filename=args.queryvcf)

    # initialize exclusion VCF if specified
    xvcf = None
    if args.excludevcf is not None:
        assert args.excludevcf.endswith('vcf.gz') and os.path.exists(args.excludevcf + '.tbi')
        xvcf = vcf.Reader(filename=args.excludevcf)

    for rec in vcfin:
        output = True
        query_rec = None
        if qvcf is not None:
            output = False
            try:
                for match in qvcf.fetch(rec.CHROM, rec.start, rec.end):
                    if vcfVariantMatch(rec, match):
                        output = True
                        query_rec = match
            except:
                pass

        if qbed is not None:
            output = False
            if rec.CHROM in qbed.contigs:
                if len(list(qbed.fetch(rec.CHROM, rec.start, rec.end))) > 0:
                    output = True

        if mbed is not None:
            if rec.CHROM in mbed.contigs:
                if len(list(mbed.fetch(rec.CHROM, rec.start, rec.end))) > 0:
                    output = False

        if xvcf is not None:
            try:
                for match in xvcf.fetch(rec.CHROM, rec.start, rec.end):
                    if vcfVariantMatch(rec, match):
                        output = False 
            except:
                pass

        if args.passonly and rec.FILTER:
            output = False

        if args.failonly and not rec.FILTER:
            output = False

        if args.somaticonly and not is_somatic(rec):
            output = False

        if args.germlineonly and is_somatic(rec):
            output = False

        if output:
            if args.vcf is not None and args.switch_report:
                vcfout.write_record(query_rec)
            else:
                vcfout.write_record(rec)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='query a tabix-indexed VCF file using another VCF and/or various parameters')
    parser.add_argument(metavar='<vcf_file>', dest='vcf', nargs=1, help='tabix-indexed files in VCF format')
    parser.add_argument('-v', '--vcf', dest='queryvcf', default=None, help='query will return only results matching this (tabix-indexed) VCF file')
    parser.add_argument('-x', '--exclude', dest='excludevcf', default=None, help='query will not return any results matching this (tabix-indexed) VCF file')
    parser.add_argument('-b', '--bed', dest='querybed', default=None, help='tabix-indexed BED file (chrom, start, end, ...), query will return only results in overlapping regions')
    parser.add_argument('-m', '--mask', dest='maskbed', default=None, help='tabix-indexed BED file (chrom, start, end, ...), query will not return results in overlapping regions')
    parser.add_argument('-p', '--passonly', action='store_true', default=False, help='only return PASS records')
    parser.add_argument('-f', '--failonly', action='store_true', default=False, help='only return non-PASS records')
    parser.add_argument('-s', '--somaticonly', action='store_true', default=False, help='only return somatic records')
    parser.add_argument('-g', '--germlineonly', action='store_true', default=False, help='only return germline records')
    parser.add_argument('-r', '--switch_report', action='store_true', default=False, help='report record in query VCF (-v/--vcf) instead of input VCF')
    args = parser.parse_args()
    main(args)
