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

def queryvcf(vcf_fn, silent=False, queryvcf=None, excludevcf=None, querybed=None, maskbed=None, vtype=None, passonly=False, 
             failonly=False, somaticonly=False, germlineonly=False, switch_report=False, outfile=None, similarity=False):

    #sys.stderr.write(','.join(map(str, (vcf_fn, queryvcf, excludevcf, querybed, maskbed, vtype, passonly, failonly, somaticonly, germlineonly, switch_report, outfile, similarity))) + "\n")

    if None == queryvcf == querybed == excludevcf == vtype and not somaticonly and not germlineonly and not passonly and not failonly and not maskbed:
        sys.exit("nothing to do!")

    if True == passonly == failonly:
        sys.exit("called with both -p/--passonly and -f/--failonly, null result")

    if switch_report and queryvcf is None:
        sys.exit("cannot use -r/--switch_report without -v/--vcf")

    assert vcf_fn.endswith('.vcf') or vcf_fn.endswith('.vcf.gz') 
    vcfin = vcf.Reader(filename=vcf_fn)

    assert vcfin is not None

    # initalize query BED if specified
    qbed = None
    if querybed is not None:
        assert querybed.endswith('.gz') and os.path.exists(querybed + '.tbi')
        qbed = pysam.Tabixfile(querybed)

    # initalize mask BED if specified
    mbed = None
    if maskbed is not None:
        assert maskbed.endswith('.gz') and os.path.exists(maskbed + '.tbi')
        mbed = pysam.Tabixfile(maskbed) 

    # initalize query VCF if specified
    qvcf = None
    if queryvcf is not None:
        assert queryvcf.endswith('vcf.gz') and os.path.exists(queryvcf + '.tbi')
        qvcf = vcf.Reader(filename=queryvcf)

    # initialize exclusion VCF if specified
    xvcf = None
    if excludevcf is not None:
        assert excludevcf.endswith('vcf.gz') and os.path.exists(excludevcf + '.tbi')
        xvcf = vcf.Reader(filename=excludevcf)

    n_query = 0
    if similarity and not qvcf:
        sys.exit("cannot compute similarity without -v/--vcf")
    else:
        if qvcf is not None:
            for rec in qvcf:
                n_query += 1

    # variant type
    if vtype is not None:
        assert vtype in ('SNV', 'INDEL', 'SV')

    vcfout = None
    if not silent:
        if switch_report:
            if outfile is not None:
                vcfout = vcf.Writer(file(outfile, 'w'), qvcf)
            else:
                vcfout = vcf.Writer(sys.stdout, qvcf)
        else:
            if outfile is not None:
                vcfout = vcf.Writer(file(outfile, 'w'), vcfin)
            else:
                vcfout = vcf.Writer(sys.stdout, vcfin)

    n_output = 0
    n_total  = 0

    for rec in vcfin:
        n_total += 1
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

        if vtype == 'SNV' and (not rec.is_snp or (rec.is_snp and rec.INFO.get('VT') == 'LOH')):
            output = False

        if vtype == 'INDEL' and not rec.is_indel:
            output = False

        if vtype == 'SV' and not rec.is_sv:
            output = False

        if passonly and rec.FILTER:
            output = False

        if failonly and not rec.FILTER:
            output = False

        if somaticonly and not is_somatic(rec):
            output = False

        if germlineonly and is_somatic(rec):
            output = False

        if output:
            n_output += 1
            if not silent:
                if vcf is not None and switch_report:
                    vcfout.write_record(query_rec)
                else:
                    vcfout.write_record(rec)

    if similarity:
        return n_total, n_query, n_output, (float(n_output) / float(n_query + n_total - n_output))
    else:
        return n_total, n_query, n_output, float(n_output)

def main(args):
    s = queryvcf(args.vcf[0], queryvcf=args.queryvcf, excludevcf=args.excludevcf, querybed=args.querybed, 
                 maskbed=args.maskbed, vtype=args.vtype, passonly=args.passonly, failonly=args.failonly, 
                 somaticonly=args.somaticonly, germlineonly=args.germlineonly, switch_report=args.switch_report,
                 outfile=args.outfile, similarity=args.return_sim)

    if args.return_sim:
        print "similarity:", s
    else:
        print "count:", s 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='query a tabix-indexed VCF file using another VCF and/or various parameters')
    parser.add_argument(metavar='<vcf_file>', dest='vcf', nargs=1, help='tabix-indexed files in VCF format')
    parser.add_argument('-v', '--vcf', dest='queryvcf', default=None, help='query will return only results matching this (tabix-indexed) VCF file')
    parser.add_argument('-x', '--exclude', dest='excludevcf', default=None, help='query will not return any results matching this (tabix-indexed) VCF file')
    parser.add_argument('-b', '--bed', dest='querybed', default=None, help='tabix-indexed BED file (chrom, start, end, ...), query will return only results in overlapping regions')
    parser.add_argument('-m', '--mask', dest='maskbed', default=None, help='tabix-indexed BED file (chrom, start, end, ...), query will not return results in overlapping regions')
    parser.add_argument('-t', '--vtype', dest='vtype', default=None, help='only include variants of vtype where vtype is SNV, INDEL, or SV')
    parser.add_argument('-o', '--outfile', dest='outfile', default=None, help='VCF output to file')
    parser.add_argument('-p', '--passonly', action='store_true', default=False, help='only return PASS records')
    parser.add_argument('-f', '--failonly', action='store_true', default=False, help='only return non-PASS records')
    parser.add_argument('-s', '--somaticonly', action='store_true', default=False, help='only return somatic records')
    parser.add_argument('-g', '--germlineonly', action='store_true', default=False, help='only return germline records')
    parser.add_argument('-r', '--switch_report', action='store_true', default=False, help='report record in query VCF (-v/--vcf required) instead of input VCF')
    parser.add_argument('-i', '--return_sim', action='store_true', default=False, help='return meet/min similarity based on queryvcf (-v/--vcf required)')
    args = parser.parse_args()
    main(args)
