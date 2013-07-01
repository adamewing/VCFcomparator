#!/usr/bin/env python

import argparse
import pysam
import vcf
import sys
from os.path import basename


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

def basecount(bam,chrom,pos):
    start = int(pos)-1
    end = int(pos)
    baselist = [] 

    for pcol in bam.pileup(chrom,start,end):
        if pcol.pos >= start and pcol.pos <= end-1:
            for pread in pcol.pileups:
                base = pread.alignment.seq[pread.qpos]
                baselist.append(base)

    bcount = {}
    for b in baselist:
        if b in bcount:
            bcount[b] += 1
        else:
            bcount[b] = 1

    return bcount

def main(args):
    invcf  = vcf.Reader(filename=args.vcffile)
    outvcf = vcf.Writer(sys.stdout, invcf)

    vtype = None
    if args.vtype is not None:
        assert args.vtype in ('SNV', 'INDEL', 'SV')
        vtype = args.vtype

    bam = pysam.Samfile(args.bamfile, 'rb')

    for rec in invcf:
        output = True
        bc = basecount(bam, rec.CHROM, rec.POS)
        for alt in rec.ALT:
            if alt in bc.keys():
                if bc[str(alt)] < int(args.minreads):
                    output = False
            else:
                output = False

        if vtype == 'SNV' and (not rec.is_snp or (rec.is_snp and rec.INFO.get('VT') == 'LOH')):
            output = False

        if vtype == 'INDEL' and not rec.is_indel:
            output = False

        if vtype == 'SV' and not rec.is_sv:
            output = False

        if args.passonly and rec.FILTER:
            output = False

        if args.failonly and not rec.FILTER:
            output = False

        if args.somaticonly and not is_somatic(rec):
            output = False

        if args.germlineonly and is_somatic(rec):
            output = False

        if output:
            outvcf.write_record(rec)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Grab pileup columns from any number of BAM files corresponding to locations in a VCF file')
    parser.add_argument('-v', '--vcf', dest='vcffile', required=True, help='VCF file')
    parser.add_argument('-b', '--bam', dest='bamfile', required=True, help='BAM file')
    parser.add_argument('-m', '--minreads', dest='minreads', default=2, help='minimum reads supporting alt (default=2)')
    parser.add_argument('-c', '--context', dest='context', default=0, help='bases of context on either side of VCF entry')
    parser.add_argument('-t', '--vtype', dest='vtype', default=None, help='only include variants of vtype where vtype is SNV, INDEL, or SV')
    parser.add_argument('-p', '--passonly', action='store_true', default=False, help='only return PASS records')
    parser.add_argument('-f', '--failonly', action='store_true', default=False, help='only return non-PASS records')
    parser.add_argument('-s', '--somaticonly', action='store_true', default=False, help='only return somatic records')
    parser.add_argument('-g', '--germlineonly', action='store_true', default=False, help='only return germline records')
    args = parser.parse_args()
    main(args)
