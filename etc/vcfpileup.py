#!/usr/bin/env python

import argparse
import pysam
import vcf
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

def pileup(bam,chrom,start,end):
    start = int(start)-1
    end = int(end)
    for pcol in bam.pileup(chrom,start,end):
        if pcol.pos >= start and pcol.pos <= end-1:
            baselist = [] 
            poslist  = []
            quallist = []

            for pread in pcol.pileups:
                base = pread.alignment.seq[pread.qpos]
                if pread.alignment.is_reverse:
                    base = base.lower()

                poslist.append(pread.qpos)
                quallist.append(pread.alignment.mapq)
                baselist.append(base)

            print ' '.join((basename(bam.filename), chrom, str(pcol.pos+1), ''.join(baselist), ','.join(map(str, quallist)), ','.join(map(str, poslist)) ))

def main(args):
    h_vcf = vcf.Reader(filename=args.vcffile)

    vtype = None
    if args.vtype is not None:
        assert args.vtype in ('SNV', 'INDEL', 'SV')
        vtype = args.vtype

    h_bams = []
    for bamfile in args.bams:
        h_bams.append(pysam.Samfile(bamfile, 'rb'))

    for rec in h_vcf:
        output = True
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
            print "\n" + str(rec), rec.FILTER
            for bam in h_bams:
                pileup(bam, rec.CHROM, rec.POS-int(args.context), rec.POS+int(args.context))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Grab pileup columns from any number of BAM files corresponding to locations in a VCF file')
    parser.add_argument(metavar='<bam files>', dest='bams', nargs='+', help='BAM files')
    parser.add_argument('-v', '--vcf', dest='vcffile', required=True, help='VCF file')
    parser.add_argument('-c', '--context', dest='context', default=0, help='bases of context on either side of VCF entry')
    parser.add_argument('-t', '--vtype', dest='vtype', default=None, help='only include variants of vtype where vtype is SNV, INDEL, or SV')
    parser.add_argument('-p', '--passonly', action='store_true', default=False, help='only return PASS records')
    parser.add_argument('-f', '--failonly', action='store_true', default=False, help='only return non-PASS records')
    parser.add_argument('-s', '--somaticonly', action='store_true', default=False, help='only return somatic records')
    parser.add_argument('-g', '--germlineonly', action='store_true', default=False, help='only return germline records')
    args = parser.parse_args()
    main(args)
