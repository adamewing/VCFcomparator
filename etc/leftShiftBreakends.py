#!/usr/bin/env python

'''
leftShiftBreakend.py

Author: Adam Ewing (ewingad@soe.ucsc.edu)

Implements left-shifting of precise breakends as described in the VCF 4.1 spec:
http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41

'''


import pysam
import vcf
import sys
import string
import argparse


def rc(seq):
    ''' reverse complement '''
    seq = seq.translate(string.maketrans("atgcATGC","tacgTACG"))
    return seq[::-1] # reverse


def fetch_bnd_seq(rec, bnd, ref, leftwindow, rightwindow, verbose=False):
    ''' get sequence for breakend +/- window'''
    bnd_seq_right  = bnd.orientation # get breakend seq to right?
    mate_seq_right = bnd.remoteOrientation # get mate seq to right?

    mate_flip = bnd_seq_right == mate_seq_right

    bnd_start = bnd_end = mate_start = mate_end = 0

    # breakend 5' end is connected to genome/3' faces mate (t[p[ or t]p])
    if not bnd_seq_right:
        bnd_end   = rec.POS
        bnd_start = bnd_end - leftwindow

    # mate 5' end is connected to genome/3' faces mate (t]p] or ]p]t)
    if not mate_seq_right:
        mate_end   = bnd.pos
        mate_start = mate_end - leftwindow
        if mate_flip:
            mate_start = mate_end - rightwindow

    # breakend 3' end is connected to genome/5' faces mate (]p]t or [p[t)
    if bnd_seq_right:
        bnd_start = rec.POS
        bnd_end   = bnd_start + rightwindow

    # mate 3' end is connected to genome/5' faces bnd (t[p[ or [p[t)
    if mate_seq_right:
        mate_start = bnd.pos
        mate_end   = mate_start + rightwindow
        if mate_flip:
            mate_end = mate_start + leftwindow

    assert bnd_start < bnd_end and mate_start < mate_end
    
    bnd_seq  = ref.fetch(rec.CHROM, bnd_start, bnd_end)
    mate_seq = ref.fetch(bnd.chr, mate_start, mate_end)

    if mate_flip:
        mate_seq = rc(mate_seq)

    # output bnd and mate in correct order
    if bnd_seq_right:
        if verbose:
            sys.stderr.write(mate_seq + "|" + bnd_seq + "\n")
        return mate_seq + bnd_seq 
    else: # bnd is left
        if verbose:
            sys.stderr.write(bnd_seq + "|" + mate_seq + "\n")
        return bnd_seq + mate_seq 

def step_left(rec,bnd):
    ''' shift the breakend left on both sides '''
    bnd_seq_right  = bnd.orientation # get breakend seq to right?
    mate_seq_right = bnd.remoteOrientation # get mate seq to right?

    mate_flip = bnd_seq_right == mate_seq_right

    if mate_flip: # shift mate right if reverse complementing
        rec.POS -= 1
        bnd.pos += 1
    else:
        rec.POS -= 1
        bnd.pos -= 1

    return rec,bnd

def step_right(rec,bnd):
    ''' shift the breakend right on both sides '''
    bnd_seq_right  = bnd.orientation # get breakend seq to right?
    mate_seq_right = bnd.remoteOrientation # get mate seq to right?

    mate_flip = bnd_seq_right == mate_seq_right

    if mate_flip: # shift mate right if reverse complementing
        rec.POS += 1
        bnd.pos -= 1
    else:
        rec.POS += 1
        bnd.pos += 1

    return rec,bnd

def shift_bnd(rec, ref, verbose=False):
    ''' shifts precise breakend record (rec) to left if possible 
        ref is a pysam.Fastafile handle to a ref genome '''
    bnd = rec.ALT[0]
    # these could pose a problem... (initial window sizes)
    l_win = 50
    r_win = 50

    if verbose:
        sys.stderr.write(" ".join(("Original:",str(rec),str(bnd),"\n")))
    orig = fetch_bnd_seq(rec, bnd, ref, l_win, r_win, verbose)

    # keep shifting left until the junctions are no longer equivalent
    while True: 
        l_win -= 1
        r_win += 1

        rec, bnd = step_left(rec, bnd)
        if orig != fetch_bnd_seq(rec, bnd, ref, l_win, r_win, verbose):
            break

    l_win += 1
    r_win -= 1
    rec, bnd = step_right(rec, bnd)
    assert orig == fetch_bnd_seq(rec, bnd, ref, l_win, r_win, verbose)

    if verbose:
        sys.stderr.write(" ".join(("Shifted:",str(rec),str(bnd),"\n\n"))) 

    return rec

def main(args): 
    ''' handle parameters, catch errors '''

    vcf_in  = None
    vcf_out = None

    if args.vcf_infile[0].endswith('.gz'):
        vcf_in = vcf.Reader(filename=args.vcf_infile[0], compressed=True)
    else:
        vcf_in = vcf.Reader(filename=args.vcf_infile[0])

    if args.vcf_outfile:
        vcf_out = vcf.Writer(file(args.vcf_outfile, 'w'), template=vcf_in)
    else:
        vcf_out = vcf.Writer(sys.stdout, template=vcf_in)

    assert vcf_in and vcf_out

    ref = pysam.Fastafile(args.ref_fasta)
    for rec in vcf_in:
        if rec.is_sv:
            if 'IMPRECISE' not in rec.INFO:
                rec = shift_bnd(rec, ref, args.v)

        vcf_out.write_record(rec)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Shifts ambiguous precise breakends to their leftmost position')
    parser.add_argument(metavar='<vcf_file>', dest='vcf_infile', nargs=1, help='VCF file')
    parser.add_argument('-r', dest='ref_fasta', required=True, help='reference genome, .fasta indexed with samtools faidx')
    parser.add_argument('-o', dest='vcf_outfile', default=None, help='output VCF (default STDOUT)')
    parser.add_argument('-v', action='store_true', default=False, help='verbose (for debugging)')
    args = parser.parse_args()
    main(args)
