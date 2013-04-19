#!/usr/bin/env python

import sys
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

def is_snp(rec):
    snp = True
    if not rec.is_snp:
        snp = False
    if rec.INFO.get('VT') == 'LOH':
        snp = False

    SS = []
    has_format_ss = False
    for sample in rec.samples:
        calldata = sample.data
        if 'SS' in calldata._fields:
            has_format_ss = True
            SS.append(calldata.SS)

    if has_format_ss:
        if not '1' in SS and not 1 in SS and not '2' in SS and not 2 in SS:
            snp = False

    return snp 

if len(sys.argv) == 3:
    assert sys.argv[2] in ('SNP', 'INDEL')
    vcfin = vcf.Reader(filename=sys.argv[1])
    n_total = 0
    n_pass  = 0
    n_fail  = 0
    n_germ  = 0
    n_som   = 0

    n_som_pass  = 0
    n_som_fail  = 0
    n_germ_pass = 0
    n_germ_fail = 0

    fail_reasons = {}

    for rec in vcfin:
        if (is_snp(rec) and sys.argv[2] == 'SNP') or (rec.is_indel and sys.argv[2] == 'INDEL'): 
            n_total += 1
            passed   = False
            failed   = False
            somatic  = False
            germline = False

            if is_somatic(rec):
                n_som += 1
                somatic = True
            else:
                n_germ += 1
                germline = True

            if rec.FILTER:
                n_fail += 1
                failed = True
                for flag in rec.FILTER:
                    if flag in fail_reasons:
                        fail_reasons[flag] += 1
                    else:
                        fail_reasons[flag] = 1
            else:
                n_pass += 1
                passed = True

            assert germline != somatic
            assert passed != failed

            if somatic and passed:
                n_som_pass += 1
            if somatic and failed:
                n_som_fail += 1
            if germline and passed:
                n_germ_pass += 1
            if germline and failed:
                n_germ_fail += 1

    print "VCF", basename(sys.argv[1])
    print "VTYPE", sys.argv[2]
    print "Total", n_total
    print "Total_Germline",  n_germ
    print "Total_Somatic",   n_som
    print "Total_Passed",    n_pass
    print "Total_Failed",    n_fail
    print "Somatic_Passed",  n_som_pass
    print "Somatic_Failed",  n_som_fail
    print "Germline_Passed", n_germ_pass
    print "Germline_Failed", n_germ_fail
else:
    print "usage:", sys.argv[0], "<VCF> <VTYPE (SNP/INDEL)>"
