#!/usr/bin/env python

import vcf
import sys
import traceback
import gzip

'''
validatevcf.py
Contact: Adam Ewing (ewingad@soe.ucsc.edu)

Attempts to validate a VCF that may contain SNVs, INDELs, and SVs
The main test is whether the VCF can be successfully parsed by PyVCF,
a few specific fields are check for as well (VCFs should define somatic,
imprecise breakends should have confidence intervals)

'''


class ValidationError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def context(vcf, recnum):
    vcf_h = None
    ctbuf = 5
    if vcf.endswith('.gz'):
        vcf_h = gzip.open(vcf, 'rb')
    else:
        vcf_h = open(vcf, 'r')

    n = 0
    for line in vcf_h:
        if not line.startswith('#'):
            n += 1
            if n >= recnum-ctbuf and n <= recnum+ctbuf:
                print n,':',line.strip()
                found = True
            if n > recnum + ctbuf:
                vcf_h.close()
                return 
 
if len(sys.argv) == 2:
    assert sys.argv[1].endswith('.vcf') or sys.argv[1].endswith('.vcf.gz')

    vcfin = vcf.Reader(filename=sys.argv[1])

    recnum = 0
    try:
        use_info_somatic = False
        use_info_ss = False
        use_info_loh = False
        use_fmt_ss = False
        indel_count = 0
        snv_count = 0
        sv_count = 0

        for rec in vcfin:
            ''' try to detemine how germline vs. somatic is specified '''
            recnum += 1
            if rec.FILTER == 'GERMLINE' or rec.FILTER == 'SOMATIC':
                raise ValidationError('GERMLINE or SOMATIC in FILTER field: invalid VCF')
            if rec.INFO.get('SOMATIC'):
                use_info_somatic=True
            if str(rec.INFO.get('SS')).upper() == 'SOMATIC':
                use_info_ss=True
            if str(rec.INFO.get('SS')).upper() == 'LOH':
                use_info_loh=True

            if rec.is_snp:
                snv_count += 1
            if rec.is_indel:
                indel_count += 1
            if rec.is_sv:
                sv_count += 1
                if rec.INFO.get('IMPRECISE'):
                    if not (rec.INFO.get('CIPOS')):
                        raise ValidationError('Imprecise SV record without CIPOS: ' + str(rec))

                    if rec.INFO.get('END') and not (rec.INFO.get('CIEND')):
                        raise ValidationError('Imprecise SV record using END without CIEND: ' + str(rec))

            for call in rec.samples:
                data = call.data
                if 'SS' in data._fields:
                    use_fmt_ss = True

        if not (use_info_somatic or use_info_ss or use_fmt_ss):
            raise ValidationError('No somatic records found using INFO or FORMAT fields')

        print "Validation complete."
        print "Total records:",recnum
        print "-"*60
        print "uses INFO/SOMATIC:", use_info_somatic
        print "uses INFO/SS=Somatic/Germline:", use_info_ss
        print "uses INFO/SS=LOH:",use_info_loh
        print "uses FORMAT/SS=0,1,2,...:", use_fmt_ss
        print "-"*60
        print "SNV count:",snv_count
        print "INDEL count:",indel_count
        print "SV count:",sv_count
        print "-"*60

    except:
        raise
        recnum += 1
        print "parse error in VCF on line",recnum
        print '-'*60
        traceback.print_exc(file=sys.stdout)
        print '-'*60
        print "context:"
        context(sys.argv[1], recnum)

else:
    print "usage:",sys.argv[0],"<vcf or vcf.gz>"
