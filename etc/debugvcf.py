#!/usr/bin/env python

import vcf
import sys
import traceback
import gzip

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
        use_fmt_ss = False
        use_filter_somatic = False
        indel_count = 0
        snv_count = 0
        sv_count = 0

        for rec in vcfin:
            recnum += 1
            if rec.FILTER == 'GERMLINE' or rec.FILTER == 'SOMATIC':
                use_filter_somatic = True
            if rec.INFO.get('SOMATIC'):
                use_info_somatic=True
            if rec.INFO.get('SS'):
                use_info_ss=True

            assert not (rec.is_snp and rec.is_indel and rec.is_sv)
            if rec.is_snp:
                snv_count += 1
            if rec.is_indel:
                indel_count += 1
            if rec.is_sv:
                sv_count += 1

            for call in rec.samples:
                data = call.data
                if 'SS' in data._fields:
                    use_fmt_ss = True

        print "Total records:",recnum
        print "-"*60
        print "uses INFO/SOMATIC:", use_info_somatic
        print "uses INFO/SS=Somatic/Germline:", use_info_ss
        print "uses FORMAT/SS=0,1,2,...:", use_fmt_ss
        print "uses FILTER/SOMATIC (please correct if true, filter should be PASS or filter name):", use_filter_somatic
        print "-"*60
        print "SNV count:",snv_count
        print "INDEL count:",indel_count
        print "SV count:",sv_count
        print "-"*60

    except:
        recnum += 1
        print "parse error in VCF on line",recnum
        print '-'*60
        traceback.print_exc(file=sys.stdout)
        print '-'*60
        print "context:"
        context(sys.argv[1], recnum)

else:
    print "usage:",sys.argv[0],"<vcf or vcf.gz>"
