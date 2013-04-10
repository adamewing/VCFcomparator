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

        for rec in vcfin:
            recnum += 1
            if rec.FILTER == 'GERMLINE' or rec.FILTER == 'SOMATIC':
                use_filter_somatic = True
            if rec.INFO.get('SOMATIC'):
                use_info_somatic=True
            if rec.INFO.get('SS'):
                use_info_ss=True

            for call in rec.samples:
                data = call.data
                if 'SS' in data._fields:
                    use_fmt_ss = True

        print "use_info_somatic:", use_info_somatic
        print "use_info_ss", use_info_ss
        print "use_fmt_ss", use_fmt_ss
        print "use_filter_somatic", use_filter_somatic

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
