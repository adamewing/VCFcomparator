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
        for rec in vcfin:
            recnum += 1
            for call in rec.samples:
                data = call.data
    except:
        print "parse error in VCF on line",recnum
        print '-'*60
        traceback.print_exc(file=sys.stdout)
        print '-'*60
        print "context:"
        context(sys.argv[1], recnum)

else:
    print "usage:",sys.argv[0],"<vcf or vcf.gz>"
