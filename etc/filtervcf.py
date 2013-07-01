#!/usr/bin/env python

import vcf
import sys
import argparse

def get_val(a):
    # is it iterable?
    try:
        for val in a:
            if isinstance(val, (int, long, float)):
                return val

    except TypeError:
        if isinstance(a, (int, long, float)):
            return a

class Filter:
    def __init__(self, filtline):
        self.field, self.tag, self.dir, self.value = filtline.strip().split()[:4]
        self.field = self.field.upper()
        self.dir   = self.dir.upper()
        self.value = float(self.value)

        assert self.field in ('INFO', 'FORMAT')
        assert self.dir in ('LTE', 'GT')

        self.fmtname = None
        if self.field == 'FORMAT':
            self.fmtname = filtline.strip().split()[4]

    def is_filtered(self, rec, quiet=False):
        ''' return True if rec is filtered, False otherwise '''
        if self.field == 'INFO':
            if rec.INFO.get(self.tag):
                recval = float(get_val(rec.INFO.get(self.tag)))
                if self.dir == 'LTE':
                    return recval <= self.value
                else:
                    return recval > self.value 
            else:
                if not quiet:
                    sys.stderr.write("warning: INFO tag " + self.tag + " not in record " + str(rec) + "\n")
                return False

        if self.field == 'FORMAT':
            all_samples = []
            for sample in rec.samples:
                all_samples.append(sample.sample)
                if sample.sample == self.fmtname:
                    fmt = sample.data._asdict()
                    if self.tag in fmt and get_val(fmt[self.tag]) is not None:
                        recval = float(get_val(fmt[self.tag]))
                        if self.dir == 'LTE':
                            return recval <= self.value
                        else:
                            return recval > self.value 

                    else:
                        if not quiet:
                            sys.stderr.write("warning: FORMAT tag " + self.tag + " not in record " + str(rec) + "\n")
                    return False

            assert self.fmtname in all_samples


def main(args):
    assert args.vcf[0].endswith('.vcf') or args.vcf[0].endswith('.vcf.gz')
    vcfin  = vcf.Reader(filename=args.vcf[0])
    vcfout = vcf.Writer(sys.stdout, vcfin)

    filters = []
    with open(args.ff, 'r') as ff:
        for line in ff:
            filters.append(Filter(line))

    for rec in vcfin:
        if not rec.FILTER:
            filtered = False
            for filter in filters:
                if not filtered:
                    filtered = filter.is_filtered(rec, quiet=args.quiet)

            if filtered:
                rec.FILTER = ['autofilter']

        vcfout.write_record(rec)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Apply filters to VCF')
    parser.add_argument(metavar='<vcf_file>', dest='vcf', nargs=1, help='file in VCF format')
    parser.add_argument('-f', '--filterfile', dest='ff', required=True, 
                        help='File with filter parameters. Format (whitespace delimited): <field>(INFO or FORMAT) <tag> <direction>(LTE or GT) <value> <sample name>(FORMAT only)')
    parser.add_argument('--quiet', action='store_true', default=False, help='supress warning messages')
    args = parser.parse_args()
    main(args)
