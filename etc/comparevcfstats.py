#!/usr/bin/env python

import argparse
import vcf
import pprint
import numpy as np
import scipy.stats as ss
from os.path import basename
from textwrap import wrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def get_val(a):
    # is it iterable?
    try:
        for val in a:
            if isinstance(val, (int, long, float)):
                return val

    except TypeError:
        if isinstance(a, (int, long, float)):
            return a

    return None

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

def get_stats(args, h_vcf, shared_info_keys, shared_fmt_keys):
    info_vcf = {}
    fmt_vcf  = {}

    vtype = None
    if args.vtype is not None:
        assert args.vtype in ('SNV', 'INDEL', 'SV')
        vtype = args.vtype

    for k_info in shared_info_keys:
        info_vcf[k_info] = []

    for rec in h_vcf:
        selected = True

        if vtype == 'SNV' and (not rec.is_snp or (rec.is_snp and rec.INFO.get('VT') == 'LOH')):
            output = False

        if vtype == 'INDEL' and not rec.is_indel:
            output = False

        if vtype == 'SV' and not rec.is_sv:
            output = False

        if args.passonly and rec.FILTER:
            selected = False

        if args.failonly and not rec.FILTER:
            selected = False

        if args.somaticonly and not is_somatic(rec):
            selected = False

        if args.germlineonly and is_somatic(rec):
            selected = False

        for k_info in shared_info_keys:
            if k_info in rec.INFO:
                value = get_val(rec.INFO.get(k_info))
                if selected and value is not None:
                    info_vcf[k_info].append(value)
        
        for sample in rec.samples:
            name = sample.sample
            fmt  = sample.data._asdict() # fmt is a collections.namedtuple
            if name not in fmt_vcf:
                fmt_vcf[name] = {}

            for k_fmt in shared_fmt_keys:
                if k_fmt in fmt:
                    if k_fmt not in fmt_vcf[name]:
                        fmt_vcf[name][k_fmt] = []
                    value = get_val(fmt[k_fmt])
                    if selected and value is not None:
                        fmt_vcf[name][k_fmt].append(value)

    return info_vcf, fmt_vcf

def cutoffs(tp, fp):
    tp.sort()
    fp.sort()

    min_range = None
    max_range = None

    if min(tp) < min(fp):
        for i in tp:
            if i >= min(fp):
                break
        max_range = i

    if min(fp) < min(tp):
        for i in fp:
            if i >= min(tp):
                break
        min_range = i

    if max(tp) > max(fp):
        for i in tp[::-1]:
            if i <= max(fp):
                break
        min_range = i

    if max(fp) > max(tp):
        for i in fp[::-1]:
            if i <= max(tp):
                break
        max_range = i

    return min_range, max_range

def main(args):
    h_vcf1 = vcf.Reader(filename=args.vcf[0])
    h_vcf2 = vcf.Reader(filename=args.vcf[1])

    if args.failonly == args.passonly == True:
        sys.exit("Error: specifying both -f/--failonly and -p/--passonly yields no results.")

    if args.somaticonly == args.germlineonly == True:
        sys.exit("Error: specifying both -s/--somaticonly and -g/--germlineonly yields no results.")

    shared_info_keys = set(h_vcf1.infos.keys()).intersection(h_vcf2.infos.keys())
    shared_fmt_keys  = set(h_vcf1.formats.keys()).intersection(h_vcf2.formats.keys())

    info_vcf1, fmt_vcf1 = get_stats(args, h_vcf1, shared_info_keys, shared_fmt_keys)
    info_vcf2, fmt_vcf2 = get_stats(args, h_vcf2, shared_info_keys, shared_fmt_keys)

    labels = map(str.upper, (args.label1, args.label2))
    vcf_TP = None
    vcf_FP = None

    if 'TP' in labels and 'FP' in labels:
        if args.label1 == 'TP':
            vcf_TP = 1
            vcf_FP = 2
        else:
            vcf_TP = 2 
            vcf_FP = 1

    filtout = None
    if args.filteroutfile is not None and 'TP' in labels and 'FP' in labels:
        filtout = open(args.filteroutfile, 'w')
        print "filterout:", args.filteroutfile

    for tag, values in info_vcf1.iteritems():
        # only compare things worth comparing
        if values and info_vcf2[tag] and len(set(values)) > 1 and len(set(info_vcf2[tag])) > 1:
            vcf1_values = np.asarray(values)
            vcf2_values = np.asarray(info_vcf2[tag])

            print '-'*60
            print 'INFO', tag, ':', h_vcf1.infos[tag].desc
            print basename(args.vcf[0]), 'INFO', tag, ','.join(map(str, vcf1_values))
            print basename(args.vcf[1]), 'INFO', tag, ','.join(map(str, vcf2_values))
            mwu = ss.mannwhitneyu(vcf1_values, vcf2_values)
            mwstring = "INFO (" + tag + ") Mann-Whitney U: " + "%0.1f" % mwu[0] + " P=" + "%0.3f" % mwu[1]
            print mwstring 

            if filtout is not None and mwu[1] < 0.05:
                assert vcf_TP != vcf_FP
                mincut = None
                maxcut = None
                if vcf_TP == 1:
                    (mincut, maxcut) = cutoffs(vcf1_values, vcf2_values)
                elif vcf_TP == 2:
                    (mincut, maxcut) = cutoffs(vcf2_values, vcf1_values)

                print "vcf_TP, vcf_FP:", vcf_TP, vcf_FP
                print "mincut, maxcut:", mincut, maxcut

                if mincut is not None:
                    filtout.write(' '.join(('INFO', tag, 'LTE', str(float(mincut)))) + "\n")
                if maxcut is not None:
                    filtout.write(' '.join(('INFO', tag, 'GT', str(float(maxcut)))) + "\n")

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title("INFO: " + tag + "\n(" + "\n".join(wrap(h_vcf1.infos[tag].desc)) + ")" + "\n" + mwstring)
            range = (min(vcf1_values.min(), vcf2_values.min()), max(vcf1_values.max(), vcf2_values.max()))
            ax.hist(vcf1_values, range=range, bins=20, alpha=0.3, label=args.label1, normed=True)
            ax.hist(vcf2_values, range=range, bins=20, alpha=0.3, label=args.label2, normed=True)
            ax.legend()
            plt.savefig(args.plotname + "_INFO_" + tag + "_" + ".png", bbox_inches='tight')

    for sample_name, fmt in fmt_vcf1.iteritems():
        for tag, values in fmt.iteritems():
            if values and fmt_vcf2[sample_name][tag] and len(set(values)) > 1 and len(set(fmt_vcf2[sample_name][tag])) > 1:
                vcf1_values = np.asarray(values)
                vcf2_values = np.asarray(fmt_vcf2[sample_name][tag])

                print '-'*60
                print 'FORMAT', tag, sample_name, ':', h_vcf1.formats[tag].desc
                print basename(args.vcf[0]), 'FORMAT', sample_name, tag, ','.join(map(str, vcf1_values))
                print basename(args.vcf[1]), 'FORMAT', sample_name, tag, ','.join(map(str, vcf2_values))
                mwu = ss.mannwhitneyu(vcf1_values, vcf2_values)
                mwstring = "FORMAT (" + tag + ") Mann-Whitney U: " + "%0.1f" % mwu[0] + " P=" + "%0.3f" % mwu[1]
                print mwstring 


                if filtout is not None and mwu[1] < 0.05:
                    assert vcf_TP != vcf_FP
                    mincut = None
                    maxcut = None
                    if vcf_TP == 1:
                        (mincut, maxcut) = cutoffs(vcf1_values, vcf2_values)
                    elif vcf_TP == 2:
                        (mincut, maxcut) = cutoffs(vcf2_values, vcf1_values)

                    print "vcf_TP, vcf_FP:", vcf_TP, vcf_FP
                    print "mincut, maxcut:", mincut, maxcut

                    if mincut is not None:
                        filtout.write(' '.join(('FORMAT', tag, 'LTE', str(float(mincut)), sample_name)) + "\n")
                    if maxcut is not None:
                        filtout.write(' '.join(('FORMAT', tag, 'GT', str(float(maxcut)), sample_name)) + "\n")

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_title("FORMAT: " + tag + " " + sample_name + "\n(" + "\n".join(wrap(h_vcf1.formats[tag].desc)) + ")" + "\n" + mwstring) 
                range = (min(vcf1_values.min(), vcf2_values.min()), max(vcf1_values.max(), vcf2_values.max()))
                ax.hist(vcf1_values, range=range, bins=20, alpha=0.3, label=args.label1, normed=True)
                ax.hist(vcf2_values, range=range, bins=20, alpha=0.3, label=args.label2, normed=True)
                ax.legend()
                plt.savefig(args.plotname + "_FORMAT_" + tag + "_" + sample_name + ".png", bbox_inches='tight')

    if filtout is not None:
        filtout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare info and format values from two VCF files')
    parser.add_argument(metavar='<vcf_file>', dest='vcf', nargs=2, help='files in VCF format')
    parser.add_argument('--label1', dest='label1', default='vcf1', help='label for vcf1') 
    parser.add_argument('--label2', dest='label2', default='vcf2', help='label for vcf2')
    parser.add_argument('--name', dest='plotname', default='plot', help='basename for plots')
    parser.add_argument('--filterout', dest='filteroutfile', default=None, help='output filters for filtervcf.py, labels must be TP and FP') 
    parser.add_argument('-t', '--vtype', dest='vtype', default=None, help='only include variants of vtype where vtype is SNV, INDEL, or SV')
    parser.add_argument('-p', '--passonly', action='store_true', default=False, help='only return PASS records')
    parser.add_argument('-f', '--failonly', action='store_true', default=False, help='only return non-PASS records')
    parser.add_argument('-s', '--somaticonly', action='store_true', default=False, help='only return somatic records')
    parser.add_argument('-g', '--germlineonly', action='store_true', default=False, help='only return germline records')
    args = parser.parse_args()
    main(args)
