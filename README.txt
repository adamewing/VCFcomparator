VCF comparator: compares mutation calls in VCF formatted files

requires: 
PyVCF (https://github.com/jamescasbon/PyVCF)
numpy (http://numpy.scipy.org/)

usage: vcfcomparator.py [-h] [-m MASKFILE] [-w] <vcf_file> <vcf_file>

Compares two sorted VCF files and (optionally) masks regions.

positional arguments:
  <vcf_file>            tabix-indexed files in VCF format

optional arguments:
  -h, --help            show this help message and exit
  -m MASKFILE, --mask MASKFILE
                        BED file of masked intervals
  -w, --weight_intervals
                        apply a normally-distributed weight to interval match
                        scores

