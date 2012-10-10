VCF comparator: compares mutation calls in VCF formatted files

usage: vcfcomparator.py [-h] [-m MASKFILE] <vcf_file> <vcf_file>

Compares two sorted VCF files and (optionally) masks regions.

positional arguments:
  <vcf_file>            tabix-indexed files in VCF format

optional arguments:
  -h, --help            show this help message and exit
  -m MASKFILE, --mask MASKFILE
                        BED file of masked interval
