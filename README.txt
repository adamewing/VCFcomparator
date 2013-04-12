VCF comparator: compares mutation calls in VCF formatted files

requires: 
PyVCF (https://github.com/jamescasbon/PyVCF)
numpy (http://numpy.scipy.org/)

usage: vcfcomparator.py [-h] [-m MASKFILE] [-o OUTDIR] [-t TRUTH] [-c CHROM]
                        [-s START] [-e END] [-v]
                        <vcf_file> <vcf_file>

Compares two sorted VCF files and (optionally) masks regions.

positional arguments:
  <vcf_file>            tabix-indexed files in VCF format

optional arguments:
  -h, --help            show this help message and exit
  -m MASKFILE, --mask MASKFILE
                        tabix-indexed BED file of masked intervals
  -o OUTDIR, --outdir OUTDIR
                        directory for output
  -t TRUTH, --truth TRUTH
                        also compare results to a "truth" VCF (should be
                        sorted and tabix-indexed)
  -c CHROM, --chrom CHROM
                        limit to one chromosome
  -s START, --start START
                        start position
  -e END, --end END     end position
  -v, --verbose         verbose mode for debugging

