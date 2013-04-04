#!/bin/sh

if [ $# -ne 2 ]
then
    echo "usage: $0 <directory containing VCF files> <reference fasta>"
    exit 65;
fi

if [ ! -e $1 ]
then
    echo "$1 does not exist"
    exit 65;
fi

if [ ! -e $2 ]
then
    echo "$2 does not exist"
    exit 65;
    if [ ! -e $2.fai ]
    then
        echo "$2.fai not found: reference should in indexed with samtools faidx"
    fi
fi

tmpfile=tmp.$RANDOM.vcf
echo "using tmp file: $tmpfile"

echo "un-gzipping..."
ls -1 $1/*.vcf.gz | xargs -n 1 -P 8 gzip -d

for vcf in `ls -1 $1/*.vcf`
do
    echo "refsorting $vcf..."
    ./sort_vcf_by_ref.py $vcf $2.fai > $tmpfile 
    mv $tmpfile $vcf
done

echo "left shift indels..."
for vcf in `ls -1 $1/*.vcf`
do
    java -Xmx2g -jar /cluster/home/ewingad/GATK/GenomeAnalysisTK.jar -R $2 -T LeftAlignVariants --variant $vcf -o $tmpfile 
    mv $tmpfile $vcf
done

echo "left shift breakends..."
for vcf in `ls -1 $1/*.vcf`
do
    ./leftShiftBreakends.py -r $2 -o $tmpfile -v $vcf
    mv $tmpfile $vcf
done

echo "sorting VCFs..."
for vcf in `ls -1 $1/*.vcf`
do
    echo $vcf
    vcf-sort $vcf > $tmpfile
    mv $tmpfile $vcf
done

echo "bgzipping..."
ls -1 $1/*.vcf | xargs -n 1 -P 8 bgzip

echo "indexing..."
ls -1 $1/*.vcf.gz | xargs -n 1 -P 8 tabix -f -p vcf
