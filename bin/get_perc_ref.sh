#!/usr/bin/env bash

id=$1
refseg=$2
tupleid=$3

bcftools mpileup -Ou --per-sample-mF -d 10000 -B -f ${refseg}.fasta ${id}.${refseg}.bam | bgzip > ${id}.${refseg}.region.pileup.gz
tabix ${id}.${refseg}.region.pileup.gz
bcftools call ${id}.${refseg}.region.pileup.gz -T ${tupleid}.${refseg}.posmatrix -P 1.0e-2 -m -Ov | bcftools norm -m +any >  ${id}.${refseg}.region.vcf
bcftools filter ${id}.${refseg}.region.vcf -e "DP=0" | bcftools stats ${id}.${refseg}.region.vcf> ${id}.${tupleid}.${refseg}.stat
num=`grep "number of no-ALTs:" ${id}.${tupleid}.${refseg}.stat | awk '{print $6}'`
den=`wc -l ${tupleid}.${refseg}.posmatrix | awk '{print $1}'`
if [ ${den} -eq 0 ]; then
  perc=0
else
  perc=`echo "scale=4;$num/$den" | bc`
fi
echo "${perc}"
