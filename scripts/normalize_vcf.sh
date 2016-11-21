#!bin/bash

DIRFREE=~/NGS_ANALYSIS/PROCESSING/6_Variant/FreeBayes
DIRVARSCAN=~/NGS_ANALYSIS/PROCESSING/6_Variant/VarScan
DIRGATK=~/NGS_ANALYSIS/PROCESSING/6_Variant/GATK
REFERENCE=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
SCRIPTS_DIR=~/git/CMG/scripts
TOOLS_DIR=~/NGS_TOOLS
BCFTOOLS=~/NGS_TOOLS/bcftools-1.3.1/bcftools


for i in $DIRFREE/*.vcf;
do

	filename=$(basename "$i")
	extension="${filename##*.}"
	filename_no_ext="${filename%.*}"

	python $SCRIPTS_DIR/header_fix.py \
	-f $i \
	-v F \
	> $DIRFREE/$filename_no_ext.fixed.vcf
		
	$BCFTOOLS norm -m -both \
	-f $REFERENCE \
	$DIRFREE/$filename_no_ext.fixed.vcf \
	> $DIRFREE/$filename_no_ext.split.vcf	

done


for i in $DIRVARSCAN/*.vcf;
do

	filename=$(basename "$i")
	extension="${filename##*.}"
	filename_no_ext="${filename%.*}"

	python $SCRIPTS_DIR/header_fix.py \
	-f $i \
	-v V \
	> $DIRVARSCAN/$filename_no_ext.fixed.vcf
		
	$BCFTOOLS norm -m -both \
	-f $REFERENCE \
	$DIRVARSCAN/$filename_no_ext.fixed.vcf \
	> $DIRVARSCAN/$filename_no_ext.split.vcf	

done


for i in $DIRGATK/*.vcf;
do

	filename=$(basename "$i")
	extension="${filename##*.}"
	filename_no_ext="${filename%.*}"

	python $SCRIPTS_DIR/header_fix.py \
	-f $i \
	-v G \
	> $DIRGATK/$filename_no_ext.fixed.vcf
		
	$BCFTOOLS norm -m -both \
	-f $REFERENCE \
	$DIRGATK/$filename_no_ext.fixed.vcf \
	> $DIRGATK/$filename_no_ext.split.vcf	

done



