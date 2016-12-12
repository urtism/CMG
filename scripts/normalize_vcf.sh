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

<<<<<<< HEAD
if [[ $filename =~ .*FreeBayes.* ]bas]
	then
		python $SCRIPTS_DIR/header_fix.py \
		-f $i \
		-v F \
		> $DIR_OUT/$filename_no_ext.fixed.vcf
elif [[ $filename =~ .*VarScan.* ]]
	then 
		cp $i $DIR_OUT/$filename_no_ext.fixed.vcf
=======
	python $SCRIPTS_DIR/header_fix.py \
	-f $i \
	-v F \
	> $DIRFREE/NORM_FREE/$filename_no_ext.fixed.vcf
>>>>>>> branch 'devel' of https://github.com/urtism/CMG.git
		
	$BCFTOOLS norm -m -both \
	-f $REFERENCE \
	$DIRFREE/NORM_FREE/$filename_no_ext.fixed.vcf \
	> $DIRFREE/NORM_FREE/$filename_no_ext.split.vcf	

done


for i in $DIRVARSCAN/*.vcf;
do

	filename=$(basename "$i")
	extension="${filename##*.}"
	filename_no_ext="${filename%.*}"

	python $SCRIPTS_DIR/header_fix.py \
	-f $i \
	-v V \
	> $DIRVARSCAN/NORM_VARSCAN/$filename_no_ext.fixed.vcf
		
	$BCFTOOLS norm -m -both \
	-f $REFERENCE \
	$DIRVARSCAN/NORM_VARSCAN/$filename_no_ext.fixed.vcf \
	> $DIRVARSCAN/NORM_VARSCAN/$filename_no_ext.split.vcf	

done


for i in $DIRGATK/*.vcf;
do

	filename=$(basename "$i")
	extension="${filename##*.}"
	filename_no_ext="${filename%.*}"

	python $SCRIPTS_DIR/header_fix.py \
	-f $i \
	-v G \
	> $DIRGATK/NORM_GATK/$filename_no_ext.fixed.vcf
		
	$BCFTOOLS norm -m -both \
	-f $REFERENCE \
	$DIRGATK/NORM_GATK/$filename_no_ext.fixed.vcf \
	> $DIRGATK/NORM_GATK/$filename_no_ext.split.vcf	

done



