#!bin/bash

DIR=~/Scrivania/TEST/VCF
DIR_OUT=~/Scrivania/TEST/VCF/OUT
REFERENCE=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
SCRIPTS_DIR=~/git/CMG/scripts
TOOLS_DIR=~/NGS_TOOLS

mkdir -p $DIR_OUT

#cd $DIR
for i in $DIR/*.vcf; do

filename=$(basename "$i")
echo $filename
extension="${filename##*.}"
filename_no_ext="${filename%.*}"

echo $filename_no_ext

if [[ $filename =~ .*FreeBayes.* ]]
	then
		python $SCRIPTS_DIR/header_fix.py \
		-f $i \
		-v F \
		> $DIR_OUT/$filename_no_ext.fixed.vcf
elif [[ $filename =~ .*VarScan.* ]]
	then 
		cp $i $DIR_OUT/$filename_no_ext.fixed.vcf
		
elif [[ $filename =~ .*GATK.* ]]
	then
		python $SCRIPTS_DIR/header_fix.py \
		-f $i \
		-v G \
		> $DIR_OUT/$filename_no_ext.fixed.vcf
fi

$TOOLS_DIR/bcftools-1.3.1/bcftools norm -m -both \
-f $REFERENCE \
$DIR_OUT/$filename_no_ext.fixed.vcf \
> $DIR_OUT/$filename_no_ext.split.vcf


    
done



