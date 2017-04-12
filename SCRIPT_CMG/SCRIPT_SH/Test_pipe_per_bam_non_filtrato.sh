#!bin/bash

FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.15
PICARD=~/NGS_TOOLS/picard-tools-2.9.0/picard.jar
GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_138.hg19.vcf
INPUT=~/NGS_ANALYSIS/INPUT_DATA/CARDIO
INPUTBRCA=~/NGS_ANALYSIS/INPUT_DATA/BRCA
INPUTCANCER=~/NGS_ANALYSIS/INPUT_DATA/CANCER
INPUTEXOME=~/NGS_ANALYSIS/INPUT_DATA/EXOME
PROCESSING=~/NGS_ANALYSIS/PROCESSING
TARGET=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.list
TARGETCARDIOBED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.bed
TARGETMETRICS=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a.list
TARGETBRCA=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.list
TARGETBRCABED=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.bed
TARGBRCAFREE=~/NGS_ANALYSIS/TARGET/BRCA_FreeBayes_amplicon.bed
TARGETEXOME=~/NGS_ANALYSIS/TARGET/TruSight_One_v1.1_ESTESO+-1000.list
TARGETEXOMEBED=~/NGS_ANALYSIS/TARGET/TruSight_One_v1.1_ESTESO+-1000.bed
TARGETCANCER=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.list
TARGETCANCERBED=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.bed
OUTVCF=~/NGS_ANALYSIS/OUTPUT_DATA
STORAGE=~/NGS_ANALYSIS/STORAGE
VEP=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/filter_vep.pl

DataRun=20161125

# cd $INPUT
# for file in *\_L001_R1_001.fastq.gz
# do

# extract=$(echo $file| cut -d'_' -f 1,2,3)

# $BWA/bwa mem $REF \
# -M $INPUT/$extract\_L001_R1_001.fastq.gz \
# -t 4 \
# $INPUT/$extract\_L001_R2_001.fastq.gz > $PROCESSING/1_Align/Aligned_Sam/$extract.sam

# java -Xmx8g -jar $PICARD SamFormatConverter \
# I=$PROCESSING/1_Align/Aligned_Sam/$extract.sam \
# O=/home/minime/NGS_ANALYSIS/PROCESSING/1_Align/Converted_bam/$extract\_Converted.bam

# printf $"\n~~~>	Sample $extract => Sam Format Converter: DONE\n\n"

# printf $"\n~~~>	Sample $extract => Sort Sam\n\n"

# java -Xmx8g -jar $PICARD SortSam \
# I=/home/minime/NGS_ANALYSIS/PROCESSING/1_Align/Converted_bam/$extract\_Converted.bam \
# O=/home/minime/NGS_ANALYSIS/PROCESSING/1_Align/Sort_bam/$extract\_Sort.bam \
# SORT_ORDER=coordinate
# done


# cd $PROCESSING/1_Align/Sort_bam/

# for filename in *.bam
# do
# java -Xmx8g -jar $PICARD AddOrReplaceReadGroups \
# I=${filename%.*}.bam \
# O=$PROCESSING/2_Add/${filename%.*}_Add.bam \
# RGID=${filename%.*} \
# RGPL=ILLUMINA \
# RGSM=${filename%.*} \
# RGPU=ILLUMINA_$DATE \
# RGLB=Cardio \
# VALIDATION_STRINGENCY=LENIENT

# java -Xmx8g -jar $PICARD MarkDuplicates \
# I=$PROCESSING/2_Add/${filename%.*}_Add.bam \
# O=$PROCESSING/3_Mark/${filename%.*}_debup.bam \
# METRICS_FILE=$PROCESSING/3_Mark/${filename%.*}_Metrics.txt \
# READ_NAME_REGEX=null \
# VALIDATION_STRINGENCY=LENIENT \
# ASSUME_SORTED=true

# java -Xmx8g -jar $PICARD BuildBamIndex \
# I=$PROCESSING/3_Mark/${filename%.*}_debup.bam \
# O=$PROCESSING/3_Mark/${filename%.*}_debup.bai \
# VALIDATION_STRINGENCY=LENIENT

# printf $"\n~~~>	Sample ${filename%.*} => Realigner Target Creator: DONE\n\n"

# printf $"\n~~~>	Sample ${filename%.*} => Indel Realigner\n\n"

# java -Xmx8g -jar $GATK -T RealignerTargetCreator \
# -R $REF \
# -I $PROCESSING/3_Mark/${filename%.*}_debup.bam \
# -o $PROCESSING/4_Indel/${filename%.*}_IndelRealigner.intervals \
# --known $MILLS \
# -L $TARGET

# printf $"\n~~~>	Sample ${filename%.*} => Realigner Target Creator: DONE\n\n"

# printf $"\n~~~>	Sample ${filename%.*} => Indel Realigner\n\n"

# java -Xmx8g -jar $GATK -T IndelRealigner \
# -R $REF \
# -I $PROCESSING/3_Mark/${filename%.*}_debup.bam \
# -targetIntervals $PROCESSING/4_Indel/${filename%.*}_IndelRealigner.intervals \
# -o $PROCESSING/4_Indel/${filename%.*}_Realigned.bam \
# -known $MILLS

# echo $'\n =========>	Indel Realigner: DONE\n'
		
# printf $"\n~~~>	Sample ${filename%.*} => Base Quality Score Recalibration (BQSR): Base Recalibrator\n\n"

# cd $PROCESSING/4_Indel/
# for filename in *.bam
# do

# java -Xmx8g -jar $GATK -T BaseRecalibrator \
# -R $REF \
# -I $PROCESSING/4_Indel/${filename%.*}.bam \
# -knownSites $MILLS \
# -o $PROCESSING/5_BQSR/${filename%.*}_Recal.table \
# -L $TARGET

# printf $"\n~~~>	Sample ${filename%.*} => Base Quality Score Recalibration (BQSR): Base Recalibrator: DONE\n\n"
# cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
# printf $"\n~~~>	Sample ${filename%.*} => Base Quality Score Recalibration (BQSR): Print reads\n\n"

# java -Xmx8g -jar $GATK -T PrintReads \
# -R $REF \
# -I $PROCESSING/4_Indel/${filename%.*}.bam \
# -BQSR $PROCESSING/5_BQSR/${filename%.*}_Recal.table  \
# -o $PROCESSING/5_BQSR/${filename%.*}.bam \
# -L $TARGET

# done

# #---------------GATK

cd $PROCESSING/5_BQSR/
for filename in *.bam
do

	cat ~/Scrivania/SCRIPT_PIPELINE/logo_variant.txt
	printf $"\n =========>	Sample ${filename%.*} => Variant Calling: Haplotype Caller\n\n"

	java -Xmx8g -jar $GATK -T HaplotypeCaller \
	-R $REF \
	-I $PROCESSING/5_BQSR/${filename%.*}.bam \
	-o $PROCESSING/6_Variant/GATK/${filename%.*}.g.vcf \
	-ERC GVCF \
	-mbq 1 \
	-minReadsPerAlignStart 1 \
	-stand_call_conf 1 \
	--doNotRunPhysicalPhasing \
	-mmq 1 \
	-drf DuplicateRead \
	-bamout $PROCESSING/6_Variant/GATK/${filename%.*}.g.vcf.bam \
	-L $TARGET

	printf $"\n =========>	Sample ${filename%.*} => Variant Calling: Haplotype Caller: DONE\n\n"


done

#-------------FREEB

# 	cd $PROCESSING/5_BQSR/
# 	ls *.bam > Bam_list.txt
# 	ls *.bam > Sample_list.txt
# 	sed -i -e "s/.bam//g" Sample_list.txt

# 	~/NGS_TOOLS/freebayes/bin/freebayes -f $REF \
# 	-L Bam_list.txt \
# 	-K \
# 	-J \
# 	-s Sample_list.txt \
# 	--genotype-qualities \
# 	--report-genotype-likelihood-max \
# 	--allele-balance-priors-off \
# 	-t $TARGETCARDIOBED > $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes.vcf

# 	echo $'\n =========>	Variant Calling with FreeBayes: DONE\n'

# 	echo $'\n =========>	Variant Calling with VarScan2\n\n'

# #----------varsc


# 	cd $PROCESSING/5_BQSR
# 	samtools mpileup -B -q 1 -f $REF \
# 	-l $TARGETCARDIOBED -b $PROCESSING/5_BQSR/Bam_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup

# 	cp Bam_list.txt Sample_list.txt
# 	sed -i -e 's/.bam//g' Sample_list.txt

# 	java -jar -Xmx8g $VARSCAN mpileup2snp $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup \
# 	--min-coverage 10 \
# 	--min-var-freq 0.20 \
# 	--pvalue 0.05 \
# 	--output-vcf 1 \
# 	--vcf-sample-list Sample_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf

# 	java -jar -Xmx8g $VARSCAN mpileup2indel $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup \
# 	--min-coverage 10 \
# 	--min-var-freq 0.10 \
# 	--pvalue 0.1 \
# 	--output-vcf 1 \
# 	--vcf-sample-list Sample_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf
