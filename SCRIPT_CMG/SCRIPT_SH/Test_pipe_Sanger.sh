#!bin/bash

cat ~/Scrivania/SCRIPT_PIPELINE/logo.txt 

FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.15
PICARD=~/NGS_TOOLS/picard-tools-2.7.1/picard.jar
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


$BWA/bwa mem $REF \
/home/jarvis/Scrivania/Test_sanger.fastq \
-t 10 > /home/jarvis/Scrivania/Test_Sanger/Test_sanger.sam

java -Xmx64g -jar $PICARD SamFormatConverter \
I=/home/jarvis/Scrivania/Test_Sanger/Test_sanger.sam \
O=/home/jarvis/Scrivania/Test_Sanger/Test_sanger_Converted.bam

printf $"\n~~~>	Sample $extract => Sam Format Converter: DONE\n\n"
cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
printf $"\n~~~>	Sample $extract => Sort Sam\n\n"

java -Xmx64g -jar $PICARD SortSam \
I=/home/jarvis/Scrivania/Test_Sanger/Test_sanger_Converted.bam \
O=/home/jarvis/Scrivania/Test_Sanger/Test_sanger_Sort.bam \
SORT_ORDER=coordinate

java -Xmx64g -jar $PICARD AddOrReplaceReadGroups \
I=/home/jarvis/Scrivania/Test_Sanger/Test_sanger_Sort.bam \
O=/home/jarvis/Scrivania/Test_Sanger/Test_sanger_Add.bam \
RGID=Test_Sanger_01 \
RGPL=Sanger \
RGSM=Test_01 \
RGPU=Sanger_0209 \
RGLB=Ex1 \
VALIDATION_STRINGENCY=LENIENT

java -Xmx64g -jar $PICARD BuildBamIndex \
I=/home/jarvis/Scrivania/Test_Sanger/Test_sanger_Add.bam \
O=/home/jarvis/Scrivania/Test_Sanger/Test_sanger_Add.bai \
VALIDATION_STRINGENCY=LENIENT

	java -Xmx64g -jar $GATK -T RealignerTargetCreator \
		-R $REF \
		-I /home/jarvis/Scrivania/Test_Sanger/Test_sanger_Add.bam \
		-o /home/jarvis/Scrivania/Test_Sanger/Test_sanger_IndelRealigner.intervals \
		--known $MILLS

		printf $"\n~~~>	Sample ${filename%.*} => Realigner Target Creator: DONE\n\n"
		cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
		printf $"\n~~~>	Sample ${filename%.*} => Indel Realigner\n\n"

		java -Xmx64g -jar $GATK -T IndelRealigner \
		-R $REF \
		-I /home/jarvis/Scrivania/Test_Sanger/Test_sanger_Add.bam \
		-targetIntervals /home/jarvis/Scrivania/Test_Sanger/Test_sanger_IndelRealigner.intervals \
		-o /home/jarvis/Scrivania/Test_Sanger/Test_sanger_Realigned.bam \
		-known $MILLS

		echo $'\n =========>	Indel Realigner: DONE\n'
		cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
		printf $"\n~~~>	Sample ${filename%.*} => Base Quality Score Recalibration (BQSR): Base Recalibrator\n\n"


samtools view /home/jarvis/Scrivania/Test_Sanger/Test_sanger_Realigned.bam > /home/jarvis/Scrivania/Test_Sanger/bam_view.txt


java -Xmx64g -jar $GATK -T HaplotypeCaller \
-R $REF \
-I /home/jarvis/Scrivania/Test_Sanger/Test_sanger_Realigned.bam \
-o /home/jarvis/Scrivania/Test_Sanger/Test_sanger.g.vcf \
-ERC GVCF \
--doNotRunPhysicalPhasing \
--min_base_quality_score 1

/home/jarvis/NGS_TOOLS/freebayes/bin/freebayes -f $REF \
/home/jarvis/Scrivania/Test_Sanger/Test_sanger_Realigned.bam \
--genotype-qualities \
--report-genotype-likelihood-max \
--allele-balance-priors-off \
-m 0 \
-q 0 \
-R 0 \
-Y 0 \
-Q 1 \
-F 0.001 \
-C 1 > $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes.vcf


samtools mpileup -B -q 1 -f $REF /home/jarvis/Scrivania/Test_Sanger/Test_sanger_Realigned.bam  > /home/jarvis/Scrivania/Test_Sanger/Test_sanger.mpileup

	java -jar -Xmx64g $VARSCAN mpileup2snp /home/jarvis/Scrivania/Test_Sanger/Test_sanger.mpileup \
	--min-coverage 1 \
	--min-var-freq 0.001 \
	--min-reads2 1 \
	--min-avg-qual 0 \
	--strand-filter 0 \
	--output-vcf 1 > /home/jarvis/Scrivania/Test_Sanger/Test_sanger_Varscan_snp.vcf

	java -jar -Xmx64g $VARSCAN mpileup2indel /home/jarvis/Scrivania/Test_Sanger/Test_sanger.mpileup \
	--min-coverage 1 \
	--min-var-freq 0.001 \
	--min-reads2 1 \
	--min-avg-qual 0 \
	--strand-filter 0 \
	--output-vcf 1 > /home/jarvis/Scrivania/Test_Sanger/Test_sanger_Varscan_Indel.vcf



