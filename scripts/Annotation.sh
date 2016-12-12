#!/bin/bash

#commento

cat ~/Scrivania/SCRIPT_PIPELINE/logo.txt 

FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.12
PICARD=~/NGS_TOOLS/picard-tools-2.3.0/picard.jar
GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_135.hg19.vcf
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

	cd $VEP

	perl $VEPANN -i $PROCESSING/7_Filter/20161103_BRCA_GATK_Filter.vcf \
	-o $PROCESSING/8_Annotation/20161103_BRCA_GATK_Filter_ANN.vcf \
	--stats_file $PROCESSING/8_Annotation/20161103_BRCA_GATK_Filter_ANN.html \
	--cache \
	--everything \
	--dont_skip \
	--assembly GRCh37 \
	--offline \
	--force_overwrite \
	-v \
	--fork 10 \
	--variant_class \
	--allele_number \
	--total_length \
	--vcf_info_field ANN \
	--numbers \
	--hgvs \
	--protein \
	--canonical \
	--check_existing \
	--check_alleles \
	--check_svs \
	--gmaf \
	--tsl \
	--pubmed \
	--species homo_sapiens \
	--failed 1 \
	--vcf

