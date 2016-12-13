#!/bin/bash


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

	perl $VEPANN -i /home/jarvis/NGS_ANALYSIS/PROCESSING/7_Filter/20160520_Cardio_GATK_Filter.vcf \
	-o $PROCESSING/8_Annotation/20160520_Cardio_GATK_Filter_ANN_TABBBBB.tsv \
	--stats_file $PROCESSING/8_Annotation/20160520_Cardio_GATK_Filter_ANN.html \
	--cache \
	--dont_skip \
	--assembly GRCh37 \
	--force_overwrite \
	-v \
	--fork 200 \
	--variant_class \
	--total_length \
	--offline \
	--vcf_info_field ANN \
	--numbers \
	--sift b \
	--polyphen b \
	--maf_esp \
	--symbol \
	--hgvs \
	--canonical \
	--check_existing \
	--check_alleles \
	--check_svs \
	--gmaf \
	--maf_1kg \
	--pubmed \
	--species homo_sapiens \
	--failed 1 \
	--plugin Blosum62 \
	--plugin dbNSFP,/home/jarvis/.vep/Plugins/dbNSFP_2.9/dbNSFP.gz,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,VEST3_score, \
	--plugin dbscSNV,/home/jarvis/.vep/Plugins/dbscSNV/dbscSNV.txt.gz \
	--plugin ExAC,/home/jarvis/.vep/Plugins/ExAC/ExAC.r0.3.1.sites.vep.vcf.gz \
	--plugin GeneSplicer,/home/jarvis/.vep/Plugins/GeneSplicer/GeneSplicer/bin/alpha/genesplicer,/home/jarvis/.vep/Plugins/GeneSplicer/GeneSplicer/training_data_sets/Human,context=200 \
	--plugin Gwava,tss,/home/jarvis/.vep/Plugins/Gwava/gwava_scores.bed.gz \
	--plugin HGVSshift \
	--plugin MaxEntScan,/home/jarvis/.vep/Plugins/MaxEntScan/fordownload \
	--plugin SameCodon \
	--tab

