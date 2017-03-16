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
INPUT=~/NGS_ANALYSIS_2/INPUT_DATA/CARDIO
INPUTBRCA=~/NGS_ANALYSIS_2/INPUT_DATA/BRCA
INPUTCANCER=~/NGS_ANALYSIS_2/INPUT_DATA/CANCER
INPUTEXOME=~/NGS_ANALYSIS_2/INPUT_DATA/EXOME
PROCESSING=~/NGS_ANALYSIS_2/PROCESSING
TARGET=~/NGS_ANALYSIS_2/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.list
TARGETCARDIOBED=~/NGS_ANALYSIS_2/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.bed
TARGETMETRICS=~/NGS_ANALYSIS_2/TARGET/trusight_cardio_manifest_a.list
TARGETBRCA=~/NGS_ANALYSIS_2/TARGET/AFP2_manifest_v1.list
TARGETBRCABED=~/NGS_ANALYSIS_2/TARGET/AFP2_manifest_v1.bed
TARGBRCAFREE=~/NGS_ANALYSIS_2/TARGET/BRCA_FreeBayes_amplicon.bed
TARGETEXOME=~/NGS_ANALYSIS_2/TARGET/TruSight_One_v1.1_ESTESO+-1000.list
TARGETEXOMEBED=~/NGS_ANALYSIS_2/TARGET/TruSight_One_v1.1_ESTESO+-1000.bed
TARGETCANCER=~/NGS_ANALYSIS_2/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.list
TARGETCANCERBED=~/NGS_ANALYSIS_2/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.bed
OUTVCF=~/NGS_ANALYSIS_2/OUTPUT_DATA
STORAGE=~/NGS_ANALYSIS_2/STORAGE
VEP=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/filter_vep.pl



	# perl $VEPANN -i /home/jarvis/NGS_ANALYSIS_2/PROCESSING/6_Variant/GATK/NORM_GATK/20170217_01_T_Cardio_Sane_GATK.split.vcf \
	# -o $PROCESSING/8_Annotation/20170217_01_T_Cardio_Sane_GATK.vcf \
	# --stats_file $PROCESSING/8_Annotation/20170217_01_T_Cardio_Sane_GATK.html \
	# --vcf \
	# --cache \
	# --assembly GRCh37 \
	# --offline \
	# --force_overwrite \
	# -v \
	# --fork 2 \
	# --variant_class \
	# --sift b \
	# --poly b \
	# --vcf_info_field ANN \
	# --hgvs \
	# --protein \
	# --canonical \
	# --check_existing \
	# --gmaf \
	# --pubmed \
	# --species homo_sapiens \
	# --failed 1 \
	# --vcf \
	# --plugin ExAC,/home/jarvis/.vep/Plugins/ExAC/ExAC.r0.3.1.sites.vep.vcf.gz \
	# --pick


cd /home/jarvis/NGS_ANALYSIS_2/PROCESSING/8_Annotation/CONFERME


	for file in *.tsv
		
		do
			python /home/jarvis/git/CMG/SCRIPT_CMG/SCRIPT_PYTHON/Database_Annotation.py \
			-i $file \
			-COSMIC path \
			-clnv /home/jarvis/DB_ANNOTAZIONE/ClinVar/clinvar.target_cardio.vcf,CLNSIG,CLNDB,CLNORIGIN \
			-o /home/jarvis/NGS_ANALYSIS_2/OUTPUT_DATA/${filename%.*}_Ann.tsv \

		done