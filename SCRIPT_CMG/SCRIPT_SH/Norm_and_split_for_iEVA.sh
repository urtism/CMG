#!bin/bash
### TOOLS ###
SCRIPT_PIPELINE=~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE
FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.15
PICARD=~/NGS_TOOLS/picard-tools-2.7.1/picard.jar

GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
FREEBAYES=~/NGS_TOOLS/freebayes/bin/freebayes
VARDICT=~/NGS_TOOLS/VarDictJava-master/build/install/VarDict/bin/VarDict
BCFTOOLS=bcftools
SURECALLTRIMMER=~/NGS_TOOLS/AGeNT/SurecallTrimmer_v3.5.1.46.jar
LOCATLT=~/NGS_TOOLS/AGeNT/LocatIt_v3.5.1.46.jar

VEP=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/filter_vep.pl

### DATABASES & FILES ###
LISTAFEATURES_GERMLINE=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_germline.list
LISTAFEATURES_SOMATIC=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_somatic.list
#LISTAFEATURES_SOMATIC=/home/jarvis/NGS_ANALYSIS/TARGET/Features_lists/lista_features_somatic_CF_20170315.list
ANN_LIST_GERMLINE=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_annotazione.list
ANN_LIST_SOMATIC=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_annotazione.list
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
REF_DICT=~/NGS_TOOLS/hg19/ucsc.hg19.dict
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/home/jarvis/NGS_TOOLS/hg19/common_all_20170710.vcf
LOGHI=~/git/CMG/LOGHI
TRASCR_CARDIO=~/NGS_ANALYSIS/TARGET/Lista_trascritti_Cardio.txt
TRASCR_BRCA=~/NGS_ANALYSIS/TARGET/Lista_trascritti_BRCA.txt
TRASCR_CANCER=~/NGS_ANALYSIS/TARGET/Lista_trascritti_Cancer.txt

### TARGET ###
TARGET_CARDIO_1000=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.list
TARGET_CARDIO_1000_BED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.bed
TARGET_CARDIO=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a.list
TARGET_CARDIO_BED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a.bed
TARGET_BRCA=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.list
TARGET_BRCA_BED=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.bed
TARGET_BRCA_FREEBAYES=~/NGS_ANALYSIS/TARGET/BRCA_FreeBayes_amplicon.bed
TARGET_EXOME_1000=~/NGS_ANALYSIS/TARGET/TruSight_One_v1.1_ESTESO+-1000.list
TARGET_EXOME_1000_BED=~/NGS_ANALYSIS/TARGET/TruSight_One_v1.1_ESTESO+-1000.bed
TARGET_CF=~/NGS_ANALYSIS/TARGET/ctDNA_2_113416_AmpliconsExport.list
TARGET_CF_BED=~/NGS_ANALYSIS/TARGET/ctDNA_2_113416_AmpliconsExport.bed
TARGET_CANCER_1000=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.list
TARGET_CANCER_1000_BED=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.bed
TARGET_HALOPLEX=~/NGS_ANALYSIS/TARGET/HaloPlex_Covered.list
TARGET_HALOPLEX_BED=~/NGS_ANALYSIS/TARGET/HaloPlex_Covered.bed
TARGET_SURESELECT=~/NGS_ANALYSIS/TARGET/Meloni_SureSelect_target_files/3066331/3066331_Covered.list
TARGET_SURESELECT_BED=~/NGS_ANALYSIS/TARGET/Meloni_SureSelect_target_files/3066331/3066331_Covered.bed
AMPLICONS_HALOPLEX_BED=~/NGS_ANALYSIS/TARGET/HaloPlex_target_files/HaloPlex_Amplicons.bed

cd /home/jarvis/Scrivania/VCF-IEVA/

for file in *.vcf
	do
		extract=$(echo $file| cut -d'_' -f 1,2,3)
		#printf "\n$extract\n" 

			if [[ $extract == *"FreeBayes"* ]];
				then
					printf "\n$extract\n"
					python $SCRIPT_PIPELINE/header_fix.py -v F -f $file \
					> /home/jarvis/Scrivania/VCF-IEVA/FIX/$extract.fix.vcf

					$BCFTOOLS norm -m -both \
					-f $REF \
					/home/jarvis/Scrivania/VCF-IEVA/FIX/$extract.fix.vcf \
					> /home/jarvis/Scrivania/VCF-IEVA/NORM/$extract.norm.vcf
				
			elif [[ $extract == *"GATK"* ]];
				then
					printf "\n$extract\n"
					python $SCRIPT_PIPELINE/header_fix.py -v G -f $file \
					> /home/jarvis/Scrivania/VCF-IEVA/FIX/$extract.fix.vcf

					$BCFTOOLS norm -m -both \
					-f $REF \
					/home/jarvis/Scrivania/VCF-IEVA/FIX/$extract.fix.vcf \
					> /home/jarvis/Scrivania/VCF-IEVA/NORM/$extract.norm.vcf

			elif [[ $extract == *"VarScan"* ]];
				then
					printf "\n$extract\n"
					python $SCRIPT_PIPELINE/header_fix.py -v V -f $file \
					> /home/jarvis/Scrivania/VCF-IEVA/FIX/$extract.fix.vcf

					$BCFTOOLS norm -m -both \
					-f $REF \
					/home/jarvis/Scrivania/VCF-IEVA/FIX/$extract.fix.vcf \
					> /home/jarvis/Scrivania/VCF-IEVA/NORM/$extract.norm.vcf

			else
				printf "\n$file non contiene alcun file GATK, FREEBAYES o VARSCAN\n"

			fi

	done