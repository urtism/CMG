#!bin/bash

PipeVersion=3.0

### IMPORT ###
# source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/BRCA_function.sh
# source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/Cancer_function.sh
# source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/Cardio_function.sh
# source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/Cardio_function_all.sh
# source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/Exome_function.sh
# source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/Pipeline_CellFree_Germline.sh
# source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/Pipeline_CellFree_Somatic.sh
# source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/Pipeline_CellFree_Synthetic.sh
# source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/Pipeline_Germline_all.sh
# source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/Pipeline_Somatic.sh
source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/ALLINEAMENTO.sh
source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/VARIANT_CALLING.sh
source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/PREPROCESSING.sh

### TOOLS ###

SCRIPT_PIPELINE=~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE
FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.15
PICARD=~/NGS_TOOLS/picard-tools-2.7.1/picard.jar

GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
FREEBAYES=~/NGS_TOOLS/freebayes/bin/freebayes
#VARDICT=

VEP=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/filter_vep.pl

### DATABASES & FILES ###
LISTAFEATURES_GERMLINE=
LISTAFEATURES_SOMATIC=
ANN_LIST_GERMLINE=
ANN_LIST_SOMATIC=
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_138.hg19.vcf
LOGHI=~/git/CMG/LOGHI

### TARGET ###
TARGET_CARDIO_1000=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.list
TARGET_CARDIO_1000_BED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.bed
TARGET_CARDIO=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a.list
TARGET_BRCA=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.list
TARGET_BRCA_BED=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.bed
TARGET_BRCA_FREEBAYES=~/NGS_ANALYSIS/TARGET/BRCA_FreeBayes_amplicon.bed
TARGET_EXOME_1000=~/NGS_ANALYSIS/TARGET/TruSight_One_v1.1_ESTESO+-1000.list
TARGET_EXOME__1000_BED=~/NGS_ANALYSIS/TARGET/TruSight_One_v1.1_ESTESO+-1000.bed
TARGET_CANCER_1000=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.list
TARGET_CANCER_1000_BED=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.bed



if [ "$#" == "0" ]
then
	echo "usa -h"
	exit 1;
else
	while [[ $# -gt 0 ]]
	do
	key="$1"
	case $key in
		-a|--analisi)
		ANALISI="$2"
		shift # past argument
		;;
		-w|--workdir)
		WORKDIR="$2"
		shift # past argument
		;;
		--cfg)
		FILEPATH="$2"
		shift # past argument
		;;
		-d|--data)
		DATA="$2"
		shift # past argument
		;;
		-r|--run)
		RUN="$2"
		shift # past argument
		;;
		-p|--pannello)
		PANNELLO="$2"
		shift # past argument
		;;
		-h|--help)
		echo -e "questo è l'help"
		exit 1;
		;;
		--default)
		DEFAULT=YES
		;;
		*)
		echo "wrong"
		exit 1;
		;;
	esac
	shift 
	done
fi

DATE=$(date +"%Y%m%d");
STARTTIME=$(date +%s)

NGSDIR=~/NGS_ANALYSIS_TEMP

var=1
while [ -d $NGSDIR/$DATA\_Run_$RUN\_$PANNELLO\_$var ]; do
    ((var++))
done

if [ "$WORKDIR" == "" ]
then
	WORKDIR=$NGSDIR/$DATA\_Run_$RUN\_$PANNELLO\_$var
fi
mkdir -p $WORKDIR



PIPELINE_GERMLINE() {

	ALLINEAMENTO

	if [ "$PANNELLO" == "Cardio" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CARDIO_1000
		TRANSCR_LIST=
		echo "Cardio"

	elif [ "$PANNELLO" == "Cancer" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CANCER_1000
		TRANSCR_LIST=
		echo "Cancer"

	elif [ "$PANNELLO" == "Exome" ] 
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_EXOME_1000
		TRANSCR_LIST=
		echo "Exome"

	elif [ "$PANNELLO" == "BRCA" ]
	then
		DESIGN="AMPLICON"
		TARGET=$TARGET_BRCA
		TRANSCR_LIST=
		echo "BRCA"

	fi

	PREPROCESSING
	VARIANT_CALLING_GERMLINE

}



if [ "$ANALISI" == "Germline" ] || [ "$ANALISI" == "germline" ] || [ "$ANALISI" == "GERMLINE" ]
then
	PIPELINE_GERMLINE

elif [ "$ANALISI" == "Somatic" ] || [ "$ANALISI" == "somatic" ] || [ "$ANALISI" == "SOMATIC" ]
then
	echo "somatic"
	PIPELINE_SOMATIC
elif [ "$ANALISI" == "Cellfree" ] || [ "$ANALISI" == "cellfree" ] || [ "$ANALISI" == "CELLFREE" ]
then
	echo "cellfree"
	#PIPELINE_CELLFREE

fi

cat $FILEPATH | while read line
do

	som_name=$(echo "$line" | cut -f1)
	sane_name=$(echo "$line" | cut -f2)
	#echo $line
	#echo -e "$som_name"
	#echo -e "$sane_name"
done




