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
source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/PREPROCESSING.sh
source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/VARIANT_CALLING.sh


### TOOLS ###

SCRIPT_PIPELINE=~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE
FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.15
PICARD=~/NGS_TOOLS/picard-tools-2.7.1/picard.jar

GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
FREEBAYES=~/NGS_TOOLS/freebayes/bin/freebayes
BCFTOOLS=bcftools
#VARDICT=

VEP=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/filter_vep.pl

### DATABASES & FILES ###
LISTAFEATURES_GERMLINE=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_germline.list
LISTAFEATURES_SOMATIC=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_somatic.list
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


PIPELINE_GERMLINE () {


	if [[ "$START" == *"A"* ]]
	then
		ALLINEAMENTO $CFG
		START="P"
	fi

	if [ "$PANNELLO" == "Cardio" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CARDIO_1000
		echo "Cardio"

	elif [ "$PANNELLO" == "Cancer" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CANCER_1000
		echo "Cancer"

	elif [ "$PANNELLO" == "Exome" ] 
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_EXOME_1000
		echo "Exome"

	elif [ "$PANNELLO" == "BRCA" ]
	then
		DESIGN="AMPLICON"
		TARGET=$TARGET_BRCA
		echo "BRCA"

	fi

	if [[ "$START" == *"P"* ]]
	then
		#CFG=$WORKDIR/PostAlignment.cfg
		PREPROCESSING $CFG
		START="V"
	fi
	
	if [[ "$START" == *"V"* ]]
	then

		VARIANT_CALLING_GERMLINE $CFG
	
	fi	
}

if [ "$#" == "0" ]
then
	echo "usa -h"
	exit 1;
else
	while [[ $# -gt 1 ]]
	do
	key="$1"
	case $key in
		-a|--analisi)
		ANALISI="$2"
		shift # past argument
		;;
		--start)
		START="$2"
		shift # past argument
		;;
		-w|--workdir)
		WORKDIR="$2"
		shift # past argument
		;;
		-d|--design)
		DESIGN="$2"
		shift # past argument
		;;
		--cfg)
		CFG="$2"
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
#START="A"

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

if [ "$ANALISI" == "Germline" ] || [ "$ANALISI" == "germline" ] || [ "$ANALISI" == "GERMLINE" ]
then
	PIPELINE_GERMLINE

elif [ "$ANALISI" == "Somatic" ] || [ "$ANALISI" == "somatic" ] || [ "$ANALISI" == "SOMATIC" ]
then
	echo "somatic"
	#PIPELINE_SOMATIC
elif [ "$ANALISI" == "Cellfree" ] || [ "$ANALISI" == "cellfree" ] || [ "$ANALISI" == "CELLFREE" ]
then
	echo "cellfree"
	#PIPELINE_CELLFREE
fi



