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
VARDICT=~/NGS_TOOLS/VarDictJava-master/build/install/VarDict/bin/VarDict
BCFTOOLS=bcftools

VEP=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/filter_vep.pl

### DATABASES & FILES ###
LISTAFEATURES_GERMLINE=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_germline.list
LISTAFEATURES_SOMATIC=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_somatic.list
ANN_LIST_GERMLINE=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_annotazione.list
ANN_LIST_SOMATIC=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_annotazione.list
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_138.hg19.vcf
LOGHI=~/git/CMG/LOGHI
TRASCR_CARDIO=~/NGS_ANALYSIS/TARGET/Lista_trascritti_cardio.txt
TRASCR_BRCA=~/NGS_ANALYSIS/TARGET/Lista_trascritti_BRCA.txt

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

HELP () {

echo -e "\n\nGenomicANALYZER v.3.0 by Matteo Digiovannantonio & Mario Urtis.\n"
echo -e "USAGE: 	bash PIPELINE_3.0.sh [OPTION(s)]\n
OPTIONS:\n
-h, --help			Print this HELP.

-a, --analisi	[String]	Type of analysis: Germline, Somatic, Cellfree.

-d, --data	[YYYYMMDD]	Data of the analysis.

-r, --run	[Int]		Number of the analysis.

-p, --pannello	[String]	Type of Panel used for the analysis: Cardio, Cancer, Exome, BRCA.
	
	--cfg	[Path]		Path to the configuration file containing the paths to necessary files.
	
	--start	[String]	String that indicates which steps the analysis must do: 
				A-> Alignment, R-> AddOrReplaceReadGroups, M-> MarkDuplicates, I-> IndelRealigner, B-> BaseRecalibrator, V-> Variant Calling, F-> Features Extraction, E-> Annotation. 
				es: --start AMIBVFE indicates all steps (like --start ALL), --start MIBV starts from MarkDuplicates and ends to Variant Calling. DEFAULT: ALL   

-w, --workdir	[Path]		Path to the directory of the analysis. If not specified a new directory will be created.

-d, --design	[String]	Design of the NGS experiment: ENRICHMENT, AMPLICON. (At the moment doesn\'t work, the Design is choosen by the Panel)\n\n"

}


PIPELINE_GERMLINE () {

	if [ "$PANNELLO" == "Cardio" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CARDIO_1000
		TARGETBED=$TARGET_CARDIO_1000_BED
		TRANSCR_LIST=$TRASCR_CARDIO
		
	elif [ "$PANNELLO" == "Cancer" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CANCER_1000
		TARGETBED=$TARGET_CANCER_1000_BED
		
	elif [ "$PANNELLO" == "Exome" ] 
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_EXOME_1000
		TARGETBED=$TARGET_EXOME_1000_BED
		
	elif [ "$PANNELLO" == "BRCA" ]
	then
		DESIGN="AMPLICON"
		TARGET=$TARGET_BRCA
		TARGETBED=$TARGET_BRCA_BED
		TRANSCR_LIST=$TRASCR_BRCA
	fi


	if [[ "$START" == *"A"* ]]
	then
		ALLINEAMENTO $CFG
	fi
	if [[ "$START" == *"R"* ]] || [[ "$START" == *"M"* ]] || [[ "$START" == *"I"* ]] || [[ "$START" == *"B"* ]]
	then
		PREPROCESSING $CFG
	fi
	if [[ "$START" == *"V"* ]]
	then
		VARIANT_CALLING_GERMLINE $CFG
	fi
	if [[ "$START" == *"F"* ]]
	then
		Features_extraction_germline $CFG
	fi
	if [[ "$START" == *"E"* ]]
	then
		INPUT=/home/jarvis/NGS_ANALYSIS_TEMP/20170118_Run_69_Cancer_5/VARIANT_CALLING/FEATURES_EXTRACTION/20170118_Cancer_TOTAL.vcf
		ANNOTATION_germline $INPUT
		#ADD_ANNOTATION_germline $INPUT
	fi	
}


PIPELINE_SOMATIC () {

	if [ "$PANNELLO" == "Cardio" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CARDIO_1000
		TARGETBED=$TARGET_CARDIO_1000_BED
		TRANSCR_LIST=$TRASCR_CARDIO
		
	elif [ "$PANNELLO" == "Cancer" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CANCER_1000
		TARGETBED=$TARGET_CANCER_1000_BED
		
	elif [ "$PANNELLO" == "Exome" ] 
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_EXOME_1000
		TARGETBED=$TARGET_EXOME_1000_BED
		
	elif [ "$PANNELLO" == "BRCA" ]
	then
		DESIGN="AMPLICON"
		TARGET=$TARGET_BRCA
		TARGETBED=$TARGET_BRCA_BED
		TRANSCR_LIST=$TRASCR_BRCA
	fi

	if [[ "$START" == *"A"* ]]
	then
		ALLINEAMENTO $CFG
	fi
	if [[ "$START" == *"R"* ]] || [[ "$START" == *"M"* ]] || [[ "$START" == *"I"* ]] || [[ "$START" == *"B"* ]]
	then
		PREPROCESSING $CFG
	fi
	if [[ "$START" == *"V"* ]]
	then
		VARIANT_CALLING_SOMATIC $CFG
	fi
	if [[ "$START" == *"F"* ]]
	then
		Features_extraction_somatic $CFG
	fi
	if [[ "$START" == *"E"* ]]
	then
		cat $CFG | while read line
		do
			INPUT=$(echo "$line" | cut -f1)
			ANNOTATION_somatic $INPUT.vcf
			ADD_ANNOTATION_somatic $INPUT
		done
	fi	
}



if [ "$#" == "0" ]
then
	echo "\nERROR: No INPUT command(s). To read the HELP type [-h] [--help] option.\n"
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
		HELP
		exit 1;
		;;
		--default)
		DEFAULT=1
		;;
		*)
		echo "ERROR: Wrong command(s). To read the HELP type [-h] [--help] option."
		exit 1;
		;;
	esac
	shift 
	done
fi

DATE=$(date +"%Y%m%d");
STARTTIME=$(date +%s)

if [ "$START" == "" ] || [ "$START" == "ALL" ]
then
	START=AMIBVFE
fi

NGSDIR=~/NGS_ANALYSIS_TEMP

var=1
while [ -d $NGSDIR/$DATA\_Run_$RUN\_$ANALISI\_$PANNELLO\_$var ]; do
    ((var++))
done

if [ "$WORKDIR" == "" ]
then
	WORKDIR=$NGSDIR/$DATA\_Run_$RUN\_$ANALISI\_$PANNELLO\_$var
fi

mkdir -p $WORKDIR


if [ "$ANALISI" == "Germline" ] || [ "$ANALISI" == "germline" ] || [ "$ANALISI" == "GERMLINE" ]
then
	PIPELINE_GERMLINE
elif [ "$ANALISI" == "Somatic" ] || [ "$ANALISI" == "somatic" ] || [ "$ANALISI" == "SOMATIC" ]
then
	PIPELINE_SOMATIC
elif [ "$ANALISI" == "Cellfree" ] || [ "$ANALISI" == "cellfree" ] || [ "$ANALISI" == "CELLFREE" ]
then
	echo "cellfree"
	#PIPELINE_CELLFREE
fi

ENDTIME=$(date +%s)
echo "TEMPO TOTALE DI PROCESSING: $(($ENDTIME - $STARTTIME)) sec"
printf "\n\n"

