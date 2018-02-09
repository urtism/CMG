#!bin/bash

PipeVersion=3.0

### IMPORT ###
source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/ALLINEAMENTO.sh
source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/PREPROCESSING.sh
source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/VARIANT_CALLING.sh
source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/FEATURES_EXTRACTION.sh
source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/ANNOTATION.sh
#source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/VARIANT_CALLING_NOFILTERS.sh


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
SCALPEL=~/NGS_TOOLS/scalpel-0.5.3/scalpel-discovery
SNVER_INDIVIDUAL=~/NGS_TOOLS/SNVer-0.5.3/SNVerIndividual.jar
SNVER_POOL=~/NGS_TOOLS/SNVer-0.5.3/SNVerPool.jar
PLATYPUS=~/NGS_TOOLS/Platypus_0.8.1/Platypus.py

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
#DBSNP=~/NGS_TOOLS/hg19/common_all_20170710.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_138.hg19.vcf
LOGHI=~/git/CMG/LOGHI
TRASCR_CARDIO=~/NGS_ANALYSIS/TARGET/Lista_trascritti_Cardio.txt
TRASCR_BRCA=~/NGS_ANALYSIS/TARGET/Lista_trascritti_BRCA.txt
TRASCR_CANCER=~/NGS_ANALYSIS/TARGET/Lista_trascritti_Cancer.txt

### TARGET ###
TARGET_CARDIO_1000=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.list
TARGET_CARDIO_1000_BED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.bed
TARGET_CARDIO=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a.list
TARGET_CARDIO_BED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a.bed
TARGET_ILLUMINA_BRCA=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.list
TARGET_ILLUMINA_BRCA_BED=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.bed
TARGET_MULTIPLICOM_BRCA_BED=/home/jarvis/NGS_ANALYSIS/TARGET/MULTIPLICOM_BRCA/IFU423_BRCA_MASTR_Plus_SeqNext_v170133.noheader.bed
TARGET_MULTIPLICOM_BRCA=/home/jarvis/NGS_ANALYSIS/TARGET/MULTIPLICOM_BRCA/IFU423_BRCA_MASTR_Plus_SeqNext_v170133.noheader.list
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
TARGET_SURESELECT_BED=~/NGS_ANALYSIS/TARGET/Meloni_SureSelect_target_files/3066331/3066331_Covered_senza_browser.bed

AMPLICONS_HALOPLEX_BED=~/NGS_ANALYSIS/TARGET/HaloPlex_target_files/HaloPlex_Amplicons.bed

HELP () {

echo -e "\n\nGenomicANALYZER v.3.0 by Matteo Di Giovannantonio & Mario Urtis.\n"
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

CHECK_PANNELLO () {

	if [ "$PANNELLO" == "Cardio" ]
	then
		DESIGN="ENRICHMENT"
		#TARGET=$TARGET_CARDIO
		#TARGETBED=$TARGET_CARDIO_BED
		TARGET=$TARGET_CARDIO_1000
		TARGETBED=$TARGET_CARDIO_1000_BED
		TRANSCR_LIST=$TRASCR_CARDIO
		
	elif [ "$PANNELLO" == "Cancer" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CANCER_1000
		TARGETBED=$TARGET_CANCER_1000_BED
		TRANSCR_LIST=$TRASCR_CANCER
		
	elif [ "$PANNELLO" == "Exome" ] 
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_EXOME_1000
		TARGETBED=$TARGET_EXOME_1000_BED
		
	elif [ "$PANNELLO" == "BRCA" ]
	then
		DESIGN="AMPLICON"
		#TARGET=$TARGET_ILLUMINA_BRCA
		#TARGETBED=$TARGET_ILLUMINA_BRCA_BED
		TARGET=$TARGET_MULTIPLICOM_BRCA
		TARGETBED=$TARGET_MULTIPLICOM_BRCA_BED
		TRANSCR_LIST=$TRASCR_BRCA

	elif [ "$PANNELLO" == "CellFree" ]
	then
		DESIGN="AMPLICON"
		TARGET=$TARGET_CF
		TARGETBED=$TARGET_CF_BED

	elif [ "$PANNELLO" == "HaloPlex" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_HALOPLEX
		TARGETBED=$TARGET_HALOPLEX_BED

	elif [ "$PANNELLO" == "SureSelect" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_SURESELECT
		TARGETBED=$TARGET_SURESELECT_BED

	elif [ "$PANNELLO" == "" ] || [ "$PANNELLO" == "Sanger" ]
	then
		DESIGN="AMPLICON"
		TARGET=""
		TARGETBED=""
	
	fi
}

PIPELINE_GERMLINE () {

	cat $LOGHI/logo_cmg.txt

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
		#iEVA_germline $CFG
		#Features_extraction_germline_iEVA $CFG
		Features_extraction_germline $CFG
	fi
	if [[ "$START" == *"E"* ]]
	then
		cat $CFG | while read line
		do
			INPUT=$(echo "$line" | cut -f1)
			if [ "$PANNELLO" == "SureSelect" ]
				then
					ANNOTATION $INPUT CANONICAL
			elif [ "$PANNELLO" == "BRCA" ]
				then
				echo "brca"
					ANNOTATION $INPUT NOCANONICAL
					SPLIT_TRANSCRIPTS $INPUT $TRANSCR_LIST
			else
				ANNOTATION $INPUT NOCANONICAL
				SPLIT_TRANSCRIPTS $INPUT $TRANSCR_LIST
				ANNOTATION $INPUT2 CANONICAL
				MERGE_2VCF $INPUT1 $INPUT
			fi
			cat $WORKDIR/Sample_list.txt | while read line
			do
				SAMPLE=$(echo "$line" | cut -f1)
				ADD_ANNOTATION $INPUT $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE.tsv
			done
		done
		
	fi	
}


PIPELINE_SOMATIC () {

	cat $LOGHI/logo_cmg.txt

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
			if [ "$PANNELLO" == "HaloPlex" ]
			then
				ANNOTATION $INPUT.vcf CANONICAL
				ADD_ANNOTATION $INPUT ${INPUT%.*.*}.tsv
			else
				ANNOTATION $INPUT.vcf NOCANONICAL
				SPLIT_TRANSCRIPTS $INPUT $TRANSCR_LIST
				ANNOTATION $INPUT2 CANONICAL
				MERGE_2VCF $INPUT1 $INPUT
				ADD_ANNOTATION $INPUT ${INPUT%.*.*}.tsv
			fi
			
		done
	fi	
}

PIPELINE_SINGLESAMPLE () {

	cat $LOGHI/logo_cmg.txt

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
		VARIANT_CALLING_SINGLE_SAMPLE $CFG
	fi
	if [[ "$START" == *"F"* ]]
	then
		Features_extraction_germline $CFG
	fi
	if [[ "$START" == *"E"* ]]
	then
		cat $CFG | while read line
		do
			INPUT=$(echo "$line" | cut -f1)
			SAMPLE_NAME=$(echo "$line" | cut -f2)
			if [ "$PANNELLO" == "SureSelect" ]
				then
					ANNOTATION $INPUT CANONICAL
			elif [ "$PANNELLO" == "BRCA" ]
				then
					ANNOTATION $INPUT NOCANONICAL
					SPLIT_TRANSCRIPTS $INPUT $TRANSCR_LIST
			else
				ANNOTATION $INPUT NOCANONICAL
				SPLIT_TRANSCRIPTS $INPUT $TRANSCR_LIST
				ANNOTATION $INPUT2 CANONICAL
				MERGE_2VCF $INPUT1 $INPUT
			fi
			ADD_ANNOTATION $INPUT $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME.tsv

		done
		
	fi	
}

PIPELINE_CELLFREE () {

	cat $LOGHI/logo_cmg.txt

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
		VARIANT_CALLING_CELLFREE $CFG
	fi
	if [[ "$START" == *"F"* ]]
	then
		Features_extraction_cellfree $CFG
	fi
	if [[ "$START" == *"E"* ]]
	then
		cat $CFG | while read line
		do
			INPUT=$(echo "$line" | cut -f1)
			ANNOTATION $INPUT.vcf NOCANONICAL
			SPLIT_TRANSCRIPTS $INPUT $TRANSCR_LIST
			ANNOTATION $INPUT2 CANONICAL
			MERGE_2VCF $INPUT1 $INPUT
			ADD_ANNOTATION $INPUT ${INPUT%.*}.tsv
			
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
	#START=AMIB
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

STORAGE=$WORKDIR/STORAGE/$DATA\_Run_$RUN\_$ANALISI\_$PANNELLO
OUT=$WORKDIR/OUTPUT/$DATA\_Run_$RUN\_$ANALISI\_$PANNELLO\_$PipeVersion
DELETE=$WORKDIR/DELETE

mkdir -p $WORKDIR
mkdir -p $DELETE
mkdir -p $WORKDIR/STORAGE
mkdir -p $WORKDIR/OUTPUT
mkdir -p $STORAGE
mkdir -p $OUT

LOG=$WORKDIR/$DATA\_Run_$RUN\_$ANALISI\_$PANNELLO.log

CHECK_PANNELLO

cp $CFG $WORKDIR/original.cfg

if [ "$ANALISI" == "Germline" ]
then
	ANN_LIST=$ANN_LIST_GERMLINE
	PIPELINE_GERMLINE
elif [ "$ANALISI" == "Somatic" ]
then
	ANN_LIST=$ANN_LIST_SOMATIC
	PIPELINE_SOMATIC
elif [ "$ANALISI" == "CellFree" ] 
then
	PIPELINE_CELLFREE

elif [ "$ANALISI" == "SingleSample" ] 
then
	ANN_LIST=$ANN_LIST_GERMLINE
	PIPELINE_SINGLESAMPLE
fi

ENDTIME=$(date +%s)
echo "TEMPO TOTALE DI PROCESSING: $(($ENDTIME - $STARTTIME)) sec"
#printf "\nINVIO LE MAIL DI AVVISO\n"
#nome_analisi=$DATA\_Run_$RUN\_$ANALISI\_$PANNELLO\_$var
#echo -e "Ciao Valentina,\nl'analisi $nome_analisi e' completa. Speriamo sia andato tutto bene.\n\nBioinfo Alert" | mutt -s "BIOALERT" valentinafavalli@gmail.com
#echo -e "Ciao Matteo,\nl'analisi $nome_analisi e' completa. Speriamo sia andato tutto bene.\n\nBioinfo Alert" | mutt -s "BIOALERT" matteodeg@gmail.com
#echo -e "Ciao Mario,\nl'analisi $nome_analisi e' completa. Speriamo sia andato tutto bene.\n\nBioinfo Alert" | mutt -s "BIOALERT" mario.urtis01@universitadipavia.it