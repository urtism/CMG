#!bin/bash

PipeVersion=3.0

### IMPORT ###
source ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/VARIANT_CALLING.sh



### TOOLS ###
SCRIPT_PIPELINE=~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE

GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
FREEBAYES=~/NGS_TOOLS/freebayes/bin/freebayes
BCFTOOLS=bcftools
SCALPEL=~/NGS_TOOLS/scalpel-0.5.3/scalpel-discovery
PLATYPUS=~/NGS_TOOLS/Platypus_0.8.1/Platypus.py
SNVER=~/NGS_TOOLS/SNVer-0.5.3/SNVerIndividual.jar

### DATABASES & FILES ###
#LISTAFEATURES_SOMATIC=/home/jarvis/NGS_ANALYSIS/TARGET/Features_lists/lista_features_somatic_CF_20170315.list
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
REF_DICT=~/NGS_TOOLS/hg19/ucsc.hg19.dict
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_138.hg19.vcf
LOGHI=~/git/CMG/LOGHI
VARIANT_CALLING=~/NGS_ANALYSIS/PROCESSING/6_Variant/

### TARGET ###
TARGET_Cardio_1000=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.list
TARGET_Cardio_1000_BED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.bed
TARGET_Cardio_1000_TXT=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.txt
TARGET_Cancer_1000=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.list
TARGET_Cancer_1000_BED=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.bed

HELP () {

echo -e "\n\nGenomicANALYZER v.3.0 by Matteo Di Giovannantonio & Mario Urtis.\n"
echo -e "USAGE: 	bash PIPELINE_3.0.sh [OPTION(s)]\n
OPTIONS:\n

-h, --help			Print this HELP.

-list, --lista	[String]	Path to list of folder containing bam file to analyze. One path per row. Each row represents a run\n\n

-O, --output [String]	Path to folder where report the run folder containing vcf files\n\n"

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
		-l|--list)
		FOLDERLIST="$2"
		shift # past argument
		;;
		-O|--ouput)
		OUTDIR="$2"
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

while read p;
do
	RUN=$(echo $p| rev | cut -d'/' -f 1 | rev)
	PANNELLO=$(echo $p| rev | cut -d'/' -f 1 | rev | cut -d'_' -f 4)
	DATA=$(echo $p| rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1)

  	cd $p

  	echo $p
  	echo $RUN
  	echo $PANNELLO
  	echo $DATA

  	ls *.bam > $DATA\_Bamlist.list

	BAMLIST=$DATA\_Bamlist.list

  	if [ "$PANNELLO" == "Cardio" ]

  		then

  		BED=$TARGET_Cardio_1000_BED
  		LIST=$TARGET_Cardio_1000
  		PLATXT=$TARGET_Cardio_1000_TXT

  		mkdir $OUTDIR/$RUN\_VCF\_IEVA

  		FILENAME=$DATA\_$PANNELLO

  		RUN_OUT=$OUTDIR/$RUN\_VCF\_IEVA

  		cd $p

		#----------- SAMTOOLS ----------------
		
		samtools mpileup \
		-uf $REF \
		-l $BED \
		-d 50000 \
		-L 50000 \
		-b $BAMLIST | bcftools call -mv -Oz -f GQ,GP --threads 4 > $VARIANT_CALLING/Samtools/$FILENAME\_Samtools.vcf.gz

		gunzip $VARIANT_CALLING/Samtools/$FILENAME\_Samtools.vcf.gz

		python $SCRIPT_PIPELINE/header_fix.py \
		-f $VARIANT_CALLING/Samtools/$FILENAME\_Samtools.vcf \
		-v S > $VARIANT_CALLING/Samtools/$FILENAME\_Samtools.fix.vcf

		$BCFTOOLS norm -m -both \
		-f $REF \
		$VARIANT_CALLING/Samtools/$FILENAME\_Samtools.fix.vcf \
 		> $VARIANT_CALLING/Samtools/$FILENAME\_Samtools.norm.vcf

		rm -f $VARIANT_CALLING/Samtools/$FILENAME\_Samtools.fix.vcf
		rm -f $VARIANT_CALLING/Samtools/$FILENAME\_Samtools.vcf

		mv $VARIANT_CALLING/Samtools/$FILENAME\_Samtools.norm.vcf $RUN_OUT/

 		#---------- PLATYPUS ----------------

 	 	python $PLATYPUS callVariants --bamFiles=$BAMLIST \
 	 	--refFile=$REF \
 	 	--regions=$PLATXT \
 	 	--nCPU=6 \
 	 	--maxSize=5000 \
 	 	--output=$VARIANT_CALLING/Platypus/$FILENAME\_Platypus_Dup.vcf

 	 	#Rimuovo righe duplicate
	 	awk '!a[$0]++' $VARIANT_CALLING/Platypus/$FILENAME\_Platypus_Dup.vcf > $VARIANT_CALLING/Platypus/$FILENAME\_Platypus.vcf

	 	rm -f $VARIANT_CALLING/Platypus/$FILENAME\_Platypus_Dup.vcf

	 	python $SCRIPT_PIPELINE/header_fix.py \
	 	-f $VARIANT_CALLING/Platypus/$FILENAME\_Platypus.vcf \
	 	-v P > $VARIANT_CALLING/Platypus/$FILENAME\_Platypus.fix.vcf

	 	$BCFTOOLS norm -m -both \
		-f $REF \
		$VARIANT_CALLING/Platypus/$FILENAME\_Platypus.fix.vcf \
 		> $VARIANT_CALLING/Platypus/$FILENAME\_Platypus.norm.vcf

		mv $VARIANT_CALLING/Platypus/$FILENAME\_Platypus.norm.vcf $RUN_OUT/

		#-------------------------- SINGLE SAMPLE SNVer -------------------------

		# for file in *.bam

  # 			do

  # 				SAMPLE_NAME=$(echo $file| cut -d'.' -f 1 | cut -d'_' -f1-3 )

  # 				echo $SAMPLE_NAME

		#   		java -jar -Xmx8g $SNVER \
		# 		-i $file \
		# 		-l $BED \
		# 		-r $REF \
		# 		-o $VARIANT_CALLING/SNVer/$SAMPLE_NAME\_SNVer

		# 		rm -f ~/NGS_ANALYSIS/PROCESSING/6_Variant/SNVer/20160426_03_Conn_SNVer.failed.log

		# 		bgzip $VARIANT_CALLING/SNVer/$SAMPLE_NAME\_SNVer.filter.vcf
		# 		bgzip $VARIANT_CALLING/SNVer/$SAMPLE_NAME\_SNVer.indel.filter.vcf
		# 		tabix $VARIANT_CALLING/SNVer/$SAMPLE_NAME\_SNVer.filter.vcf.gz
		# 		tabix $VARIANT_CALLING/SNVer/$SAMPLE_NAME\_SNVer.indel.filter.vcf.gz

		# 		vcf-concat $VARIANT_CALLING/SNVer/$SAMPLE_NAME\_SNVer.filter.vcf.gz $VARIANT_CALLING/SNVer/$SAMPLE_NAME\_SNVer.indel.filter.vcf.gz > $VARIANT_CALLING/SNVer/$SAMPLE_NAME\_SNVer.concat.vcf
		# 		vcf-sort -c $VARIANT_CALLING/SNVer/$SAMPLE_NAME\_SNVer.concat.vcf > $VARIANT_CALLING/SNVer/$SAMPLE_NAME\_SNVer.sort.vcf

		# 		rm -f ~/NGS_ANALYSIS/PROCESSING/6_Variant/SNVer/$SAMPLE_NAME\_SNVer.filter.vcf.gz
		# 		rm -f ~/NGS_ANALYSIS/PROCESSING/6_Variant/SNVer/$SAMPLE_NAME\_SNVer.filter.vcf.gz.tbi
		# 		rm -f ~/NGS_ANALYSIS/PROCESSING/6_Variant/SNVer/$SAMPLE_NAME\_SNVer.indel.filter.vcf.gz
		# 		rm -f ~/NGS_ANALYSIS/PROCESSING/6_Variant/SNVer/$SAMPLE_NAME\_SNVer.indel.filter.vcf.gz.tbi
		# 		rm -f ~/NGS_ANALYSIS/PROCESSING/6_Variant/SNVer/$SAMPLE_NAME\_SNVer.failed.log
		# 		rm -f ~/NGS_ANALYSIS/PROCESSING/6_Variant/SNVer/$SAMPLE_NAME\_SNVer.raw.vcf
		# 		rm -f ~/NGS_ANALYSIS/PROCESSING/6_Variant/SNVer/$SAMPLE_NAME\_SNVer.indel.raw.vcf

		# 		$BCFTOOLS norm -m -both \
		# 		-f $REF \
		# 		$VARIANT_CALLING/SNVer/$SAMPLE_NAME\_SNVer.sort.vcf \
 	# 			> $RUN_OUT/$SAMPLE_NAME\_SNVer.norm.vcf

 	# 			rm -f ~/NGS_ANALYSIS/PROCESSING/6_Variant/SNVer/$SAMPLE_NAME\_SNVer.sort.vcf
  # 				rm -f ~/NGS_ANALYSIS/PROCESSING/6_Variant/SNVer/$SAMPLE_NAME\_SNVer.concat.vcf

		# 	done

		#VCF_SCALPEL=$WORKDIR/VARIANT_CALLING/$2\_Scalpel.vcf


	  	#for file in *.bam
	  	#	do
	  	#		echo $file
	  	#	done
	  	fi
done < $FOLDERLIST

ENDTIME=$(date +%s)
echo "TEMPO TOTALE DI PROCESSING: $(($ENDTIME - $STARTTIME)) sec"

