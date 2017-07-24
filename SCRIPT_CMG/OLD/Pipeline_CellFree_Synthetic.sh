#!/bin/bash
clear


cat ~/Scrivania/SCRIPT_PIPELINE/logo.txt 

DATA='20170111'
RUN='CF1'

FASTQC=~/NGS_TOOLS/FastQC/fastqc
#BWA=~/NGS_TOOLS/bwa-0.7.12 (aggiornato)
BWA=~/NGS_TOOLS/bwa-0.7.15
#PICARD=~/NGS_TOOLS/picard-tools-2.3.0/picard.jar (aggiornato)
PICARD=~/NGS_TOOLS/picard-tools-2.7.1/picard.jar 
GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
VARDICT=~/NGS_TOOLS/VarDictJava-master/build/install/VarDict/bin/VarDict
BCFTOOLS=~/NGS_TOOLS/bcftools-1.3.1/bcftools
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_135.hg19.vcf
INPUT=~/Scrivania/NGS_ANALYSIS_TEST/INPUT_DATA/CARDIO
INPUTBRCA=~/Scrivania/NGS_ANALYSIS_TEST/INPUT_DATA/BRCA
INPUTSOMATIC=~/Scrivania/NGS_ANALYSIS_TEST/INPUT_DATA/SOMATIC
INPUTCANCER=~/Scrivania/NGS_ANALYSIS_TEST/INPUT_DATA/CANCER
INPUTEXOME=~/Scrivania/NGS_ANALYSIS_TEST/INPUT_DATA/EXOME
PROCESSING=~/Scrivania/NGS_ANALYSIS_TEST/PROCESSING
TARGET=~/Scrivania/NGS_ANALYSIS_TEST/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.list
TARGETCELLFREEBED=~/NGS_ANALYSIS/TARGET/ctDNA_2_113416_AmpliconsExport.bed
TARGETCELLFREE=~/NGS_ANALYSIS/TARGET/ctDNA_2_113416_AmpliconsExport.list
TARGETCARDIOBED=~/Scrivania/NGS_ANALYSIS_TEST/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.bed
TARGETMETRICS=~/Scrivania/NGS_ANALYSIS_TEST/TARGET/trusight_cardio_manifest_a.list
TARGETBRCA=~/Scrivania/NGS_ANALYSIS_TEST/TARGET/AFP2_manifest_v1.list
TARGETBRCABED=~/Scrivania/NGS_ANALYSIS_TEST/TARGET/AFP2_manifest_v1.bed
TARGBRCAFREE=~/Scrivania/NGS_ANALYSIS_TEST/TARGET/BRCA_FreeBayes_amplicon.bed
TARGETEXOME=~/Scrivania/NGS_ANALYSIS_TEST/TARGET/TruSight_One_v1.1_ESTESO+-1000.list
TARGETEXOMEBED=~/Scrivania/NGS_ANALYSIS_TEST/TARGET/TruSight_One_v1.1_ESTESO+-1000.bed
TARGETCANCER=~/Scrivania/NGS_ANALYSIS_TEST/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.list
TARGETCANCERBED=~/Scrivania/NGS_ANALYSIS_TEST/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.bed
TRASBRCA=~/Scrivania/NGS_ANALYSIS_TEST/TARGET/Lista_trascritti_BRCA.txt
OUTVCF=~/Scrivania/NGS_ANALYSIS_TEST/OUTPUT_DATA
#STORAGE=~/Scrivania/NGS_ANALYSIS_TEST/STORAGE
VEP=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/filter_vep.pl
#VEP=~/NGS_TOOLS/ensembl-tools-release-84/scripts/variant_effect_predictor/ (aggiornato)
#VEPANN=~/NGS_TOOLS/ensembl-tools-release-84/scripts/variant_effect_predictor/variant_effect_predictor.pl (aggiornato)
#VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-84/scripts/variant_effect_predictor/filter_vep.pl (aggiornato)

 mkdir -p ~/Scrivania/NGS_ANALYSIS_TEST/OUTPUT_DATA/$DATA\_Run_$RUN\_Conferme_2.1
 mkdir -p ~/Scrivania/NGS_ANALYSIS_TEST/STORAGE/$DATA\_Run_$RUN\_2.1

 STORAGE=~/Scrivania/NGS_ANALYSIS_TEST/STORAGE/$DATA\_Run_$RUN\_2.1
 OUTVCF=~/Scrivania/NGS_ANALYSIS_TEST/OUTPUT_DATA/$DATA\_Run_$RUN\_Conferme_2.1
 Data=20170111

#  cd $INPUTSOMATIC
#   for file in *\_L001_R1_001.fastq.gz
#   do
#   		extract=$(echo $file| cut -d'_' -f 1,2,3,4) # Estraggo il nome del file come 20150716_01_BRCA
  		
#   		printf "\n\n"
#   		cat ~/Scrivania/SCRIPT_PIPELINE/logo_alignment.txt
#   		printf $"\n~~~>	Sample $extract => BWA MEM\n\n"
 
#   		$BWA/bwa mem $REF \
#   		-M $INPUTSOMATIC/$extract\_L001_R1_001.fastq.gz \
#   		$INPUTSOMATIC/$extract\_L001_R2_001.fastq.gz -t 10 \
#   		> $PROCESSING/1_Align/Aligned_Sam/$extract.sam
 
#   		printf $"\n~~~>	Sample $extract => BWA MEM: DONE\n\n"
#   		cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
#   		printf $"\n~~~>	Sample $extract => Sam Format Converter\n\n"
 
#   		java -Xmx64g -jar $PICARD SamFormatConverter \
#   		I=$PROCESSING/1_Align/Aligned_Sam/$extract.sam \
#   		O=$PROCESSING/1_Align/Converted_bam/$extract\_Converted.bam
 
#   		printf $"\n~~~>	Sample $extract => Sam Format Converter: DONE\n\n"
#   		cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
#   		printf $"\n~~~>	Sample $extract => Sort Sam\n\n"
 
#   		java -Xmx64g -jar $PICARD SortSam \
#   		I=$PROCESSING/1_Align/Converted_bam/$extract\_Converted.bam \
#   		O=$PROCESSING/1_Align/Sort_bam/$extract.bam \
#   		SORT_ORDER=coordinate
 
#   		printf $"\n~~~>	Sample $extract => Sort Sam: DONE\n\n"
#   done

# #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::   FUNZIONE PRE-PROCESSING     ::::::::::::::::::::::::::::::::::::::::::::::::::::::
# #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#  	cd $PROCESSING/1_Align/Sort_bam
#  	printf $"\n"
#  	cat ~/Scrivania/SCRIPT_PIPELINE/logo_processing.txt
#  	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
 	
#  for filename in *.bam
#  do
#  	printf $"\n~~~>	Sample ${filename%.*} => Add Or Replace Read Groups\n\n"
 
#  	java -Xmx64g -jar $PICARD AddOrReplaceReadGroups \
#  	I=$PROCESSING/1_Align/Sort_bam/${filename%.*}.bam \
#  	O=$PROCESSING/2_Add/${filename%.*}.bam \
#  	RGID=${filename%.*} \
#  	RGPL=ILLUMINA RGSM=${filename%.*} \
#  	RGPU=ILLUMINA_$DATE \
#  	RGLB=CELL_FREE \
#  	VALIDATION_STRINGENCY=LENIENT
 
 
#  	printf $"\n~~~>	Sample ${filename%.*} => Add Or Replace Read Groups: DONE\n\n"
#  	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
#  	printf $"\n~~~>	Sample ${filename%.*} => Build Bam Index\n\n"
 
#  	java -Xmx64g -jar $PICARD BuildBamIndex \
#  	I=$PROCESSING/2_Add/${filename%.*}.bam \
#  	O=$PROCESSING/2_Add/${filename%.*}.bai \
#  	VALIDATION_STRINGENCY=LENIENT
 
#  	printf $"\n~~~>	Sample ${filename%.*} => Build Bam Index: DONE\n\n"
 
#  	printf $"\n~~~>	Sample ${filename%.*} => Indel Realigner: DONE\n\n"
#  	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt

#  done


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::;:::::::::::::::     FUNZIONE VARIANT CALLING     :::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
paz=01

cd $PROCESSING/2_Add/

for a in 0 5 10 25 50 100
do
	b=$(expr 100 - $a)

	TUMOR=$Data\_$paz\_CF_$a\D$b.bam
	NORMAL=$Data\_$paz\_G_Rec.bam
	TUMOR_NAME=$Data\_$paz\_CF_$a\D$b
	NORMAL_NAME=$Data\_$paz\_G_Rec
 	Data=20170111

# 	printf $"\n =========>	Sample => Variant Calling: MuTect2 per T1\n\n"
#  	java -Xmx64g -jar $GATK -T MuTect2 \
#  	-R $REF \
#  	-I:tumor $PROCESSING/2_Add/$TUMOR \
#  	-I:normal $PROCESSING/2_Add/$NORMAL \
#  	-L $TARGETCELLFREE \
#  	-o $PROCESSING/6_Variant/GATK/$Data\_$paz\_G_$a\D$b\_GATK.vcf \
#  	-bamout $STORAGE/$Data\_$paz\_G_$a\D$b.vcf.bam


# 	echo $'\n =========>	Variant Calling with VarDict per $paz \n'
	
# 	$VARDICT -G $REF -f 0.01 -N $TUMOR_NAME -b "$TUMOR|$NORMAL" \
# 	-z -F 0 -c 1 -S 2 -E 3 -g 4 $TARGETCELLFREEBED | /home/jarvis/NGS_TOOLS/VarDictJava-master/VarDict/testsomatic.R | /home/jarvis/NGS_TOOLS/VarDictJava-master/VarDict/var2vcf_paired.pl \
# 	-N "$TUMOR_NAME|$NORMAL_NAME" -f 0.01 > $PROCESSING/6_Variant/VarDict/$Data\_$paz\_G_$a\D$b\_Vardict.vcf


# 	echo $'\n =========>	Variant Calling with VarScan per T1\n'

# 	samtools mpileup -f $REF -l $TARGETCELLFREEBED -q 1 -d 50000 -L 50000 -B $PROCESSING/2_Add/$NORMAL $PROCESSING/2_Add/$TUMOR \
# 	> $PROCESSING/6_Variant/VarScan/$Data\_$paz\_G_$a\D$b.mpileup
 		
# 	java -jar -Xmx64g $VARSCAN somatic $PROCESSING/6_Variant/VarScan/$Data\_$paz\_G_$a\D$b.mpileup \
# 	$PROCESSING/6_Variant/VarScan/$Data\_$paz\_G_$a\D$b\_VarScan --output-vcf 1 --mpileup 1

# done

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     NORM     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# for (( a=1; a<5; a++ ))
# do

#:::::::::::::::::::MUTECT:::::::::::::::::::
	echo "splitto GATK"

	$BCFTOOLS norm -m -both \
	-f $REF \
	$PROCESSING/6_Variant/GATK/$Data\_$paz\_G_$a\D$b\_GATK.vcf \
	> $PROCESSING/6_Variant/GATK/NORM_GATK/$Data\_$paz\_G_$a\D$b\_GATK.split.vcf



#:::::::::::::::::::VARDICT:::::::::::::::::::
	echo "splitto VARDICT"

	$BCFTOOLS norm -m -both \
	-f $REF \
	$PROCESSING/6_Variant/VarDict/$Data\_$paz\_G_$a\D$b\_Vardict.vcf \
	> $PROCESSING/6_Variant/VarDict/NORM_VARDICT/$Data\_$paz\_G_$a\D$b\_Vardict.split.vcf


#:::::::::::::::::::VARSCAN:::::::::::::::::::
	echo "splitto VARSCAN"

	$BCFTOOLS norm -m -both \
	-f $REF \
	$PROCESSING/6_Variant/VarScan/$Data\_$paz\_G_$a\D$b\_VarScan.snp.vcf \
	>$PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$Data\_$paz\_G_$a\D$b\_VarScan.snp.split.vcf

	$BCFTOOLS norm -m -both \
	-f $REF \
	$PROCESSING/6_Variant/VarScan/$Data\_$paz\_G_$a\D$b\_VarScan.indel.vcf \
	>$PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$Data\_$paz\_G_$a\D$b\_VarScan.indel.split.vcf

	#done




# 	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::     FEATURES EXTRACTION     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 	
 	printf $"Estraggo le Features per il campione $TUMOR_NAME\n"

 	python /home/jarvis/git/CMG/script_CMG/features_extraction_somatic.py \
 	-m $PROCESSING/6_Variant/GATK/NORM_GATK/$Data\_$paz\_G_$a\D$b\_GATK.split.vcf \
 	-d $PROCESSING/6_Variant/VarDict/NORM_VARDICT/$Data\_$paz\_G_$a\D$b\_Vardict.split.vcf \
 	-v $PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$Data\_$paz\_G_$a\D$b\_VarScan.snp.split.vcf \
 	-i $PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$Data\_$paz\_G_$a\D$b\_VarScan.indel.split.vcf \
 	-n $NORMAL_NAME -t $TUMOR_NAME -a \
 	-o $PROCESSING/6_Variant/Features_extraction/$Data\_$paz\_G_$a\D$b
	printf $"Features extraction----->DONE\n"

#  	printf $"\n\n"
#  	cat ~/Scrivania/SCRIPT_PIPELINE/logo_annotation.txt
#  	printf $"\n\n\n"


#  	cd $VEP    
 
#  	perl $VEPANN -i $PROCESSING/6_Variant/Features_extraction/$Data\_$paz\_G_$a\D$b \
#  	-o $PROCESSING/8_Annotation/$Data\_$paz\_G_$a\D$b.ANN.vcf \
#  	--stats_file $PROCESSING/8_Annotation/$Data\_$paz\_G_$a\D$b.ANN.html \
#  	--cache \
#  	--assembly GRCh37 \
#  	--offline \
#  	--force_overwrite \
#  	-v \
#  	--fork 10 \
#  	--variant_class \
#  	--sift b \
#  	--poly b \
#  	--vcf_info_field ANN \
#  	--hgvs \
#  	--protein \
#  	--canonical \
#  	--check_existing \
#  	--gmaf \
#  	--pubmed \
#  	--species homo_sapiens \
#  	--failed 1 \
#  	--vcf \
#	--pick

#	perl $VEPANN -i $PROCESSING/6_Variant/Features_extraction/$Data\_$paz\_G_CFT2 \
#  	-o $PROCESSING/8_Annotation/$Data\_$paz\_G_CFT2.ANN.vcf \
#  	--stats_file $PROCESSING/8_Annotation/$Data\_$paz\_G_CFT2.ANN.html \
#  	--cache \
#  	--assembly GRCh37 \
#  	--offline \
#  	--force_overwrite \
#  	-v \
#  	--fork 10 \
#  	--variant_class \
#  	--sift b \
#  	--poly b \
#  	--vcf_info_field ANN \
#  	--hgvs \
#  	--protein \
#  	--canonical \
#  	--check_existing \
#  	--gmaf \
#  	--pubmed \
#  	--species homo_sapiens \
#  	--failed 1 \
#  	--vcf \
#	--pick
 	
 	
#  	mkdir -p $OUTVCF/$Data\_$paz\_G_$a\D$b
#	mkdir -p $OUTVCF/$Data\_$paz\_G_CFT2
 	
#  	python ~/git/CMG/script_CMG/Estrai_Annotazione_Somatic.py \
#  	-i $PROCESSING/8_Annotation/$Data\_$paz\_G_$a\D$b.ANN.vcf \
#  	-f $PROCESSING/6_Variant/Features_extraction/$Data\_$paz\_G_$a\D$b.features.tsv \
#  	-l $PROCESSING/8_Annotation/Annotation_list_somatic.list \
#  	-t $TRASBRCA \
#  	-o $OUTVCF/$sample_name\_SOMATIC/$Data\_$paz\_G_$a\D$b.ANN

#  	python ~/git/CMG/script_CMG/Estrai_Annotazione_Somatic.py \
#  	-i $PROCESSING/8_Annotation/$Data\_$paz\_G_CFT2.ANN.vcf \
#  	-f $PROCESSING/6_Variant/Features_extraction/$Data\_$paz\_G_CFT2.features.tsv \
#  	-l $PROCESSING/8_Annotation/Annotation_list_somatic.list \
#  	-t $TRASBRCA \
#  	-o $OUTVCF/$sample_name\_SOMATIC/$Data\_$paz\_G_CFT2.ANN
 	
#  	echo $'\n =========>	ANNOTATION: DONE\n'
#  	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
 	
#  	sort -V $OUTVCF/$Data\_$paz\_G_$a\D$b/$Data\_$paz\_G_$a\D$b.ANN.tsv > $OUTVCF/$Data\_$paz\_G_$a\D$b/$Data\_$paz\_G_$a\D$b.sorted.ANN.tsv
#  	sort -V $OUTVCF/$Data\_$paz\_G_$a\D$b/$Data\_$paz\_G_$a\D$b.ANN.Other_transcripts.tsv > $OUTVCF/$Data\_$paz\_G_$a\D$b/$Data\_$paz\_G_$a\D$b.sorted.ANN.Other_transcripts.tsv
 
#  	sort -V $OUTVCF/$Data\_$paz\_G_CFT2/$Data\_$paz\_G_CFT2.ANN.tsv > $OUTVCF/$Data\_$paz\_G_CFT2/$Data\_$paz\_G_CFT2.sorted.ANN.tsv
#  	sort -V $OUTVCF/$Data\_$paz\_G_CFT2/$Data\_$paz\_G_CFT2.ANN.Other_transcripts.tsv > $OUTVCF/$Data\_$paz\_G_CFT2/$Data\_$paz\_G_CFT2.sorted.ANN.Other_transcripts.tsv
 
 	
 	
#  	mv $PROCESSING/2_Add/$SOMATIC.bam $STORAGE
# 	mv $PROCESSING/2_Add/$SOMATIC.bai $STORAGE
# 	mv $PROCESSING/2_Add/$SANE.bam $STORAGE
# 	mv $PROCESSING/2_Add/$SANE.bai $STORAGE
#  	mv $PROCESSING/6_Variant/GATK/$sample_name\_Somatic_Sane_GATK.vcf $STORAGE
#  	mv $PROCESSING/6_Variant/GATK/$sample_name\_Somatic_Sane_GATK.vcf.idx $STORAGE
#  	mv $PROCESSING/6_Variant/GATK/$sample_name\_Somatic_Sane.vcf.bam $STORAGE
#  	mv $PROCESSING/6_Variant/GATK/$sample_name\_Somatic_Sane.vcf.bai $STORAGE
#  	mv $PROCESSING/6_Variant/VarDict/$sample_name\_Somatic_Sane_VarDict.vcf $STORAGE
#  	mv $PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane_VarScan.target.snp.vcf $STORAGE
#  	mv $PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane_VarScan.target.indel.vcf $STORAGE
#  	rm $PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane.mpileup
#  	rm $PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane_VarScan.indel.vcf
#  	rm $PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane_VarScan.snp.vcf
 	
  done