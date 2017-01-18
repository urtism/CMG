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

 cd $INPUTSOMATIC
  for file in *\_L001_R1_001.fastq.gz
  do
  		extract=$(echo $file| cut -d'_' -f 1,2,3,4) # Estraggo il nome del file come 20150716_01_BRCA
  		
  		printf "\n\n"
  		cat ~/Scrivania/SCRIPT_PIPELINE/logo_alignment.txt
  		printf $"\n~~~>	Sample $extract => BWA MEM\n\n"
 
  		$BWA/bwa mem $REF \
  		-M $INPUTSOMATIC/$extract\_L001_R1_001.fastq.gz \
  		$INPUTSOMATIC/$extract\_L001_R2_001.fastq.gz -t 10 \
  		> $PROCESSING/1_Align/Aligned_Sam/$extract.sam
 
  		printf $"\n~~~>	Sample $extract => BWA MEM: DONE\n\n"
  		cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
  		printf $"\n~~~>	Sample $extract => Sam Format Converter\n\n"
 
  		java -Xmx64g -jar $PICARD SamFormatConverter \
  		I=$PROCESSING/1_Align/Aligned_Sam/$extract.sam \
  		O=$PROCESSING/1_Align/Converted_bam/$extract\_Converted.bam
 
  		printf $"\n~~~>	Sample $extract => Sam Format Converter: DONE\n\n"
  		cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
  		printf $"\n~~~>	Sample $extract => Sort Sam\n\n"
 
  		java -Xmx64g -jar $PICARD SortSam \
  		I=$PROCESSING/1_Align/Converted_bam/$extract\_Converted.bam \
  		O=$PROCESSING/1_Align/Sort_bam/$extract.bam \
  		SORT_ORDER=coordinate
 
  		printf $"\n~~~>	Sample $extract => Sort Sam: DONE\n\n"
  done

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::   FUNZIONE PRE-PROCESSING     ::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 	cd $PROCESSING/1_Align/Sort_bam
 	printf $"\n"
 	cat ~/Scrivania/SCRIPT_PIPELINE/logo_processing.txt
 	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
 	
 for filename in *.bam
 do
 	printf $"\n~~~>	Sample ${filename%.*} => Add Or Replace Read Groups\n\n"
 
 	java -Xmx64g -jar $PICARD AddOrReplaceReadGroups \
 	I=$PROCESSING/1_Align/Sort_bam/${filename%.*}.bam \
 	O=$PROCESSING/2_Add/${filename%.*}.bam \
 	RGID=${filename%.*} \
 	RGPL=ILLUMINA RGSM=${filename%.*} \
 	RGPU=ILLUMINA_$DATE \
 	RGLB=CELL_FREE \
 	VALIDATION_STRINGENCY=LENIENT
 
 
 	printf $"\n~~~>	Sample ${filename%.*} => Add Or Replace Read Groups: DONE\n\n"
 	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
 	printf $"\n~~~>	Sample ${filename%.*} => Build Bam Index\n\n"
 
 	java -Xmx64g -jar $PICARD BuildBamIndex \
 	I=$PROCESSING/2_Add/${filename%.*}.bam \
 	O=$PROCESSING/2_Add/${filename%.*}.bai \
 	VALIDATION_STRINGENCY=LENIENT
 
 	printf $"\n~~~>	Sample ${filename%.*} => Build Bam Index: DONE\n\n"
 
 	printf $"\n~~~>	Sample ${filename%.*} => Indel Realigner: DONE\n\n"
 	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt

 done


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::;:::::::::::::::     FUNZIONE VARIANT CALLING     :::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


 cd $PROCESSING/2_Add/

 for filename in *Somatic.bam
 do
 	cat ~/Scrivania/SCRIPT_PIPELINE/logo_variant.txt


 	sample_name=$(echo $filename| cut -d'_' -f 1,2) # Estraggo il nome del file come 20150716_01

 	SANE=''
 	SOMATIC=$sample_name\_Somatic
 	echo "Campione in analisi: $sample_name"

 	if [[ -f $sample_name\_Sane.bam ]]
 		then
 		SANE=$sample_name\_Sane
 		#echo $SANE
 	elif [[ -f $sample_name\_BRCA.bam ]]
 		then
 		SANE=$sample_name\_BRCA
 		#echo $SANE
 	else
 		echo "Campione sano di controllo non trovato"
 	fi

 	echo -e "\nCampione Sano: $SANE \nCampione Somatico: $SOMATIC\n"
	

##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     MUTECT2     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::	
	
# 	printf $"\n =========>	Sample => Variant Calling: MuTect2\n\n"
# 	java -Xmx64g -jar $GATK -T MuTect2 \
# 	-R $REF \
# 	-I:tumor $PROCESSING/2_Add/$SOMATIC.bam \
# 	-I:normal $PROCESSING/2_Add/$SANE.bam \
# 	-L $TARGETBRCA \
# 	-o $PROCESSING/6_Variant/GATK/$sample_name\_Somatic_Sane_GATK.vcf \
# 	-bamout $PROCESSING/6_Variant/GATK/$sample_name\_Somatic_Sane.vcf.bam 
	
	
# # 	python /home/jarvis/git/CMG/scripts/header_fix.py 
# # 	
# # 	$BCFTOOLS norm -m -both \
# # 	-f $REF \
# # 	$PROCESSING/6_Variant/GATK/$sample_name\_Somatic_Sane_GATK.vcf \
# # 	> $PROCESSING/6_Variant/GATK/$sample_name\_Somatic_Sane_GATK.split.vcf	

# 	printf $"\n =========>	Sample => Variant Calling: MuTect2: DONE\n\n"
	
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARDICT     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	
# 	echo $'\n =========>	Variant Calling with VarDict\n'
	
# 	#AF_THR="0.01" # minimum allele frequency
# 	#<path_to_vardict_folder>/build/install/VarDict/bin/VarDict -G /path/to/hg19.fa -f $AF_THR -N tumor_sample_name -b "/path/to/tumor.bam|/path/to/normal.bam" -z -F -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | VarDict/testsomatic.R | VarDict/var2vcf_paired.pl -N "tumor_sample_name|normal_sample_name" -f $AF_THR

# 	$VARDICT -G $REF -f 0.01 -N $SOMATIC -b "$PROCESSING/2_Add/$SOMATIC.bam|$PROCESSING/2_Add/$SANE.bam" \
# 	-z -F 0 -c 1 -S 2 -E 3 -g 4 $TARGETBRCABED | /home/jarvis/NGS_TOOLS/VarDictJava-master/VarDict/testsomatic.R | /home/jarvis/NGS_TOOLS/VarDictJava-master/VarDict/var2vcf_paired.pl \
# 	-N "$SOMATIC|$SANE" -f 0.01 > $PROCESSING/6_Variant/VarDict/$sample_name\_Somatic_Sane_VarDict.vcf
	
# 	echo "output: $PROCESSING/6_Variant/VarDict/$sample_name\_Somatic_Sane_VarDict.vcf"
# 	printf $"\n =========>	Sample => Variant Calling: VarDict: DONE\n\n"	

# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARSCAN2     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 	echo $'\n =========>	Variant Calling with Varscan2\n\n'
	
# 	samtools mpileup -f $REF -l $TARGETBRCABED -q 1 -B $SANE.bam $SOMATIC.bam > $PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane.mpileup 
 		
# 	java -jar -Xmx64g $VARSCAN somatic $PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane.mpileup \
# 	$PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane_VarScan --output-vcf 1 --mpileup 1
	
# # 	bedtools intersect -a $PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane_VarScan.indel.vcf -b $TARGETBRCABED -wa -header \
# # 	>$PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane_VarScan.target.indel.vcf
# # 	
# # 	bedtools intersect -a $PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane_VarScan.snp.vcf -b $TARGETBRCABED -wa -header \
# # 	>$PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane_VarScan.target.snp.vcf
	
# 	printf $"\n =========>	Sample => Variant Calling: Varscan2: DONE\n\n"	
	
	

# 	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::     FEATURES EXTRACTION     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 	
#  	printf $"Estraggo le Features...\n"
#  	python /home/jarvis/git/CMG/script_CMG/features_extraction_somatic.py \
#  	-m $PROCESSING/6_Variant/GATK/$sample_name\_Somatic_Sane_GATK.vcf \
#  	-d $PROCESSING/6_Variant/VarDict/$sample_name\_Somatic_Sane_VarDict.vcf \
#  	-v $PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane_VarScan.target.snp.vcf \
#  	-i $PROCESSING/6_Variant/VarScan/$sample_name\_Somatic_Sane_VarScan.target.indel.vcf \
#  	-n $sample_name\_Sane -t $sample_name\_Somatic -a \
#  	-o $PROCESSING/6_Variant/Features_extraction/$sample_name\_Somatic
# 	printf $"Features extraction----->DONE"
 	
 	
#  	printf $"\n\n"
#  	cat ~/Scrivania/SCRIPT_PIPELINE/logo_annotation.txt
#  	printf $"\n\n\n"
#  	cd $VEP    
 
#  	perl $VEPANN -i $PROCESSING/6_Variant/Features_extraction/$sample_name\_Somatic.vcf \
#  	-o $PROCESSING/8_Annotation/$sample_name\_Somatic.ANN.vcf \
#  	--stats_file $PROCESSING/8_Annotation/$sample_name\_Somatic.ANN.html \
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
#  	--vcf
 	
 	
#  	mkdir -p $OUTVCF/$sample_name\_Somatic
 	
#  	python ~/git/CMG/script_CMG/Estrai_Annotazione_Somatic.py \
#  	-i $PROCESSING/8_Annotation/$sample_name\_Somatic.ANN.vcf \
#  	-f $PROCESSING/6_Variant/Features_extraction/$sample_name\_Somatic.features.tsv \
#  	-l $PROCESSING/8_Annotation/Annotation_list_somatic.list \
#  	-t $TRASBRCA \
#  	-o $OUTVCF/$sample_name\_SOMATIC/$sample_name\_Somatic.ANN
 	
#  	echo $'\n =========>	ANNOTATION: DONE\n'
#  	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
 	
#  	sort -V $OUTVCF/$sample_name\_Somatic/$sample_name\_Somatic.ANN.tsv > $OUTVCF/$sample_name\_Somatic/$sample_name\_Somatic.sorted.ANN.tsv
#  	sort -V $OUTVCF/$sample_name\_Somatic/$sample_name\_Somatic.ANN.Other_transcripts.tsv > $OUTVCF/$sample_name\_Somatic/$sample_name\_Somatic.sorted.ANN.Other_transcripts.tsv
 	
 	
 	
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
 	
#  done
