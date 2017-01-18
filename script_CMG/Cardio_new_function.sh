#!/bin/bash

cat ~/SCRIPT_PIPELINE/logo_cmg.txt

FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.15
PICARD=~/NGS_TOOLS/picard-tools-2.7.1/picard.jar
GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
BCFTOOLS=~/NGS_TOOLS/bcftools-1.3.1/bcftools
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_138.hg19.vcf
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
TARGETCARDIOPARALLEL=~/NGS_ANALYSIS/TARGET/Parallel_regions_Cardio+-1000.regions
DIAGNOSE=~/NGS_ANALYSIS/PROCESSING/9_Diagnose
OUTVCF=~/NGS_ANALYSIS/OUTPUT_DATA
STORAGE=~/NGS_ANALYSIS/STORAGE
VEP=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/filter_vep.pl

DataRun=20151113


cd $PROCESSING/5_BQSR/

		ls $PROCESSING/5_BQSR/*.bam > $PROCESSING/5_BQSR/Path_bam.list

		java -jar -Xmx64g $GATK -T DiagnoseTargets \
		-R ~/NGS_TOOLS/hg19/ucsc.hg19.fasta \
		-I $PROCESSING/5_BQSR/Path_bam.list \
		-L ~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a.list \
		-o $DIAGNOSE/$DataRun\_Target_Diagnosis.vcf \
		--missing_intervals $DIAGNOSE/$DataRun\_Missing.tsv


cd $PROCESSING/5_BQSR/

	for filename in *.bam
	do
 
 		cat ~/SCRIPT_PIPELINE/logo_variant.txt
 		printf $"\n =========>	Sample ${filename%.*} => Variant Calling: Haplotype Caller\n\n"
 
 		java -Xmx64g -jar $GATK -T HaplotypeCaller \
 		-R $REF \
 		-I ${filename%.*}.bam \
 		-o $PROCESSING/6_Variant/GATK/${filename%.*}.g.vcf \
 		-ERC GVCF \
 		--doNotRunPhysicalPhasing \
 		--heterozygosity 0.001 \
 		--indel_heterozygosity 1.25E-4 \
 		--maxReadsInRegionPerSample 50000 \
		--min_base_quality_score 10 \
		--minReadsPerAlignmentStart 1 \
 		--max_alternate_alleles 6 \
 		--standard_min_confidence_threshold_for_calling 1.0 \
		-L $TARGET
#		-bamout $PROCESSING/6_Variant/GATK/${filename%.*}.g.vcf.bam \

 		printf $"\n =========>	Sample ${filename%.*} => Variant Calling: Haplotype Caller: DONE\n\n"
 
	done




#::::::::::::::::::::::::::::::::::::::::::     INIZIO VARIANT CALLING MULTI-SAMPLE     :::::::::::::::::::::::::::::::::::::::::::::::




	#Stampo tutti i file .g.vcf in una lista
	
	cd $PROCESSING/6_Variant/GATK/
	
	ls *.g.vcf > gvcf.list
	cat ~/SCRIPT_PIPELINE/logo_multi.txt
	echo $'\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling\n\n'

		java -jar -Xmx64g $GATK -T GenotypeGVCFs \
		-R $REF \
		-V:VCF $PROCESSING/6_Variant/GATK/gvcf.list \
		-o $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK.vcf


	# Salvo nella cartella storage
	# cp $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK.vcf $STORAGE/$Name_Dir
	
	echo $'\n\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling: DONE\n'
	cat ~/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	GATK vcf norm\n\n'

	python ~/git/CMG/scripts/header_fix.py \
	-f $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK.vcf \
	-v G \
	> $PROCESSING/6_Variant/GATK/NORM_GATK/$DataRun\_Cardio_GATK.fixed.vcf
		
	bcftools norm -m -both \
	-f $REF \
	$PROCESSING/6_Variant/GATK/NORM_GATK/$DataRun\_Cardio_GATK.fixed.vcf \
	> $PROCESSING/6_Variant/GATK/NORM_GATK/$DataRun\_Cardio_GATK.split.vcf	

	echo $'\n =========>	GATK vcf norm: DONE\n'
	cat ~/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	Variant Calling with VarScan2\n\n'




#::::::::::::::::::::::::::::::::::::::::::::::::::::::::     FILTRAGGIO GATK     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



	
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     FREEBAYES     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::




	cd $PROCESSING/5_BQSR/

	ls *.bam > Bam_list.txt
	ls *.bam > Sample_list.txt
	sed -i -e "s/.bam//g" Sample_list.txt

	~/NGS_TOOLS/freebayes/scripts/freebayes-parallel $TARGETCARDIOPARALLEL 10 -f $REF \
	-L $PROCESSING/5_BQSR/Bam_list.txt \
	-K \
	-J \
	-s $PROCESSING/5_BQSR/Sample_list.txt \
	--pooled-continuous \
	--pooled-discrete \
	--genotype-qualities \
	--report-genotype-likelihood-max \
	--allele-balance-priors-off \
	--min-mapping-quality 1 \
	--min-base-quality 1 \
	--mismatch-base-quality-threshold 1 \
	--min-alternate-fraction 0.001 \
	--min-alternate-count 1 \
	--min-coverage 10 > $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes.vcf

	echo $'\n =========>	Variant Calling with FreeBayes: DONE\n'
	cat ~/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	FreeBayes vcf norm\n\n'

	python ~/git/CMG/scripts/header_fix.py \
	-f $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes.vcf \
	-v F \
	> $PROCESSING/6_Variant/FreeBayes/NORM_FREE/$DataRun\_Cardio_FreeBayes.fixed.vcf
		
	bcftools norm -m -both \
	-f $REF \
	$PROCESSING/6_Variant/FreeBayes/NORM_FREE/$DataRun\_Cardio_FreeBayes.fixed.vcf \
	> $PROCESSING/6_Variant/FreeBayes/NORM_FREE/$DataRun\_Cardio_FreeBayes.split.vcf	

	echo $'\n =========>	FreeBayes vcf norm: DONE\n'
	cat ~/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	Variant Calling with VarScan2\n\n'





#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARSCAN2     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



 	samtools mpileup -B -q 1 -f --max-depth 30000 -L 30000 --min-BQ 1 $REF \
 	-l $TARGETCARDIOBED -b $PROCESSING/5_BQSR/Bam_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup
 
 	java -jar -Xmx64g $VARSCAN mpileup2snp $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup \
 	--min-coverage 1 \
 	--min-var-freq 0.001 \
 	--min-reads2 1 \
 	--min-avg-qual 1 \
 	--strand-filter 0 \
 	--output-vcf 1 \
 	--vcf-sample-list $PROCESSING/5_BQSR/Sample_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf
 
 	java -jar -Xmx64g $VARSCAN mpileup2indel $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup \
 	--min-coverage 1 \
 	--min-var-freq 0.001 \
 	--min-reads2 1 \
 	--min-avg-qual 1 \
 	--strand-filter 0 \
 	--output-vcf 1 \
 	--vcf-sample-list $PROCESSING/5_BQSR/Sample_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf
 
 	echo $'\n =========>	Variant Calling with VarScan2: DONE\n'
 	cat ~/SCRIPT_PIPELINE/logo_cornice.txt
 	echo $'\n =========>	VarScan2 vcf merge and sort\n\n'
 
 	#Copio in Intersect
 	bgzip $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf
 	bgzip $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf
 	tabix $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf.gz
 	tabix $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf.gz
 
 	#Effettuo il merging dei file
 	
 	vcf-concat $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf.gz \
 	$PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf.gz > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Merge.vcf
 	
 	vcf-sort -c $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Merge.vcf > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Merge_Sort.vcf
 	mv $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Merge_Sort.vcf $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan.vcf
 	rm $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Merge.vcf
 	rm $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf.gz.tbi
 	rm $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf.gz.tbi
 	rm $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf.gz
 	rm $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf.gz
 	
 	echo $'\n =========>	VarScan2 vcf merge and sort: DONE\n'
	cat ~/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	FreeBayes vcf norm\n\n'

	python ~/git/CMG/scripts/header_fix.py \
	-f $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan.vcf \
	-v V \
	> $PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$DataRun\_Cardio_VarScan.fixed.vcf
		
	bcftools norm -m -both \
	-f $REF \
	$PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$DataRun\_Cardio_VarScan.fixed.vcf \
	> $PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$DataRun\_Cardio_VarScan.split.vcf	

	echo $'\n =========>	FreeBayes vcf norm: DONE\n'
	cat ~/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	Variant Calling with VarScan2\n\n'


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::    GATK ANNOTATION     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	#Partire da questi
	#$PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$DataRun\_Cardio_VarScan.split.vcf
	#$PROCESSING/6_Variant/FreeBayes/NORM_FREE/$DataRun\_Cardio_FreeBayes.split.vcf
	#$PROCESSING/6_Variant/GATK/NORM_GATK/$DataRun\_Cardio_GATK.split.vcf	



# 	printf "\n\n"
# 	cat ~/SCRIPT_PIPELINE/logo_annotation.txt
# 	printf $"\n\n\n"
 
# 	cd $VEP
# 
# 	perl $VEPANN -i $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK_Filter.vcf \
# 	-o $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_ANN.vcf \
# 	--stats_file $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_ANN.html \
# 	--cache \
# 	--assembly GRCh37 \
# 	--offline \
# 	--force_overwrite \
# 	-v \
# 	--fork 10 \
# 	--variant_class \
# 	--sift b \
# 	--poly b \
# 	--vcf_info_field ANN \
# 	--hgvs \
# 	--protein \
# 	--canonical \
# 	--check_existing \
# 	--gmaf \
# 	--pubmed \
# 	--species homo_sapiens \
# 	--failed 1 \
# 	--vcf
# 
# 	echo $'\n =========>	GATK ANNOTATION: DONE\n'
# 	cat ~/SCRIPT_PIPELINE/logo_cornice.txt
# 	echo $'\n =========>	INTERSECTION ANNOTATION\n\n'
# 
# 	rm $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK_Filter.vcf
# 	rm $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK_Filter.vcf.idx
# 	cp $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_ANN.vcf $STORAGE/$Name_Dir
# 	printf "\n${DataRun}_Cardio_GATK_ANN.vcf : contiene il file ${DataRun}_Cardio_GATK_Filter.vcf annotato\n" >> $STORAGE/$Name_Dir/README.txt
# 	mv $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_ANN.html $STORAGE/$Name_Dir
# 	printf "\n${DataRun}_Cardio_GATK_ANN.html : contiene le statistiche del file ${DataRun}_Cardio_GATK_Filter.vcf\n" >> $STORAGE/$Name_Dir/README.txt
# 
# 
# 
# #:::::::::::::::::::::::::::::::::::::::::::::::::::    INTERSECTION  ANNOTATION     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 
# 
# 
# 	cd $VEP
# 
# 	perl $VEPANN -i $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf \
# 	-o $PROCESSING/8_Annotation/$DataRun\_Cardio_Total_Intersect_ANN.vcf \
# 	--stats_file $PROCESSING/8_Annotation/$DataRun\_Cardio_Total_Intersect_ANN.html \
# 	--cache \
# 	--assembly GRCh37 \
# 	--offline \
# 	--force_overwrite \
# 	-v \
# 	--fork 10 \
# 	--variant_class \
# 	--sift b \
# 	--poly b \
# 	--vcf_info_field ANN \
# 	--hgvs \
# 	--protein \
# 	--canonical \
# 	--check_existing \
# 	--gmaf \
# 	--pubmed \
# 	--species homo_sapiens \
# 	--failed 1 \
# 	--vcf
# 
# 	rm $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf
# 	cp $PROCESSING/8_Annotation/$DataRun\_Cardio_Total_Intersect_ANN.vcf $STORAGE/$Name_Dir/$DataRun\_Cardio_Intersect_ANN.vcf
# 	printf "\n${DataRun}_Cardio_Intersect_ANN.vcf : contiene il file ${DataRun}_Cardio_Intersect.vcf annotato\n" >> $STORAGE/$Name_Dir/README.txt
# 	mv $PROCESSING/8_Annotation/$DataRun\_Cardio_Total_Intersect_ANN.html $STORAGE/$Name_Dir
# 	printf "\n${DataRun}_Cardio_Intersect_ANN.html : contiene le statistiche del file ${DataRun}_Cardio_Intersect_ANN.vcf\n" >> \
# 	$STORAGE/$Name_Dir/README.txt
# 
# 	echo $'\n =========>	INTERSECTION ANNOTATION: DONE\n'
# 	cat ~/SCRIPT_PIPELINE/logo_cornice.txt
# 	echo $'\n =========>	FILE MODIFICATION FOR CSV FORMAT \n\n'

