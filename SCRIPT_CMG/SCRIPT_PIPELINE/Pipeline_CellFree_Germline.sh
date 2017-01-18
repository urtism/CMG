#!/bin/bash
clear


cat ~/Scrivania/SCRIPT_PIPELINE/logo.txt 

DATA='20170111'
RUN='CF1_Germline'

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
LISTAFEATURES=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_germline.list
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




#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::;:::::::::::::::     FUNZIONE VARIANT CALLING     :::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


cd $PROCESSING/2_Add/

# for (( a=1; a<5; a++ ))
# do
#   		Data=20170111


# 		cat ~/Scrivania/SCRIPT_PIPELINE/logo_variant.txt
# 		printf $"\n =========>	Sample $a Rec => Variant Calling: Haplotype Caller\n\n"

# 		java -Xmx64g -jar $GATK -T HaplotypeCaller \
# 		-R $REF \
# 		-I $PROCESSING/2_Add/$Data\_0$a\_G_Rec.bam \
# 		-o $PROCESSING/6_Variant/GATK/$Data\_0$a\_G_Rec.g.vcf \
# 		-ERC GVCF \
# 		--doNotRunPhysicalPhasing \
#  		--heterozygosity 0.001 \
#  		--indel_heterozygosity 1.25E-4 \
#  		--maxReadsInRegionPerSample 50000 \
# 		--min_base_quality_score 10 \
# 		--minReadsPerAlignmentStart 5 \
#  		--max_alternate_alleles 6 \
# 		-bamout $STORAGE/$Data\_0$a\_G_Rec.g.vcf.bam \
# 		-L $TARGETCELLFREE

# 		printf $"\n =========>	Sample $a Rec => Variant Calling: Done\n\n"
# 		cat ~/Scrivania/SCRIPT_PIPELINE/logo_variant.txt
# 		printf $"\n =========>	Sample $a Donor => Variant Calling: Haplotype Caller\n\n"

# 		java -Xmx64g -jar $GATK -T HaplotypeCaller \
# 		-R $REF \
# 		-I $PROCESSING/2_Add/$Data\_0$a\_G_Donor.bam \
# 		-o $PROCESSING/6_Variant/GATK/$Data\_0$a\_G_Donor.g.vcf \
# 		-ERC GVCF \
# 		--doNotRunPhysicalPhasing \
#  		--heterozygosity 0.001 \
#  		--indel_heterozygosity 1.25E-4 \
#  		--maxReadsInRegionPerSample 50000 \
# 		--min_base_quality_score 10 \
# 		--minReadsPerAlignmentStart 5 \
#  		--max_alternate_alleles 6 \
# 		-bamout $STORAGE/$Data\_0$a\_G_Donor.g.vcf.bam \
# 		-L $TARGETCELLFREE

# 		printf $"\n =========>	Sample $a Donor => Variant Calling: Done\n\n"
# 		cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt

# done

	

	Data=20170111

for (( a=1; a<5; a++ ))
do
	# cd $PROCESSING/6_Variant/GATK
	# ls 20170111_0$a\_*.g.vcf > sample_0$a.list #stampo il nome di tutti  files in una lista .list
	# echo $'\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling\n\n'

	# java -jar -Xmx64g $GATK -T GenotypeGVCFs \
	# -R $REF \
	# -V:VCF sample_0$a.list \
	# -o $PROCESSING/6_Variant/GATK/$Data\_0$a\_Donor_Rec_GATK.vcf

	# printf $"\n =========>	Genotype GVCF & Multi-sample variant calling: Done\n\n"
	# cat ~/Scrivania/SCRIPT_PIPELINE/logo_variant.txt
	# printf $"\n =========>	Sample $a => Variant Calling: FreeBayes\n\n"

	# /home/jarvis/NGS_TOOLS/freebayes/bin/freebayes -f $REF \
	# -b $PROCESSING/2_Add/$Data\_0$a\_G_Rec.bam \
	# -b $PROCESSING/2_Add/$Data\_0$a\_G_Donor.bam \
	# -K \
	# -J \
	# --pooled-continuous \
	# --pooled-discrete \
	# --genotype-qualities \
	# --report-genotype-likelihood-max \
	# --allele-balance-priors-off \
	# --min-mapping-quality 20 \
	# --min-base-quality 10 \
	# --min-alternate-fraction 0.1 \
	# --min-alternate-count 2 \
	# --min-coverage 10 \
	# --genotype-qualities \
	# --report-genotype-likelihood-max \
	# --allele-balance-priors-off \
	# -t $TARGETCELLFREEBED > $PROCESSING/6_Variant/FreeBayes/$Data\_0$a\_Donor_Rec_FreeBayes.vcf

	# printf $"\n =========>	Genotype GVCF & Multi-sample variant calling FreeBayes: Done\n\n"
	# cat ~/Scrivania/SCRIPT_PIPELINE/logo_variant.txt
	# printf $"\n =========>	Sample $a => Variant Calling: VarScan\n\n"
	# cd $PROCESSING/2_Add/

	# samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
	# -l $TARGETCELLFREEBED -b /home/jarvis/Scrivania/NGS_ANALYSIS_TEST/PROCESSING/Bam_list_0$a.list > $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec.mpileup

	# sed -e 's/.bam//g' /home/jarvis/Scrivania/NGS_ANALYSIS_TEST/PROCESSING/Bam_list_0$a.list > /home/jarvis/Scrivania/NGS_ANALYSIS_TEST/PROCESSING/Sample_list_0$a.list

	# java -jar -Xmx64g $VARSCAN mpileup2snp $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec.mpileup \
	# --min-coverage 10 \
	# --min-var-freq 0.20 \
	# --pvalue 0.05 \
	# --output-vcf 1 \
	# --vcf-sample-list /home/jarvis/Scrivania/NGS_ANALYSIS_TEST/PROCESSING/Sample_list_0$a.list > $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_snp.vcf

	# java -jar -Xmx64g $VARSCAN mpileup2indel $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec.mpileup \
	# --min-coverage 10 \
	# --min-var-freq 0.10 \
	# --pvalue 0.1 \
	# --output-vcf 1 \
	# --vcf-sample-list /home/jarvis/Scrivania/NGS_ANALYSIS_TEST/PROCESSING/Sample_list_0$a.list > $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_Indel.vcf


#done

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     NORM     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#  	echo "splitto GATK"
#  	python /home/jarvis/git/CMG/scripts/header_fix.py -v G -f $PROCESSING/6_Variant/GATK/$Data\_0$a\_Donor_Rec_GATK.vcf \
# 	> $PROCESSING/6_Variant/GATK/NORM_GATK/$Data\_0$a\_Donor_Rec_GATK.fix.vcf

# 	$BCFTOOLS norm -m -both \
# 	-f $REF \
# 	$PROCESSING/6_Variant/GATK/NORM_GATK/$Data\_0$a\_Donor_Rec_GATK.fix.vcf \
# 	> $PROCESSING/6_Variant/GATK/NORM_GATK/$Data\_0$a\_Donor_Rec_GATK.split.vcf


# #:::::::::::::::::::FREEBAYES:::::::::::::::::::

# 	echo "splitto FreeBayes"
# 	#python /home/jarvis/git/CMG/scripts/header_fix.py -v F -f $PROCESSING/6_Variant/FreeBayes/$Data\_0$a\_Donor_Rec_FreeBayes.vcf \
# 	#> $PROCESSING/6_Variant/FreeBayes/NORM_FREE/$Data\_0$a\_Donor_Rec_FreeBayes.fix.vcf

# 	$BCFTOOLS norm -m -both \
# 	-f $REF \
# 	$PROCESSING/6_Variant/FreeBayes/$Data\_0$a\_Donor_Rec_FreeBayes.vcf \
# 	> $PROCESSING/6_Variant/FreeBayes/NORM_FREE/$Data\_0$a\_Donor_Rec_FreeBayes.split.vcf



# #:::::::::::::::::::VARSCAN:::::::::::::::::::
# #Copio in Intersect
# echo "splitto VarScan"

# 	$BCFTOOLS norm -m -both \
# 	-f $REF \
# 	$PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_snp.vcf \
# 	> $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_snp.split.vcf

# 	$BCFTOOLS norm -m -both \
# 	-f $REF \
# 	$PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_Indel.vcf \
# 	> $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_Indel.split.vcf

# 	bgzip $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_snp.split.vcf
# 	bgzip $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_Indel.split.vcf
# 	tabix -f $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_snp.split.vcf.gz
# 	tabix -f $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_Indel.split.vcf.gz

# 	#Effettuo il merging dei file
# 	vcf-concat $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_snp.split.vcf.gz \
# 	$PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan_Indel.split.vcf.gz > $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan.merge.vcf
# 	vcf-sort -c $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan.merge.vcf > $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan.merge.sort.vcf

# 	echo "splitto VarScan"
	

# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::     FEATURES EXTRACTION     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 	
 	printf $"Estraggo le Features...\n"
 	python /home/jarvis/git/CMG/script_CMG/features_extraction_germline.py \
 	-g $PROCESSING/6_Variant/GATK/NORM_GATK/$Data\_0$a\_Donor_Rec_GATK.split.vcf \
 	-f $PROCESSING/6_Variant/FreeBayes/NORM_FREE/$Data\_0$a\_Donor_Rec_FreeBayes.split.vcf \
 	-v $PROCESSING/6_Variant/VarScan/$Data\_0$a\_Donor_Rec_VarScan.merge.sort.vcf \
 	-l $LISTAFEATURES \
 	-a -F \
 	-o $PROCESSING/6_Variant/Feature_extraction \
    -G $PROCESSING/6_Variant/GATK
	printf $"Features extraction----->DONE"
 	
# 	mv $PROCESSING/6_Variant/Feature_extraction/TOTAL.vcf $PROCESSING/6_Variant/Feature_extraction/$Data\_0$a\_TOTAL.vcf
#  	printf $"\n\n"
#  	cat ~/Scrivania/SCRIPT_PIPELINE/logo_annotation.txt
#  	printf $"\n\n\n"
#  	cd $VEP    
 
#  	perl $VEPANN -i $PROCESSING/6_Variant/Feature_extraction/$Data\_0$a\_TOTAL.vcf \
#  	-o $PROCESSING/8_Annotation/$Data\_0$a\_TOTAL.ANN.vcf \
#  	--stats_file $PROCESSING/8_Annotation/Feature_extraction/$Data\_0$a\_TOTAL.html \
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
# 	--pick
 	
 	
#  	mkdir -p $OUTVCF/$Data\_0$a\_Donor_Rec
 	
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
 	
 done