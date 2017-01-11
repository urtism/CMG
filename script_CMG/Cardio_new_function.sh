#!/bin/bash

cat ~/Scrivania/SCRIPT_PIPELINE/logo.txt 

FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.15
PICARD=~/NGS_TOOLS/picard-tools-2.7.1/picard.jar
GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_138.hg19.vcf
INPUT=~/Scrivania/NGS_ANALYSIS_TEST/INPUT_DATA/CARDIO
INPUTBRCA=~/Scrivania/NGS_ANALYSIS_TEST/INPUT_DATA/BRCA
INPUTCANCER=~/Scrivania/NGS_ANALYSIS_TEST/INPUT_DATA/CANCER
INPUTEXOME=~/Scrivania/NGS_ANALYSIS_TEST/INPUT_DATA/EXOME
PROCESSING=~/Scrivania/NGS_ANALYSIS_TEST/PROCESSING
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
OUTVCF=~/Scrivania/NGS_ANALYSIS_TEST/OUTPUT_DATA
STORAGE=~/Scrivania/NGS_ANALYSIS_TEST/STORAGE
VEP=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/filter_vep.pl


	for filename in $PROCESSING/5_BQSR/*.bam
	do

		cat ~/Scrivania/SCRIPT_PIPELINE/logo_variant.txt
		printf $"\n =========>	Sample ${filename%.*} => Variant Calling: Haplotype Caller\n\n"

		java -Xmx64g -jar $GATK -T HaplotypeCaller \
		-R $REF \
		-I $PROCESSING/5_BQSR/${filename%.*}.bam \
		-o $PROCESSING/6_Variant/GATK/${filename%.*}.g.vcf \
		-ERC GVCF \
		--doNotRunPhysicalPhasing \
		--doNotRunPhysicalPhasing \
		--heterozygosity 0.001 \
		--indel_heterozygosity 1.25E-4 \
		--maxReadsInRegionPerSample 50000 \
		--min_base_quality_score 10 \
		--minReadsPerAlignmentStart 5 \
		--max_alternate_alleles 6 \
		#-bamout $PROCESSING/6_Variant/GATK/${filename%.*}.g.vcf.bam \
		-L $TARGET

		printf $"\n =========>	Sample ${filename%.*} => Variant Calling: Haplotype Caller: DONE\n\n"


	done




#::::::::::::::::::::::::::::::::::::::::::     INIZIO VARIANT CALLING MULTI-SAMPLE     :::::::::::::::::::::::::::::::::::::::::::::::





	#Stampo tutti i file .g.vcf in una lista
	ls $PROCESSING/6_Variant/GATK/*.g.vcf > samples.list #stampo il nome di tutti  files in una lista .list
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_multi.txt
	echo $'\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling\n\n'

		java -jar -Xmx64g $GATK -T GenotypeGVCFs \
		-R $REF \
		-V:VCF samples.list \
		-A HomopolymerRun \
#		-newQual \
		-o $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK.vcf

	# Salvo nella cartella storage
	#cp $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK.vcf $STORAGE/$Name_Dir

	echo $'\n\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling: DONE\n'




#::::::::::::::::::::::::::::::::::::::::::::::::::::::::     FILTRAGGIO GATK     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



	
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     FREEBAYES     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



	ls $PROCESSING/5_BQSR/*.bam > $PROCESSING/5_BQSR/Bam_list.txt
	ls $PROCESSING/5_BQSR/*.bam > $PROCESSING/5_BQSR/Sample_list.txt
	sed -i -e "s/.bam//g" $PROCESSING/5_BQSR/Sample_list.txt

	/home/jarvis/NGS_TOOLS/freebayes/bin/freebayes -f $REF \
	-L $PROCESSING/5_BQSR/Bam_list.txt \
	-K \
	-J \
	-s $PROCESSING/5_BQSR/Sample_list.txt \
	--pooled-continuous \
	--pooled-discrete \
	--genotype-qualities \
	--report-genotype-likelihood-max \
	--allele-balance-priors-off \
	--min-mapping-quality 20 \
	--min-base-quality 10 \
	--min-alternate-fraction 0.1 \
	--min-alternate-count 2 \
	--min-coverage 10 \
	-t $TARGETCARDIOBED > $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes.vcf

	echo $'\n =========>	Variant Calling with FreeBayes: DONE\n'
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	FreeBayes vcf manipulation for intersection\n\n'

	#Rimuovo le righe duplicate:
	#awk '!a[$0]++' $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes.vcf > $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes_RD.vcf
	
	#Converto il vcf alla versione 4.2
	#sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes_RD.vcf

	java -Xmx64g -jar $PICARD SortVcf \
	I=$PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes_RD.vcf \
	O=$PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes_Sort.vcf

	cp $DataRun\_Cardio_FreeBayes_Sort.vcf $STORAGE/$Name_Dir/$DataRun\_Cardio_FreeBayes.vcf

	#Copio il vcf nella cartella $PROCESSING/6_Variant/Intersect
	cp $DataRun\_Cardio_FreeBayes_Sort.vcf $PROCESSING/6_Variant/Intersect/
	bgzip $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_FreeBayes_Sort.vcf
	tabix $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_FreeBayes_Sort.vcf.gz
	rm $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes_RD.vcf
	rm $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes_Sort.vcf
	rm $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes.vcf
	rm $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes_Sort.vcf.idx
	rm $PROCESSING/5_BQSR/Sample_list.txt

	echo $'\n =========>	FreeBayes vcf manipulation for intersection: DONE\n'
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	Variant Calling with VarScan2\n\n'



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARSCAN2     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



	samtools mpileup -B -q 1 -f $REF \
	-l $TARGETCARDIOBED -b $PROCESSING/5_BQSR/Bam_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup

	cp $PROCESSING/5_BQSR/Bam_list.txt $PROCESSING/5_BQSR/Sample_list.txt
	sed -i -e 's/.bam//g' $PROCESSING/5_BQSR/Sample_list.txt

	java -jar -Xmx64g $VARSCAN mpileup2snp $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup \
	--min-coverage 10 \
	--min-var-freq 0.20 \
	--pvalue 0.05 \
	--output-vcf 1 \
	--vcf-sample-list $PROCESSING/5_BQSRSample_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf

	java -jar -Xmx64g $VARSCAN mpileup2indel $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup \
	--min-coverage 10 \
	--min-var-freq 0.10 \
	--pvalue 0.1 \
	--output-vcf 1 \
	--vcf-sample-list $PROCESSING/5_BQSRSample_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf

	echo $'\n =========>	Variant Calling with VarScan2: DONE\n'
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	VarScan2 vcf manipulation for intersection\n\n'
	
	#Modifico la versione del vcf
	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf
	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf

	#Copio in Intersect
	cp $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf $PROCESSING/6_Variant/Intersect/
	cp $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf $PROCESSING/6_Variant/Intersect/
	bgzip $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_snp.vcf
	bgzip $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Indel.vcf
	tabix $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_snp.vcf.gz
	tabix $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Indel.vcf.gz

	#Effettuo il merging dei file
	cd $PROCESSING/6_Variant/Intersect/
	
	vcf-concat $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_snp.vcf.gz \
	$PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Indel.vcf.gz > $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge.vcf
	
	vcf-sort -c $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge.vcf > $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge_Sort.vcf

	cp $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge_Sort.vcf $STORAGE/$Name_Dir/$DataRun\_Cardio_VarScan.vcf

	bgzip $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge_Sort.vcf
	tabix $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge_Sort.vcf.gz

	#Rimuovo i file inutili:
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_snp.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Indel.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_snp.vcf.gz.tbi
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Indel.vcf.gz.tbi
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge.vcf
	rm $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf
	rm $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf
	rm $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup
	
	echo $'\n =========>	VarScan2 vcf manipulation for intersection: DONE\n'
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	GATK + FreeBayes + VarScan2 INTERSECTION\n\n'



#::::::::::::::::::::::::::::::::::::::::::::::::::::::::    GATK ANNOTATION     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



	printf "\n\n"
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_annotation.txt
	printf $"\n\n\n"

	cd $VEP

	perl $VEPANN -i $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK_Filter.vcf \
	-o $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_ANN.vcf \
	--stats_file $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_ANN.html \
	--cache \
	--assembly GRCh37 \
	--offline \
	--force_overwrite \
	-v \
	--fork 10 \
	--variant_class \
	--sift b \
	--poly b \
	--vcf_info_field ANN \
	--hgvs \
	--protein \
	--canonical \
	--check_existing \
	--gmaf \
	--pubmed \
	--species homo_sapiens \
	--failed 1 \
	--vcf

	echo $'\n =========>	GATK ANNOTATION: DONE\n'
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	INTERSECTION ANNOTATION\n\n'

	rm $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK_Filter.vcf
	rm $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK_Filter.vcf.idx
	cp $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_ANN.vcf $STORAGE/$Name_Dir
	printf "\n${DataRun}_Cardio_GATK_ANN.vcf : contiene il file ${DataRun}_Cardio_GATK_Filter.vcf annotato\n" >> $STORAGE/$Name_Dir/README.txt
	mv $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_ANN.html $STORAGE/$Name_Dir
	printf "\n${DataRun}_Cardio_GATK_ANN.html : contiene le statistiche del file ${DataRun}_Cardio_GATK_Filter.vcf\n" >> $STORAGE/$Name_Dir/README.txt



#:::::::::::::::::::::::::::::::::::::::::::::::::::    INTERSECTION  ANNOTATION     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::



	cd $VEP

	perl $VEPANN -i $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf \
	-o $PROCESSING/8_Annotation/$DataRun\_Cardio_Total_Intersect_ANN.vcf \
	--stats_file $PROCESSING/8_Annotation/$DataRun\_Cardio_Total_Intersect_ANN.html \
	--cache \
	--assembly GRCh37 \
	--offline \
	--force_overwrite \
	-v \
	--fork 10 \
	--variant_class \
	--sift b \
	--poly b \
	--vcf_info_field ANN \
	--hgvs \
	--protein \
	--canonical \
	--check_existing \
	--gmaf \
	--pubmed \
	--species homo_sapiens \
	--failed 1 \
	--vcf

	rm $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf
	cp $PROCESSING/8_Annotation/$DataRun\_Cardio_Total_Intersect_ANN.vcf $STORAGE/$Name_Dir/$DataRun\_Cardio_Intersect_ANN.vcf
	printf "\n${DataRun}_Cardio_Intersect_ANN.vcf : contiene il file ${DataRun}_Cardio_Intersect.vcf annotato\n" >> $STORAGE/$Name_Dir/README.txt
	mv $PROCESSING/8_Annotation/$DataRun\_Cardio_Total_Intersect_ANN.html $STORAGE/$Name_Dir
	printf "\n${DataRun}_Cardio_Intersect_ANN.html : contiene le statistiche del file ${DataRun}_Cardio_Intersect_ANN.vcf\n" >> \
	$STORAGE/$Name_Dir/README.txt

	echo $'\n =========>	INTERSECTION ANNOTATION: DONE\n'
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	FILE MODIFICATION FOR CSV FORMAT \n\n'

