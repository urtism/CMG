

FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.12
PICARD=~/NGS_TOOLS/picard-tools-2.3.0/picard.jar
GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_135.hg19.vcf

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

VEP=~/NGS_TOOLS/ensembl-tools-release-84/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-84/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-84/scripts/variant_effect_predictor/filter_vep.pl


INPUT=/home/jarvis/TEST_TOOLS/INPUT_TEST
PROCESSING=/home/jarvis/TEST_TOOLS/PROCESSING_TEST
OUTPUT=/home/jarvis/TEST_TOOLS/OUTPUT_TEST

DataRun=20160930

#GATK

#	cd $INPUT
#	for filename in *.bam
#	do

#		java -Xmx64g -jar $GATK -T HaplotypeCaller \
#		-R $REF \
#		-I $INPUT/${filename%.*}.bam \
#		-o $PROCESSING/${filename%.*}.g.vcf \
#		-ERC GVCF \
#		--doNotRunPhysicalPhasing \
#		--heterozygosity 0.001 \
#		--indel_heterozygosity 1.25E-4 \
#		--maxReadsInRegionPerSample 50000 \
#		--min_base_quality_score 10 \
#		--minReadsPerAlignmentStart 5 \
#		--max_alternate_alleles 6 \
#		-L $TARGET

#	done

#MULTI-SAMPLE

#	cd $PROCESSING
	#Stampo tutti i file .g.vcf in una lista
#	ls *.g.vcf > samples.list #stampo il nome di tutti  files in una lista .list

#		java -jar -Xmx64g $GATK -T GenotypeGVCFs \
#		-R $REF \
#		-V:VCF samples.list \
#		-o $OUTPUT/$DataRun\_Cardio_GATK.vcf

#FREEBAYES

#	cd $INPUT
#	ls *.bam > Bam_list.txt

#	freebayes -f $REF \
#	-L Bam_list.txt \
#	--pooled-continuous \
#	--pooled-discrete \
#	--genotype-qualities \
#	--report-genotype-likelihood-max \
#	--allele-balance-priors-off \
#	--min-mapping-quality 20 \
#	--min-base-quality 10 \
#	--min-alternate-fraction 0.1 \
#	--min-alternate-count 2 \
#	--min-coverage 10 \
#	-t $TARGETCARDIOBED > $OUTPUT/$DataRun\_Cardio_FreeBayes.vcf

	#Rimuovo le righe duplicate:
#	cd $OUTPUT
#	awk '!a[$0]++' $DataRun\_Cardio_FreeBayes.vcf > $DataRun\_Cardio_FreeBayes_RD.vcf
	
	#Converto il vcf alla versione 4.2
#	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $DataRun\_Cardio_FreeBayes_RD.vcf

#	java -Xmx64g -jar $PICARD SortVcf \
#	I=$DataRun\_Cardio_FreeBayes_RD.vcf \
#	O=$DataRun\_Cardio_FreeBayes_Sort.vcf

#VARSCAN

	cd $INPUT
	samtools mpileup \
	-f $REF \
	-l $TARGETCARDIOBED \
	-b $INPUT/Bam_list.txt \
	-B \
	-q 20 \
	-Q 10 \
	> $OUTPUT/$DataRun\_Cardio.mpileup

	cp Bam_list.txt Sample_list.txt
	sed -i -e 's/.bam//g' Sample_list.txt

	java -jar -Xmx64g $VARSCAN mpileup2snp $OUTPUT/$DataRun\_Cardio.mpileup \
	--min-coverage 10 \
	--min-reads2 2 \
	--min-avg-qual 10 \
	--min-var-freq 0.10 \
	--min-freq-for-hom 0.75 \
	--pvalue 0.99 \
	--output-vcf 1 \
	--vcf-sample-list Sample_list.txt > $OUTPUT/$DataRun\_Cardio_VarScan_snp.vcf

	java -jar -Xmx64g $VARSCAN mpileup2indel $OUTPUT/$DataRun\_Cardio.mpileup \
	--min-coverage 10 \
	--min-reads2 2 \
	--min-avg-qual 10 \
	--min-var-freq 0.10 \
	--min-freq-for-hom 0.75 \
	--pvalue 0.99 \
	--output-vcf 1 \
	--vcf-sample-list $INPUT/Sample_list.txt > $OUTPUT/$DataRun\_Cardio_VarScan_Indel.vcf
	
#	cd $PROCESSING/6_Variant/VarScan/
	#Modifico la versione del vcf
#	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $DataRun\_Cardio_VarScan_Indel.vcf
#	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $DataRun\_Cardio_VarScan_snp.vcf


#	#Effettuo il merging dei file
#	cd $PROCESSING/6_Variant/Intersect/
#	vcf-concat $DataRun\_Cardio_VarScan_snp.vcf.gz $DataRun\_Cardio_VarScan_Indel.vcf.gz > $DataRun\_Cardio_VarScan_Merge.vcf
#	vcf-sort -c $DataRun\_Cardio_VarScan_Merge.vcf > $DataRun\_Cardio_VarScan_Merge_Sort.vcf
