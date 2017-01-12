#!bin/bash

##TEMPO TOTALE DI PROCESSING: 2173 sec (PARALLELIZZATO)

SAMTOOLS=~/NGS_TOOLS/samtools-1.3.1/samtools
FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.12
PICARD=~/NGS_TOOLS/picard-tools-2.3.0/picard.jar
GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.4.0.jar
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
TARGPARALLFREE=/home/minime/NGS_ANALYSIS/TARGET/Parallel_regions_Cardio+-1000.regions


VEP=~/NGS_TOOLS/ensembl-tools-release-84/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-84/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-84/scripts/variant_effect_predictor/filter_vep.pl

DATE=$(date +"%Y%m%d");
STARTTIME=$(date +%s)

	/home/minime/NGS_TOOLS/freebayes/scripts/freebayes-parallel $TARGPARALLFREE 6 -f $REF \
	-L /home/minime/Scrivania/TEST/20161215/20161215_bamlist.list \
	--genotype-qualities \
	--report-genotype-likelihood-max \
	--allele-balance-priors-off \
	--min-mapping-quality 20 \
	--min-base-quality 10 \
	--min-alternate-fraction 0.1 \
	--min-alternate-count 2 \
	--min-coverage 10 >/home/minime/Scrivania/TEST/Cardio_FreeBayes.vcf
ENDTIME=$(date +%s)
echo "TEMPO TOTALE DI PROCESSING: $(($ENDTIME - $STARTTIME)) sec"

DATE=$(date +"%Y%m%d");
STARTTIME=$(date +%s)

	freebayes -f $REF \
	-L /home/minime/Scrivania/TEST/20161215/20161215_bamlist.list \
	--genotype-qualities \
	--report-genotype-likelihood-max \
	--allele-balance-priors-off \
	--min-mapping-quality 20 \
	--min-base-quality 10 \
	--min-alternate-fraction 0.1 \
	--min-alternate-count 2 \
	--min-coverage 10 \
	-t $TARGETCARDIOBED > /home/minime/Scrivania/TEST/Cardio_FreeBayes.2.vcf

ENDTIME=$(date +%s)
echo "TEMPO TOTALE DI PROCESSING: $(($ENDTIME - $STARTTIME)) sec"
printf "\n\n"
