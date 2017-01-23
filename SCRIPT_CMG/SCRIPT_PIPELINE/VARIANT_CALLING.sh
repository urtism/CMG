
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     GATK     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



HaplotypeCaller () {

	printf $"\n =========>	Sample $SAMPLE_NAME => Variant Calling: Haplotype Caller\n\n"

	java -Xmx64g -jar $GATK -T HaplotypeCaller \
	-R $REF \
	-I $1.bam \
	-o $WORKDIR/VARIANT_CALLING/GVCF/$SAMPLE_NAME.g.vcf \
	-ERC GVCF \
	--doNotRunPhysicalPhasing \
	-bamout $WORKDIR/VARIANT_CALLING/GVCF/$SAMPLE_NAME.g.vcf.bam \
	-L $TARGET

	printf $"\n =========>	Sample $SAMPLE_NAME => Variant Calling: Haplotype Caller: DONE\n\n"
	
}


GenotypeGVCFs () {
	
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_multi.txt

	echo $'\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling\n\n'

	java -jar -Xmx64g $GATK -T GenotypeGVCFs \
	-R $REF \
	-V:VCF $1 \
	-o $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.vcf

	python $SCRIPT_PIPELINE/header_fix.py -v G -f $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.vcf \
	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.fix.vcf

	$BCFTOOLS norm -m -both \
	-f $REF \
	$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.fix.vcf \
	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.split.vcf

	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.fix.vcf
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.vcf



	echo $'\n\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling: DONE\n'

{


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     FREEBAYES     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

FreeBayes () {

	echo $'\n =========>	Variant Calling with FreeBayes\n'

	$FREEBAYES -f $REF \
	-L $WORKDIR/Bam_list.txt \
	-K \
	-J \
	-s $1 \
	--genotype-qualities \
	--report-genotype-likelihood-max \
	--allele-balance-priors-off \
	-t $TARGET > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf


	python  $SCRIPT_PIPELINE/header_fix.py -v F -f $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf \
	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.fix.vcf

 	$BCFTOOLS norm -m -both \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.fix.vcf \
 	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.split.vcf

 	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.fix.vcf
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf

	echo $'\n =========>	Variant Calling with FreeBayes: DONE\n'
	
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARSCAN2     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

VarScan2_germline () {	

	echo $'\n =========>	Variant Calling with VarScan2\n\n'

	samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
	-l $TARGET -b $1 > $WORKDIR/PREPROCESSING/$DATA\_$PANNELLO.mpileup

	java -jar -Xmx64g $VARSCAN mpileup2snp $WORKDIR/PREPROCESSING/$DATA\_$PANNELLO.mpileup \
	--min-coverage 10 \
	--min-var-freq 0.20 \
	--pvalue 0.05 \
	--output-vcf 1 \
	--vcf-sample-list $2 > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf

	java -jar -Xmx64g $VARSCAN mpileup2indel $WORKDIR/PREPROCESSING/$DATA\_$PANNELLO.mpileup \
	--min-coverage 10 \
	--min-var-freq 0.10 \
	--pvalue 0.1 \
	--output-vcf 1 \
	--vcf-sample-list $2 > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf

	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf
	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf

	bgzip $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf
	bgzip $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf
	tabix $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf.gz
	tabix $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf.gz

	vcf-concat $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf.gz $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf.gz >$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.Merge.vcf
	vcf-sort -c $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.Merge.vcf > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.sort.vcf

	$BCFTOOLS norm -m -both \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.sort.vcf \
 	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.split.vcf


	rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.Merge.vcf
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.sort.vcf
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf.gz
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf.gz
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf.gz.tbi
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf.gz.tbi
	#rm $WORKDIR/PREPROCESSING/$DATA\_$PANNELLO.mpileup

	echo $'\n =========>	Variant Calling with VarScan2: DONE\n'

}

	
Features_extraction_germline () {


	printf $"Estraggo le Features...\n"

	if [ "$DESIGN" == "ENRICHMENT" ]
	then

		python $SCRIPT_PIPELINE/features_extraction_germline.py \
	 	-g $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.split.vcf \
	 	-f $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.split.vcf \
	 	-v $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.split.vcf \
	 	-l $LISTAFEATURES_GERMLINE \
	 	-F -s \
	 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION \
	    -G $WORKDIR/VARIANT_CALLING/GVCF

	elif [ "$DESIGN" == "AMPLICON" ]
	then

		python $SCRIPT_PIPELINE/features_extraction_germline.py \
	 	-g $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.split.vcf \
	 	-f $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.split.vcf \
	 	-v $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.split.vcf \
	 	-l $LISTAFEATURES_GERMLINE \
	 	-F -a -s \
	 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION \
	    -G $WORKDIR/VARIANT_CALLING/GVCF

	fi

	printf $"\n=========> Features extraction: DONE"

	mv $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/TOTAL.vcf $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.vcf

}

ANNOTATION_germline () {

	printf "\n\n"
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_annotation.txt
	printf $"\n\n\n"
	
	perl $VEPANN -i $1 \
	-o ${1%.*}.ANN.vcf \
	--stats_file ${1%.*}.ANN.html \
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

	echo $'\n =========>	ANNOTATION: DONE\n'

}


ADD_ANNOTATION_germline () {

	python $SCRIPT_PIPELINE/Estrai_Annotazione_Somatic.py \
  	-i $1 \
  	-f $2 \
  	-l $ANN_LIST_GERMLINE \
  	-t $TRANSCR_LIST \
  	-o ${2%.*}.ANN
 	
	
 	sort -V ${2%.*}.ANN.tsv > ${2%.*}.ANN.sort.tsv
  	sort -V ${2%.*}.ANN.Other_transcripts.tsv > ${2%.*}.ANN.Other_transcripts.sort.tsv
 	
  	rm ${2%.*}.ANN.tsv
  	rm ${2%.*}.ANN.Other_transcripts.tsv


}
	

VARIANT_CALLING_GERMLINE () {
	
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_variant.txt

	mkdir -p $WORKDIR/VARIANT_CALLING
	mkdir -p $WORKDIR/VARIANT_CALLING/GVCF
	mkdir -p $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION

	cat $FILEPATH | while read line
	do

		FASTQ1=$(echo "$line" | cut -f1)
		FASTQ2=$(echo "$line" | cut -f2)
		SAMPLE_NAME=$(echo "$line" | cut -f3)
				
		HaplotypeCaller $SAMPLE_NAME.bam

	done

	ls $WORKDIR/GVCF/*.g.vcf > $WORKDIR/gvcf.list

	GenotypeGVCFs $WORKDIR/gvcf.list

	ls $WORKDIR/PREPROCESSING/*.bam > $WORKDIR/Bam_list.txt
	
	cd $WORKDIR/PREPROCESSING
	ls *.bam > $WORKDIR/Sample_list.txt
	sed -i -e "s/.bam//g" $WORKDIR/Sample_list.txt
	cd 

	FreeBayes $WORKDIR/Bam_list.txt

	VarScan2_germline $WORKDIR/Bam_list.txt $WORKDIR/Sample_list.txt

	Features_extraction_germline

	ANNOTATION_germline $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.vcf

	
	cat $FILEPATH | while read line
	do

		FASTQ1=$(echo "$line" | cut -f1)
		FASTQ2=$(echo "$line" | cut -f2)
		SAMPLE_NAME=$(echo "$line" | cut -f3)
				
		ADD_ANNOTATION_germline $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.ANN.vcf $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME.features.tsv

	done





}