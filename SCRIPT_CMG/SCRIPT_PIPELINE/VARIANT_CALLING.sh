save_in_storage () {
	cat $1 | while read line
	do
		mv $line $STORAGE
	done
}




#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     GATK     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Mutect2 () {
	printf $"\n =========>	Sample $3 => Variant Calling: MuTect2\n\n"
	
	java -Xmx64g -jar $GATK -T MuTect2 \
	-R $REF \
	-I:tumor $1 \
	-I:normal $2 \
	-L $TARGET \
	-o $WORKDIR/VARIANT_CALLING/$3\_Sane_GATK.vcf \
	-bamout $WORKDIR/VARIANT_CALLING/$3\_Sane_GATK.vcf.bam

	$BCFTOOLS norm -m -both \
	-f $REF \
	$WORKDIR/VARIANT_CALLING/$3\_SANE_GATK.vcf \
	> $WORKDIR/VARIANT_CALLING/$3\_Sane_GATK.split.vcf 

	VCF_MUTECT=$WORKDIR/VARIANT_CALLING/$3\_Sane_GATK.split.vcf

	printf $"\n =========>	Sample $3 => Variant Calling: MuTect2: DONE\n\n"
}


HaplotypeCaller () {

	printf $"\n =========>	Sample $SAMPLE_NAME => Variant Calling: Haplotype Caller\n\n"

	java -Xmx64g -jar $GATK -T HaplotypeCaller \
	-R $REF \
	-I $1 \
	-o $WORKDIR/VARIANT_CALLING/GVCF/$SAMPLE_NAME.g.vcf \
	-ERC GVCF \
	--doNotRunPhysicalPhasing \
	-bamout $WORKDIR/VARIANT_CALLING/GVCF/$SAMPLE_NAME.g.vcf.bam \
	-L $TARGET

	printf $"$WORKDIR/VARIANT_CALLING/GVCF/$SAMPLE_NAME.g.vcf\n" >> $WORKDIR/gvcf.list
	printf $"\n =========>	Sample $SAMPLE_NAME => Variant Calling: Haplotype Caller: DONE\n\n"
	
}


GenotypeGVCFs () {
	
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_multi.txt

	printf $'\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling\n\n'

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

	VCF_GATK=$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.split.vcf

	printf $'\n\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling: DONE\n'

}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     FREEBAYES     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

FreeBayes () {

	printf $'\n =========>	Variant Calling with FreeBayes\n'

	$FREEBAYES -f $REF \
	-L $1 \
	-K \
	-J \
	-s $2 \
	--genotype-qualities \
	--report-genotype-likelihood-max \
	--allele-balance-priors-off \
	-t $TARGETBED > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf

	python  $SCRIPT_PIPELINE/header_fix.py -v F -f $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf \
	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.fix.vcf

 	$BCFTOOLS norm -m -both \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.fix.vcf \
 	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.split.vcf

 	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.fix.vcf
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf

	VCF_FREEBAYES=$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.split.vcf

	printf $'\n =========>	Variant Calling with FreeBayes: DONE\n\n'
	
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARSCAN2     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

VarScan2_germline () {	

	printf $'\n =========>	Variant Calling with VarScan2\n\n'

	samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
	-l $TARGETBED -b $1 > $WORKDIR/PREPROCESSING/$DATA\_$PANNELLO.mpileup

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

 	VCF_VARSCAN=$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.split.vcf

	rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.Merge.vcf
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.sort.vcf
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf.gz
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf.gz
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf.gz.tbi
	#rm $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf.gz.tbi
	#rm $WORKDIR/PREPROCESSING/$DATA\_$PANNELLO.mpileup

	printf $'\n =========>	Variant Calling with VarScan2: DONE\n'

}


VarScan2_somatic () {
	
	printf $"\n =========>	Sample $3 => Variant Calling: VarScan2\n\n"
	
	samtools mpileup -f $REF -l $TARGETBED -q 1 -B $2 $1 > $WORKDIR/VARIANT_CALLING/$3\_Sane.mpileup 
 		
	java -jar -Xmx64g $VARSCAN somatic $WORKDIR/VARIANT_CALLING/$3\_Sane.mpileup \
	$WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan --output-vcf 1 --mpileup 1

	VCF_VARSCAN_SNP=$WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.snp.vcf
	VCF_VARSCAN_INDEL=$WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.indel.vcf
	
	printf $"\n =========>	Sample $3 => Variant Calling: Varscan2: DONE\n\n"	
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARDICT     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

VarDict () {	

	printf $"\n =========>	Sample $3 => Variant Calling: VarDict\n\n"
	
	$VARDICT -G $REF -f 0.01 -N $3 -b "$1|$2" \
	-z -F 0 -c 1 -S 2 -E 3 -g 4 $TARGETBED | ~/NGS_TOOLS/VarDictJava-master/VarDict/testsomatic.R | ~/NGS_TOOLS/VarDictJava-master/VarDict/var2vcf_paired.pl \
	-N "$3|$4" -f 0.01 > $WORKDIR/VARIANT_CALLING/$3\_Sane_VarDict.vcf

	VCF_VARDICT=$WORKDIR/VARIANT_CALLING/$3\_SANE_VarDict.vcf
	
	printf $"\n =========>	Sample $3 => Variant Calling: VarDict: DONE\n\n"	
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::     FEATURES EXTRACTION     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Features_extraction_germline () {

	rm -f $WORKDIR/PostFeatureExtraction.cfg
	CFG=$WORKDIR/PostFeatureExtraction.cfg
	
	cat $1 | while read line
	do
		VCF_GATK=$(echo "$line" | cut -f1)
		VCF_FREEBAYES=$(echo "$line" | cut -f2)
		VCF_VARSCAN=$(echo "$line" | cut -f3)
		
		printf $"Estraggo le Features...\n"

		if [ "$DESIGN" == "ENRICHMENT" ]
		then

			python $SCRIPT_PIPELINE/features_extraction_germline.py \
		 	-g $VCF_GATK \
		 	-f $VCF_FREEBAYES \
		 	-v $VCF_VARSCAN \
		 	-l $LISTAFEATURES_GERMLINE \
		 	-F -s \
		 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION \
		    -G $WORKDIR/VARIANT_CALLING/GVCF

		elif [ "$DESIGN" == "AMPLICON" ]
		then

			python $SCRIPT_PIPELINE/features_extraction_germline.py \
		 	-g $VCF_GATK \
		 	-f $VCF_FREEBAYES \
		 	-v $VCF_VARSCAN \
		 	-l $LISTAFEATURES_GERMLINE \
		 	-F -a -s \
		 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION \
		    -G $WORKDIR/VARIANT_CALLING/GVCF

		fi

		mv $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/TOTAL.vcf $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.vcf
		cp $VCF_GATK $OUT
		cp $VCF_FREEBAYES $OUT
		cp $VCF_VARSCAN $OUT
		mv $VCF_GATK $STORAGE
		mv $VCF_FREEBAYES $STORAGE
		mv $VCF_VARSCAN $STORAGE
		save_in_storage $WORKDIR/gvcf.list

		INPUT=$WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.vcf
		#printf $"$WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.vcf\n" >> $CFG
		printf $"\n=========> Features extraction: DONE"
	done

}


Features_extraction_somatic () {

	rm -f $WORKDIR/PostFeatureExtraction.cfg
	CFG=$WORKDIR/PostFeatureExtraction.cfg

	cat $1 | while read line
	do
		VCF_MUTECT=$(echo "$line" | cut -f1)
		VCF_VARDICT=$(echo "$line" | cut -f2)
		VCF_VARSCAN_SNP=$(echo "$line" | cut -f3)
		VCF_VARSCAN_INDEL=$(echo "$line" | cut -f4)
		SAMPLE_NAME_SOM=$(echo "$line" | cut -f5)
		SAMPLE_NAME_NORM=$(echo "$line" | cut -f6)

		printf $"\n =========>	Sample $SAMPLE_NAME_SOM => Features extraction\n\n"

		if [ "$DESIGN" == "ENRICHMENT" ]
		then

		 	python $SCRIPT_PIPELINE/features_extraction_somatic.py \
		 	-m $VCF_MUTECT \
		 	-d $VCF_VARDICT \
		 	-v $VCF_VARSCAN_SNP \
		 	-i $VCF_VARSCAN_INDEL \
		 	-n $SAMPLE_NAME_NORM -t $SAMPLE_NAME_SOM \
		 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME_SOM
			
		elif [ "$DESIGN" == "AMPLICON" ]
		then

			python $SCRIPT_PIPELINE/features_extraction_somatic.py \
		 	-m $VCF_MUTECT \
		 	-d $VCF_VARDICT \
		 	-v $VCF_VARSCAN_SNP \
		 	-i $VCF_VARSCAN_INDEL \
		 	-n $SAMPLE_NAME_NORM -t $SAMPLE_NAME_SOM -a \
		 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME_SOM
		fi
		
		cp $VCF_MUTECT $OUT
		cp $VCF_VARDICT $OUT
		cp $VCF_VARSCAN_SNP $OUT
		cp $VCF_VARSCAN_INDEL $OUT


		mv $VCF_MUTECT $STORAGE
		mv $VCF_VARDICT $STORAGE
		mv $VCF_VARSCAN_SNP $STORAGE
		mv $VCF_VARSCAN_INDEL $STORAGE

		printf $"\n =========>	Sample $SAMPLE_NAME_SOM => Features extraction: DONE\n\n"
		printf $"$WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME_SOM\n" >> $CFG
	done
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::     ANNOTATION     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ANNOTATION_germline () {

	printf "\n\n"
	cat $LOGHI/logo_annotation.txt
	printf $"\n\n\n"
	
	perl $VEPANN -i $1 \
	-o ${1%.*}.ANN.vcf \
	--stats_file ${1%.*}.ANN.html \
	--cache \
	--assembly GRCh37 \
	--offline \
	--force_overwrite \
	-v \
	--fork 2 \
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

	mv $1 $DELETE
	INPUT=${1%.*}.ANN.vcf
	printf $'\n =========>	ANNOTATION: DONE\n'

}


ADD_ANNOTATION_germline () {

	python $SCRIPT_PIPELINE/Estrai_Annotazione_Somatic.py \
  	-i $1 \
  	-f ${1%.*.*}.features.tsv \
  	-l $ANN_LIST_GERMLINE \
  	-t $TRANSCR_LIST \
  	-o  ${1%.*.*}.ANN
 	
 	sort -V ${1%.*.*}.ANN.tsv > ${1%.*.*}.ANN.sort.tsv
  	sort -V ${1%.*.*}.ANN.Other_transcripts.tsv > ${1%.*.*}.ANN.Other_transcripts.sort.tsv
 	
 	cp ${1%.*.*}.ANN.sort.tsv $STORAGE
 	cp ${1%.*.*}.ANN.Other_transcripts.sort.tsv $STORAGE
 	mv ${1%.*.*}.ANN.sort.tsv $OUT
 	mv ${1%.*.*}.ANN.Other_transcripts.sort.tsv $OUT

  	mv ${2%.*}.ANN.tsv $DELETE
  	mv ${2%.*}.ANN.Other_transcripts.tsv $DELETE
}


ANNOTATION_somatic () {

	printf "\n\n"
	cat $LOGHI/logo_annotation.txt
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

 	mv $1 $DELETE
 	INPUT=${1%.*}.ANN.vcf
 	printf $'\n =========>	ANNOTATION: DONE\n'
}


ADD_ANNOTATION_somatic () {
	
	python $SCRIPT_PIPELINE/Estrai_Annotazione_Somatic.py \
  	-i $1 \
  	-f ${1%.*.*}.features.tsv \
  	-l $ANN_LIST_SOMATIC \
 	-t $TRANSCR_LIST \
  	-o ${1%.*.*}.ANN
 	
 	sort -V ${1%.*.*}.ANN.tsv > ${1%.*.*}.ANN.sort.tsv
  	sort -V ${1%.*.*}.ANN.Other_transcripts.tsv > ${1%.*.*}.ANN.Other_transcripts.sort.tsv
 	
 	cp ${1%.*.*}.ANN.sort.tsv $STORAGE
 	cp ${1%.*.*}.ANN.Other_transcripts.sort.tsv $STORAGE
 	mv ${1%.*.*}.ANN.sort.tsv $OUT
 	mv ${1%.*.*}.ANN.Other_transcripts.sort.tsv $OUT
 	
  	mv ${2%.*}.ANN.tsv $DELETE
  	mv ${2%.*}.ANN.Other_transcripts.tsv $DELETE
}


VARIANT_CALLING_GERMLINE () {
	
	cat $LOGHI/logo_variant.txt

	mkdir -p $WORKDIR/VARIANT_CALLING
	mkdir -p $WORKDIR/VARIANT_CALLING/GVCF
	mkdir -p $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION
	
	rm -f $WORKDIR/Bam_list.txt
	rm -f $WORKDIR/Sample_list.txt
	rm -f $WORKDIR/gvcf.list
	rm -f $WORKDIR/PostVariantCalling.cfg
	CFG=$WORKDIR/PostVariantCalling.cfg

	cat $1 | while read line
	do

		BAM=$(echo "$line" | cut -f1)
		SAMPLE_NAME=$(echo "$line" | cut -f2)
			
		HaplotypeCaller $BAM

		printf $"$BAM\n" >> $WORKDIR/Bam_list.txt
		printf $"$SAMPLE_NAME\n" >> $WORKDIR/Sample_list.txt

	done

	GenotypeGVCFs $WORKDIR/gvcf.list
	
	FreeBayes $WORKDIR/Bam_list.txt $WORKDIR/Sample_list.txt

	VarScan2_germline $WORKDIR/Bam_list.txt $WORKDIR/Sample_list.txt

	save_in_storage $WORKDIR/Bam_list.txt

	printf $"$VCF_GATK\t$VCF_FREEBAYES\t$VCF_VARSCAN\n" >> $CFG
	
	# cat $1 | while read line
	# do

	# 	FASTQ1=$(echo "$line" | cut -f1)
	# 	FASTQ2=$(echo "$line" | cut -f2)
	# 	SAMPLE_NAME=$(echo "$line" | cut -f3)
				
	# 	ADD_ANNOTATION_germline $INPUT $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME.features.tsv;

	# done
}



VARIANT_CALLING_SOMATIC () {
	
	cat $LOGHI/logo_variant.txt

	mkdir -p $WORKDIR/VARIANT_CALLING
	mkdir -p $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION
	
	rm -f $WORKDIR/PostVariantCalling.cfg
	CFG=$WORKDIR/PostVariantCalling.cfg

	cat $1 | while read line
	do

		BAM_SOM=$(echo "$line" | cut -f1)
		SAMPLE_NAME_SOM=$(echo "$line" | cut -f2)
		BAM_NORM=$(echo "$line" | cut -f3)
		SAMPLE_NAME_NORM=$(echo "$line" | cut -f4)
			
		Mutect2 $BAM_SOM $BAM_NORM $SAMPLE_NAME_SOM $SAMPLE_NAME_NORM

		VarScan2_somatic $BAM_SOM $BAM_NORM $SAMPLE_NAME_SOM $SAMPLE_NAME_NORM
		
		VarDict $BAM_SOM $BAM_NORM $SAMPLE_NAME_SOM $SAMPLE_NAME_NORM

		mv $BAM_SOM $STORAGE
		mv $BAM_NORM $STORAGE

		printf $"$VCF_MUTECT\t$VCF_VARDICT\t$VCF_VARSCAN_SNP\t$VCF_VARSCAN_INDEL\t$SAMPLE_NAME_SOM\t$SAMPLE_NAME_NORM\n" >> $CFG
	done

}