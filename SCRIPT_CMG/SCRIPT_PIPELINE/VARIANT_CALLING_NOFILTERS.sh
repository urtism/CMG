save_in_storage () {

	cat $1 | while read line
	do
		mv $line $STORAGE
		mv $line.idx $STORAGE
		mv ${line%.*}.bai $STORAGE
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
	--min_base_quality_score 30 \
	-o $WORKDIR/VARIANT_CALLING/$3\_Sane_GATK.vcf 

	$BCFTOOLS norm -m -both \
	-f $REF \
	$WORKDIR/VARIANT_CALLING/$3\_Sane_GATK.vcf \
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
	--min_base_quality_score 1 \
	--minReadsPerAlignmentStart 1 \
	--standard_min_confidence_threshold_for_calling 5 \
	--max_alternate_alleles 20 \
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

	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.fix.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.vcf.idx $DELETE

	VCF_GATK=$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_GATK.split.vcf

	printf $'\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling: DONE\n'

}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     FREEBAYES     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

FreeBayes () {

	printf $'\n =========>	Variant Calling with FreeBayes\n\n'

	$FREEBAYES -f $REF \
	-L $1 \
	-K \
	-J \
	-s $2 \
	--genotype-qualities \
	--report-genotype-likelihood-max \
	--allele-balance-priors-off \
	--min-mapping-quality 1 \
	--min-supporting-allele-qsum 0 \
	--min-supporting-mapping-qsum 0 \
	--min-alternate-fraction 0.01 \
	--min-alternate-count 1 \
	--min-alternate-qsum 0 \
	--min-base-quality 0 \
	-t $TARGETBED > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf

	python  $SCRIPT_PIPELINE/header_fix.py -v F -f $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf \
	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.fix.vcf

 	$BCFTOOLS norm -m -both \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.fix.vcf \
 	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.split.vcf

 	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.fix.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf $DELETE

	VCF_FREEBAYES=$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.split.vcf

	printf $'\n =========>	Variant Calling with FreeBayes: DONE\n\n'
	
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARSCAN2     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

VarScan2_germline () {	

	printf $'\n =========>	Variant Calling with VarScan2\n\n'

	samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
	-l $TARGETBED -b $1 > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO.mpileup

	java -jar -Xmx64g $VARSCAN mpileup2snp $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO.mpileup \
	--min-coverage 1 \
	--min-reads2 1 \
	--min-avg-qual 1 \
	--min-var-freq 0.01 \
	--pvalue 0.05 \
	--strand-filter 0 \
	--output-vcf 1 \
	--vcf-sample-list $2 > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf

	java -jar -Xmx64g $VARSCAN mpileup2indel $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO.mpileup \
	--min-coverage 1 \
	--min-reads2 1 \
	--min-avg-qual 1 \
	--min-var-freq 0.01 \
	--pvalue 0.05 \
	--strand-filter 0 \
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

	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.Merge.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan.sort.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf.gz.tbi $DELETE
	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf.gz.tbi $DELETE
	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO.mpileup $DELETE

	printf $'\n =========>	Variant Calling with VarScan2: DONE\n'

}


VarScan2_somatic () {
	
	printf $"\n =========>	Sample $3 => Variant Calling: VarScan2\n\n"
	
	samtools mpileup -f $REF -l $TARGETBED -Q 30 -d 50000 -L 50000 -q 1 -B $2 $1 > $WORKDIR/VARIANT_CALLING/$3\_Sane.mpileup 
 		
	java -jar -Xmx64g $VARSCAN somatic $WORKDIR/VARIANT_CALLING/$3\_Sane.mpileup \
	$WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan \
	--min-var-freq 0.0005 \
	--output-vcf 1 \
	--mpileup 1

	$BCFTOOLS norm -m -both \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.snp.vcf \
 	> $WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.norm.snp.vcf

 	python $SCRIPT_PIPELINE/header_fix.py -v V -f $WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.indel.vcf \
	> $WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.fix.indel.vcf

 	$BCFTOOLS norm -m -both \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.fix.indel.vcf \
 	> $WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.norm.indel.vcf

	VCF_VARSCAN_SNP=$WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.norm.snp.vcf
	VCF_VARSCAN_INDEL=$WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.norm.indel.vcf
	
	mv $WORKDIR/VARIANT_CALLING/$3\_Sane.mpileup $DELETE
	printf $"\n =========>	Sample $3 => Variant Calling: Varscan2: DONE\n\n"	
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARDICT     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

VarDict () {	

	printf $"\n =========>	Sample $3 => Variant Calling: VarDict\n\n"
	
	$VARDICT -G $REF -f 0.0005 -N $3 -b "$1|$2" \
	-z -F 0 -c 1 -S 2 -E 3 -q 30 -g 4 $TARGETBED | ~/NGS_TOOLS/VarDictJava-master/VarDict/testsomatic.R | ~/NGS_TOOLS/VarDictJava-master/VarDict/var2vcf_paired.pl \
	-N "$3|$4" -f 0.0005 > $WORKDIR/VARIANT_CALLING/$3\_Sane_VarDict.vcf

	VCF_VARDICT=$WORKDIR/VARIANT_CALLING/$3\_Sane_VarDict.vcf
	
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
		
		printf $"\n=========> Features extraction\n\n"

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

		#printf $"$WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.vcf\n" >> $CFG
		printf $"\n=========> Features extraction: DONE"
	done
	INPUT=$WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.vcf


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
		 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME_SOM \
		 	-l $LISTAFEATURES_SOMATIC
			
		elif [ "$DESIGN" == "AMPLICON" ]
		then

			python $SCRIPT_PIPELINE/features_extraction_somatic.py \
		 	-m $VCF_MUTECT \
		 	-d $VCF_VARDICT \
		 	-v $VCF_VARSCAN_SNP \
		 	-i $VCF_VARSCAN_INDEL \
		 	-n $SAMPLE_NAME_NORM -t $SAMPLE_NAME_SOM -a \
		 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME_SOM \
		 	-l $LISTAFEATURES_SOMATIC
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


Features_extraction_cellfree () {

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
		 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME_SOM \
		 	-l $LISTAFEATURES_SOMATIC
			
		elif [ "$DESIGN" == "AMPLICON" ]
		then

			python $SCRIPT_PIPELINE/features_extraction_somatic.py \
		 	-m $VCF_MUTECT \
		 	-d $VCF_VARDICT \
		 	-v $VCF_VARSCAN_SNP \
		 	-i $VCF_VARSCAN_INDEL \
		 	-n $SAMPLE_NAME_NORM -t $SAMPLE_NAME_SOM -a \
		 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME_SOM \
		 	-l $LISTAFEATURES_SOMATIC
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

VARIANT_CALLING_CELLFREE () {
	
	cat $LOGHI/logo_variant.txt

	mkdir -p $WORKDIR/VARIANT_CALLING
	mkdir -p $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION
	
	rm -f $WORKDIR/PostVariantCalling.cfg
	CFG=$WORKDIR/PostVariantCalling.cfg

	cat $1 | while read line
	do
		IFS=$':' DIRS=(${line//$'\t'/:})
		dastampare=()

		for (( i=0 ; i<${#DIRS[@]}-2 ; i++ ))
		do
			BAM_SOM=${DIRS[i]}
			SAMPLE_NAME_SOM=${DIRS[i+1]}
			BAM_NORM=${DIRS[-2]}
			SAMPLE_NAME_NORM=${DIRS[-1]}
			((i+=1))
			
			Mutect2 $BAM_SOM $BAM_NORM $SAMPLE_NAME_SOM $SAMPLE_NAME_NORM

			VarScan2_somatic $BAM_SOM $BAM_NORM $SAMPLE_NAME_SOM $SAMPLE_NAME_NORM
		
			VarDict $BAM_SOM $BAM_NORM $SAMPLE_NAME_SOM $SAMPLE_NAME_NORM

			mv $BAM_SOM $STORAGE

			printf $"$VCF_MUTECT\t$VCF_VARDICT\t$VCF_VARSCAN_SNP\t$VCF_VARSCAN_INDEL\t$SAMPLE_NAME_SOM\t$SAMPLE_NAME_NORM\n" >> $CFG

			#dastampare+=("$VCF_MUTECT\t$VCF_VARDICT\t$VCF_VARSCAN_SNP\t$VCF_VARSCAN_INDEL\t$SAMPLE_NAME_SOM")

		done

		cp $BAM_NORM $STORAGE

		#dastampare+=("$SAMPLE_NAME_SOM")
		#bar=$(IFS=$'\t' ; echo "${dastampare[*]}")
		#echo -e "$bar">>$CFG

		
	done

}