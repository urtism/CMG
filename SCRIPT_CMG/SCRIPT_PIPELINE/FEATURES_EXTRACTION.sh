#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::      iEVA      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.

iEVA_germline(){
	
	rm -f $WORKDIR/PostiEVA.cfg
	CFG=$WORKDIR/PostiEVA.cfg
	
	cat $1 | while read line
	do
		printf $"\n=========> Extracting Variant Info: iEVA\n\n"

		VCF_MERGED=$(echo "$line" | cut -f1)
		VCF_iEVA=${VCF_MERGED%.*}.iEVA.vcf

		python iEVA \
		--input $VCF_MERGED \
		--reference $REF \
		--list $WORKDIR/Bam_list.txt \
		--outfile $VCF_iEVA \
		--SimpleRepeat \
		--SimpleRepeatLength \
		--PseudoNucleotidesComposition \
		--RepeatMasker \
		--gcContent \
		--VariantClass \
		--StrandBiasReads \
		--UnMappedReads \
		--MappingQualityZero \
		--NotPrimaryAlignment \
		--SupplementaryAlignment \
		--NotPairedReads \
		--NotProperPairedReads \
		--AlignmentScore \
		--NumberTotalDupReads \
		--NumberReadDupRef \
		--NumberReadDupAlt \
		--DuplicateReference \
		--DuplicateAlternate \
		--DeltaDuplicate \
		--iEvaDepth \
		--iAlleleDepth \
		--ReadRef \
		--ReadAlt \
		--MeanRefQscore \
		--MeanAltQscore \
		--TotalDPUnfilter \
		--NumberClippedReadsRef \
		--NumberClippedReadsAlt \
		--ClippedReadsRef \
		--ClippedReadsAlt

		printf $"$VCF_iEVA\n">>$CFG

		printf $"\n=========> Extracting Variant Info: DONE\n\n"

	done
	#INPUT=${VCF_MERGED%.*}.iEVA.vcf
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
		SAMPLE_NAME=$(echo "$line" | cut -f4)

		printf $"\n=========> Features extraction: $SAMPLE_NAME\n\n"

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

		mv $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/TOTAL.vcf $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME\_TOTAL.vcf
		cp $VCF_GATK $OUT
		cp $VCF_FREEBAYES $OUT
		cp $VCF_VARSCAN $OUT
		mv $VCF_GATK $STORAGE
		mv $VCF_FREEBAYES $STORAGE
		mv $VCF_VARSCAN $STORAGE
		
		save_in_storage $WORKDIR/gvcf.list

		printf $"$WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME\_TOTAL.vcf\t$SAMPLE_NAME\n" >> $CFG
		printf $"\n=========> Features extraction: DONE"
	done
	INPUT=$WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.vcf


}


features_extraction_germline_iEVA () {

	rm -f $WORKDIR/PostFeatureExtraction.cfg
	CFG=$WORKDIR/PostFeatureExtraction.cfg
	
	cat $1 | while read line
	do
		VCF=$(echo "$line" | cut -f1)
		
		printf $"\n=========> Features extraction\n\n"

		if [ "$DESIGN" == "ENRICHMENT" ]
		then

			python $SCRIPT_PIPELINE/features_extraction_germline_iEVA.py \
		 	-m $VCF \
		 	-l $LISTAFEATURES_GERMLINE \
		 	-F -s \
		 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION \
		    -G $WORKDIR/VARIANT_CALLING/GVCF

		elif [ "$DESIGN" == "AMPLICON" ]
		then

			python $SCRIPT_PIPELINE/features_extraction_germline_iEVA.py \
		 	-m $VCF \
		 	-l $LISTAFEATURES_GERMLINE \
		 	-F -a -s \
		 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION \
		    -G $WORKDIR/VARIANT_CALLING/GVCF

		fi

		mv $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/TOTAL.vcf $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.vcf
		#cp $VCF_GATK $OUT
		#cp $VCF_FREEBAYES $OUT
		#cp $VCF_VARSCAN $OUT
		#mv $VCF_GATK $STORAGE
		#mv $VCF_FREEBAYES $STORAGE
		#mv $VCF_VARSCAN $STORAGE
		
		save_in_storage $WORKDIR/gvcf.list

		#printf $"$WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.vcf\n" >> $CFG
		printf $"\n=========> Features extraction: DONE"
	done
	INPUT=$WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$DATA\_$PANNELLO\_TOTAL.vcf


}


Features_extraction_sanger () {

	rm -f $WORKDIR/PostFeatureExtraction.cfg
	CFG=$WORKDIR/PostFeatureExtraction.cfg
	
	cat $1 | while read line
	do
		VCF_SAMTOOLS=$(echo "$line" | cut -f1)
		VCF_FREEBAYES=$(echo "$line" | cut -f2)
		VCF_VARSCAN=$(echo "$line" | cut -f3)
		
		printf $"\n=========> Features extraction\n\n"
		
		python $SCRIPT_PIPELINE/features_extraction_sanger.py \
	 	-g $VCF_SAMTOOLS \
	 	-f $VCF_FREEBAYES \
	 	-v $VCF_VARSCAN \
	 	-F \
	 	-o $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION

		cp $VCF_SAMTOOLS $OUT
		cp $VCF_FREEBAYES $OUT
		cp $VCF_VARSCAN $OUT
		mv $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION/$SAMPLE_NAME.tsv $OUT
		#mv $VCF_SAMTOOLS $STORAGE
		#mv $VCF_FREEBAYES $STORAGE
		#mv $VCF_VARSCAN $STORAGE

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
