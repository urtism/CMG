save_in_storage () {

	cat $1 | while read line
	do
		cp $line $STORAGE
		cp $line.idx $STORAGE
		cp ${line%.*}.bai $STORAGE
	done
}

bam_to_bed_to_list (){
	
	bedtools bamtobed -i $1 > $TARGETBED
	java -Xmx64g -jar $PICARD BedToIntervalList I=$TARGETBED O=$TARGET SD=$REF_DICT
}

samtools_mpileup (){


	if [ "$TARGETBED" == "" ]
	then
		samtools mpileup -B -q 1 -d 50000 -L 50000 $3 -f $REF \
		$1> $2
	else
		samtools mpileup -B -q 1 -d 50000 -L 50000 $3 -f $REF \
		-l $TARGETBED $1 $3> $2
	fi
}

samtools_vc () {
	
	printf $"\n =========>	Sample $2 => Variant Calling: Samtools\n\n"

	samtools_mpileup $1 $WORKDIR/VARIANT_CALLING/$2\_Samtools.tobcftools.vcf -u

	$BCFTOOLS call -mv -Ov $WORKDIR/VARIANT_CALLING/$2\_Samtools.tobcftools.vcf >  $WORKDIR/VARIANT_CALLING/$2\_Samtools.vcf

	$BCFTOOLS norm -m -both \
	-f $REF \
	$WORKDIR/VARIANT_CALLING/$2\_Samtools.vcf \
	> $WORKDIR/VARIANT_CALLING/$2\_Samtools.split.vcf

	VCF_SAMTOOLS=$WORKDIR/VARIANT_CALLING/$2\_Samtools.split.vcf

	mv $WORKDIR/VARIANT_CALLING/$2\_Samtools.tobcftools.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$2\_Samtools.vcf $DELETE

	printf $"\n =========>	Sample $2 => Variant Calling: Samtools: DONE\n\n"

}



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     GATK     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Mutect2 () {
	printf $"\n =========>	Sample $3 => Variant Calling: MuTect2\n\n"
	
	java -Xmx64g -jar $GATK -T MuTect2 \
	-R $REF \
	-I:tumor $1 \
	-I:normal $2 \
	-L $TARGET \
	-o $WORKDIR/VARIANT_CALLING/$3\_Sane_GATK.vcf 

	#--min_base_quality_score 30 \
	#--maxReadsInRegionPerSample 50000 \

	$BCFTOOLS norm -m -both \
	-f $REF \
	$WORKDIR/VARIANT_CALLING/$3\_Sane_GATK.vcf \
	> $WORKDIR/VARIANT_CALLING/$3\_Sane_GATK.split.vcf 

	VCF_MUTECT=$WORKDIR/VARIANT_CALLING/$3\_Sane_GATK.split.vcf

	printf $"\n =========>	Sample $3 => Variant Calling: MuTect2: DONE\n\n"
}


HaplotypeCaller () {

	printf $"\n =========>	Sample $SAMPLE_NAME => Variant Calling: Haplotype Caller\n\n"

	if [ "$TARGET" == "" ]
	then
		java -Xmx64g -jar $GATK -T HaplotypeCaller \
		-R $REF \
		-I $1 \
		-o $WORKDIR/VARIANT_CALLING/GVCF/$SAMPLE_NAME.g.vcf \
		-ERC GVCF \
		--doNotRunPhysicalPhasing \
		--min_base_quality_score 1
	else
		java -Xmx64g -jar $GATK -T HaplotypeCaller \
		-R $REF \
		-I $1 \
		-o $WORKDIR/VARIANT_CALLING/GVCF/$SAMPLE_NAME.g.vcf \
		-ERC GVCF \
		--doNotRunPhysicalPhasing \
		-L $TARGET \
		--min_base_quality_score 1 \
		--minReadsPerAlignmentStart 2 \
		--standard_min_confidence_threshold_for_calling 5.0 


	fi


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

FreeBayes_multisample () {

	printf $'\n =========>	Variant Calling with FreeBayes\n\n'
	if [ "$TARGETBED" == "" ]
		then
		$FREEBAYES -f $REF \
		-L $1 \
		-K \
		-J \
		-s $2 \
		--use-best-n-alleles 6 \
		--genotype-qualities \
		--report-genotype-likelihood-max \
		--allele-balance-priors-off \
		# -m 0 \
		# -q 0 \
		# -R 0 \
		# -Y 0 \
		# -Q 1 \
		# -F 0.001 \
		# -C 1 \
		> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf
	else
		$FREEBAYES -f $REF \
		-L $1 \
		-K \
		-J \
		-s $2 \
		--use-best-n-alleles 6 \
		--genotype-qualities \
		--report-genotype-likelihood-max \
		--allele-balance-priors-off \
		-t $TARGETBED > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf
		# -m 0 \
		# -q 0 \
		# -R 0 \
		# -Y 0 \
		# -Q 1 \
		# -F 0.001 \
		# -C 1 \
		#> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_FreeBayes.vcf
	fi

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



FreeBayes_singlesample () {

	printf $'\n =========>	Variant Calling with FreeBayes\n\n'
	if [ "$TARGETBED" == "" ]
		then
		$FREEBAYES -f $REF \
		-b $1 \
		-K \
		-J \
		--genotype-qualities \
		--report-genotype-likelihood-max \
		--allele-balance-priors-off \
		--genotype-qualities \
		--report-genotype-likelihood-max \
		--allele-balance-priors-off \
		-m 0 \
		-q 0 \
		-R 0 \
		-Y 0 \
		-Q 1 \
		-F 0.001 \
		-C 1 > $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.vcf
	else
		$FREEBAYES -f $REF \
		-b $1 \
		-K \
		-J \
		--genotype-qualities \
		--report-genotype-likelihood-max \
		--allele-balance-priors-off \
		-t $TARGETBED \
		--genotype-qualities \
		--report-genotype-likelihood-max \
		--allele-balance-priors-off \
		-m 0 \
		-q 0 \
		-R 0 \
		-Y 0 \
		-Q 1 \
		-F 0.001 \
		-C 1 > $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.vcf
	fi

	python  $SCRIPT_PIPELINE/header_fix.py -v F -f $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.vcf \
	> $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.fix.vcf

 	$BCFTOOLS norm -m -both \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.fix.vcf \
 	> $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.split.vcf

 	mv $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.fix.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.vcf $DELETE

	VCF_FREEBAYES=$WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.split.vcf

	printf $'\n =========>	Variant Calling with FreeBayes: DONE\n\n'
	
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARSCAN2     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

VarScan2_germline_multisample () {	

	printf $'\n =========>	Variant Calling with VarScan2\n\n'
	if [ "$TARGETBED" == "" ]
	then
		samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
		-b $1 > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO.mpileup
	else
		samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
		-l $TARGETBED -b $1 > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO.mpileup
	fi

	java -jar -Xmx64g $VARSCAN mpileup2snp $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO.mpileup \
	--min-coverage 1 \
	--min-var-freq 0.001 \
	--pvalue 0.5 \
	--output-vcf 1 \
	--min-reads2 1 \
	--min-avg-qual 0 \
	--strand-filter 0 \
	--vcf-sample-list $2 > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf

	java -jar -Xmx64g $VARSCAN mpileup2indel $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO.mpileup \
	--min-coverage 1 \
	--min-var-freq 0.001 \
	--pvalue 0.5 \
	--output-vcf 1 \
	--min-reads2 1 \
	--min-avg-qual 0 \
	--strand-filter 0 \
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

VarScan2_germline_singlesample () {	

	printf $'\n =========>	Variant Calling with VarScan2\n\n'
	# if [ "$TARGETBED" == "" ]
	# then
	# 	samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
	# 	-b $1 > $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME.mpileup
	# else
	# 	samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
	# 	-l $TARGETBED -b $1 > $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME.mpileup
	# fi

	# java -jar -Xmx64g $VARSCAN mpileup2snp $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO.mpileup \
	# --min-coverage 10 \
	# --min-var-freq 0.20 \
	# --pvalue 0.05 \
	# --output-vcf 1 \
	# --vcf-sample-list $2 > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_snp.vcf

	# java -jar -Xmx64g $VARSCAN mpileup2indel $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO.mpileup \
	# --min-coverage 10 \
	# --min-var-freq 0.10 \
	# --pvalue 0.1 \
	# --output-vcf 1 \
	# --vcf-sample-list $2 > $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_VarScan_Indel.vcf

	java -jar -Xmx64g $VARSCAN mpileup2snp $1 \
	--min-coverage 0 \
	--min-var-freq 0.001 \
	--pvalue 0.5 \
	--output-vcf 1 \
	--min-reads2 0 \
	--min-avg-qual 0 \
	--strand-filter 0 \
	--vcf-sample-list $2 > $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_snp.vcf

	java -jar -Xmx64g $VARSCAN mpileup2indel $1 \
	--min-coverage 0 \
	--min-var-freq 0.001 \
	--pvalue 0.5 \
	--output-vcf 1 \
	--min-reads2 0 \
	--min-avg-qual 0 \
	--strand-filter 0 \
	--vcf-sample-list $2 > $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_Indel.vcf

	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_Indel.vcf
	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_snp.vcf

	bgzip $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_snp.vcf
	bgzip $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_Indel.vcf
	tabix $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_snp.vcf.gz
	tabix $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_Indel.vcf.gz

	vcf-concat $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_snp.vcf.gz $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_Indel.vcf.gz >$WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan.Merge.vcf
	vcf-sort -c $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan.Merge.vcf > $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan.sort.vcf

	$BCFTOOLS norm -m -both \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan.sort.vcf \
 	> $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan.split.vcf

 	VCF_VARSCAN=$WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan.split.vcf

	mv $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan.Merge.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan.sort.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_snp.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_Indel.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_snp.vcf.gz.tbi $DELETE
	mv $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_VarScan_Indel.vcf.gz.tbi $DELETE
	mv $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME.mpileup $DELETE

	printf $'\n =========>	Variant Calling with VarScan2: DONE\n'

}


VarScan2_somatic () {
	
	printf $"\n =========>	Sample $3 => Variant Calling: VarScan2\n\n"
	
	samtools mpileup -f $REF -l $TARGETBED -d 50000 -L 50000 -q 1 -B $2 $1 > $WORKDIR/VARIANT_CALLING/$3\_Sane.mpileup 
 	#-Q 30	
	java -jar -Xmx64g $VARSCAN somatic $WORKDIR/VARIANT_CALLING/$3\_Sane.mpileup \
	$WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan \
	--min-var-freq 0.005 \
	--output-vcf 1 \
	--mpileup 1

	python $SCRIPT_PIPELINE/header_fix.py -v V -f $WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.snp.vcf \
	> $WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.fix.snp.vcf

	$BCFTOOLS norm -m -both \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$3\_Sane_VarScan.fix.snp.vcf \
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


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARDICT     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

VarDict () {	

	printf $"\n =========> Variant Calling: VarDict\n\n"
	
	$VARDICT -G $REF -f 0.005 -N $3 -b "$1|$2" \
	-z -F 0 -c 1 -S 2 -E 3 -g 4 $TARGETBED | ~/NGS_TOOLS/VarDictJava-master/VarDict/testsomatic.R | ~/NGS_TOOLS/VarDictJava-master/VarDict/var2vcf_paired.pl \
	-N "$3|$4" -f 0.005 > $WORKDIR/VARIANT_CALLING/$3\_Sane_VarDict.vcf
	#-q 30
	VCF_VARDICT=$WORKDIR/VARIANT_CALLING/$3\_Sane_VarDict.vcf
	
	printf $"\n =========>	Sample $3 => Variant Calling: VarDict: DONE\n\n"	
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     SCALPEL     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

scalpel_germline_singlesample () {

	printf $"\n =========>	Sample $2 => Variant Calling: Scalpel\n\n"

	$SCALPEL --single \
	--bam $1 \
	--bed $TARGETBED \
	--ref $REF \
	--numprocs 8 \
	--dir $WORKDIR/VARIANT_CALLING/SCALPEL

	mv $WORKDIR/VARIANT_CALLING/SCALPEL/variants.indel.vcf $WORKDIR/VARIANT_CALLING/$2\_Scalpel.vcf

	VCF_SCALPEL=$WORKDIR/VARIANT_CALLING/$2\_Scalpel.vcf

	printf $"\n =========>	Sample $2 => Variant Calling: Scalpel: DONE\n\n"

}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     SNVER     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

SNVer_germline_singlesample () {

	printf $"\n =========>	Sample $2 => Variant Calling: SNVer\n\n"

	java -jar -Xmx64g $SNVER_INDIVIDUAL \
	-i $1 \
	-l $TARGETBED \
	-r $REF \
	-o $WORKDIR/VARIANT_CALLING/SNVER/$2

	bgzip $WORKDIR/VARIANT_CALLING/SNVER/$2.filter.vcf
	bgzip $WORKDIR/VARIANT_CALLING/SNVER/$2.indel.filter.vcf
	tabix $WORKDIR/VARIANT_CALLING/SNVER/$2.filter.vcf.gz
	tabix $WORKDIR/VARIANT_CALLING/SNVER/$2.indel.filter.vcf.gz

	vcf-concat $WORKDIR/VARIANT_CALLING/SNVER/$2.filter.vcf.gz $WORKDIR/VARIANT_CALLING/SNVER/$2.indel.filter.vcf.gz > $WORKDIR/VARIANT_CALLING/SNVER/$2\_SNVer.vcf
	vcf-sort -c $WORKDIR/VARIANT_CALLING/SNVER/$2\_SNVer.vcf > $WORKDIR/VARIANT_CALLING/$2\_SNVer.sort.vcf

	mv $WORKDIR/VARIANT_CALLING/SNVER/$2\_SNVer.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/SNVER/$2.filter.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/SNVER/$2.indel.filter.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/SNVER/$2.filter.vcf.gz.tbi $DELETE
	mv $WORKDIR/VARIANT_CALLING/SNVER/$2.indel.filter.vcf.gz.tbi $DELETE

	VCF_SNVER=$WORKDIR/VARIANT_CALLING/$2\_SNVer.sort.vcf

	printf $"\n =========>	Sample $2 => Variant Calling: SNVer: DONE\n\n"

}

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     MERGE VCF    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Vcf_merge () {

	printf $"\n =========>	Merging VCFs: GATK,FreeBayes,VarScan\n\n"
	
	python $SCRIPT_PIPELINE/merge_vcf.py \
	-g $1 \
	-f $2 \
	-v $3 \
	-o $(dirname $1)/$DATA\_$PANNELLO.Merged.vcf

	VCF_MERGED=$(dirname $1)/$DATA\_$PANNELLO.Merged.vcf

	printf $"\n =========>	Merging VCFs: DONE\n\n"
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
	
	FreeBayes_multisample $WORKDIR/Bam_list.txt $WORKDIR/Sample_list.txt

	VarScan2_germline_multisample $WORKDIR/Bam_list.txt $WORKDIR/Sample_list.txt

	save_in_storage $WORKDIR/Bam_list.txt

	#Vcf_merge $VCF_GATK $VCF_FREEBAYES $VCF_VARSCAN 

	printf $"$VCF_GATK\t$VCF_FREEBAYES\t$VCF_VARSCAN\n" >> $CFG
	#printf $"$VCF_MERGED\n" >> $CFG
	
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
		mv ${BAM_SOM%.*}.bai $STORAGE
		mv $BAM_NORM $STORAGE
		mv ${BAM_NORM%.*}.bai $STORAGE

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
			mv ${BAM_SOM%.*}.bai $STORAGE

			printf $"$VCF_MUTECT\t$VCF_VARDICT\t$VCF_VARSCAN_SNP\t$VCF_VARSCAN_INDEL\t$SAMPLE_NAME_SOM\t$SAMPLE_NAME_NORM\n" >> $CFG

		done



		cp $BAM_NORM $STORAGE
		cp ${BAM_NORM%.*}.bai $STORAGE

	done

}



VARIANT_CALLING_SANGER () {
	
	cat $LOGHI/logo_variant.txt

	mkdir -p $WORKDIR/VARIANT_CALLING
	mkdir -p $WORKDIR/TARGETBED
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

		TARGETBED=$WORKDIR/TARGETBED/$SAMPLE_NAME.bed
		TARGET=$WORKDIR/TARGETBED/$SAMPLE_NAME.list

		bam_to_bed_to_list $BAM
		
		samtools_mpileup $BAM $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME.mpileup

		samtools_vc $BAM
		
		FreeBayes_singlesample $BAM

		VarScan2_germline_singlesample $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME.mpileup $SAMPLE_NAME

		printf $"$BAM\n" >> $WORKDIR/Bam_list.txt
		printf $"$VCF_SAMTOOLS\t$VCF_FREEBAYES\t$VCF_VARSCAN\n" >> $CFG

	done
	
	save_in_storage $WORKDIR/Bam_list.txt
	
	
}



VARIANT_CALLING_SINGLE_SAMPLE () {
	
	cat $LOGHI/logo_variant.txt

	mkdir -p $WORKDIR/VARIANT_CALLING
	mkdir -p $WORKDIR/TARGETBED
	mkdir -p $WORKDIR/VARIANT_CALLING/GVCF
	mkdir -p $WORKDIR/VARIANT_CALLING/FEATURES_EXTRACTION
	mkdir -p $WORKDIR/VARIANT_CALLING/SCALPEL
	mkdir -p $WORKDIR/VARIANT_CALLING/SNVER
	
	rm -f $WORKDIR/Bam_list.txt
	rm -f $WORKDIR/Sample_list.txt
	rm -f $WORKDIR/gvcf.list
	rm -f $WORKDIR/PostVariantCalling.cfg
	CFG=$WORKDIR/PostVariantCalling.cfg
	

	cat $1 | while read line
	do

		BAM=$(echo "$line" | cut -f1)
		SAMPLE_NAME=$(echo "$line" | cut -f2)

		samtools_vc $BAM $SAMPLE_NAME
		scalpel_germline_singlesample $BAM $SAMPLE_NAME
		SNVer_germline_singlesample $BAM $SAMPLE_NAME

		printf $"$BAM\n" >> $WORKDIR/Bam_list.txt
		printf $"$VCF_SAMTOOLS\t$VCF_SCALPEL\t$VCF_SNVER\n" >> $CFG

	done
	
	#save_in_storage $WORKDIR/Bam_list.txt
	
	
}