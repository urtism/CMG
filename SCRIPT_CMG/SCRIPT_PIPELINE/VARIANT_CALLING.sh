save_in_storage () {

	cat $1 | while read line
	do
		cp $line $STORAGE
		if [ -e $line.idx ]; then
			cp $line.idx $STORAGE
		fi
		if [ -e ${line%.*}.bai ]; then
			cp ${line%.*}.bai $STORAGE
		fi
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

	#$BCFTOOLS norm -m -both \
	#-f $REF \
	#$WORKDIR/VARIANT_CALLING/$2\_Samtools.vcf \
	#> $WORKDIR/VARIANT_CALLING/$2\_Samtools.norm.vcf

	VCF_SAMTOOLS=$WORKDIR/VARIANT_CALLING/$2\_Samtools.norm.vcf

	mv $WORKDIR/VARIANT_CALLING/$2\_Samtools.tobcftools.vcf $DELETE
	#mv $WORKDIR/VARIANT_CALLING/$2\_Samtools.vcf $DELETE

	printf $"\n =========>	Sample $2 => Variant Calling: Samtools: DONE\n\n"

}

samtools_vc_multisample () {

	printf $'\n =========>	Variant Calling with Samtools\n\n'

	if [ "$TARGETBED" == "" ]
	then
		samtools mpileup -q 1 -d 50000 -L 50000 -uf $REF -b $1 | $BCFTOOLS call -mv -f GQ,GP --threads 8  -o $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Samtools.vcf
	else
		samtools mpileup -q 1 -d 50000 -L 50000 -uf $REF -l $TARGETBED -b $1 | $BCFTOOLS call -mv -f GQ,GP --threads 8 -o $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Samtools.vcf 
	fi

	python  $SCRIPT_PIPELINE/header_fix.py -v S -f $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Samtools.vcf \
	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Samtools.fix.vcf 


	$BCFTOOLS norm -m -both \
	-f $REF \
	$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Samtools.fix.vcf \
	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Samtools.norm.vcf


	VCF_SAMTOOLS=$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Samtools.norm.vcf

	printf $'\n =========>	Variant Calling with Samtools: DONE\n\n'
}



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     GATK     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Mutect2 () {
	printf $"\n =========>	Sample $3 => Variant Calling: MuTect2\n\n"
	
	java -Xmx64g -jar $GATK -T MuTect2 \
	-R $REF \
	-I:tumor $1 \
	-I:normal $2 \
	-L $TARGET \
	-o $WORKDIR/VARIANT_CALLING/$3\_$4\_GATK.vcf 

	#--min_base_quality_score 30 \
	#--maxReadsInRegionPerSample 50000 \

	$BCFTOOLS norm -m -both \
	-f $REF \
	$WORKDIR/VARIANT_CALLING/$3\_$4\_GATK.vcf \
	> $WORKDIR/VARIANT_CALLING/$3\_$4\_GATK.norm.vcf 

	VCF_MUTECT=$WORKDIR/VARIANT_CALLING/$3\_$4\_GATK.norm.vcf

	printf $"\n =========>	Sample $3 => Variant Calling: MuTect2: DONE\n\n"
}


HaplotypeCaller () {

	printf $"\n =========>	Sample $2 => Variant Calling: Haplotype Caller\n\n"

	if [ "$TARGET" == "" ]
	then
		java -Xmx64g -jar $GATK -T HaplotypeCaller \
		-R $REF \
		-I $1 \
		-o $WORKDIR/VARIANT_CALLING/GVCF/$2.g.vcf \
		-ERC GVCF \
		--doNotRunPhysicalPhasing \
		--min_base_quality_score 1
	else
		java -Xmx64g -jar $GATK -T HaplotypeCaller \
		-R $REF \
		-I $1 \
		-o $WORKDIR/VARIANT_CALLING/GVCF/$2.g.vcf \
		-ERC GVCF \
		--doNotRunPhysicalPhasing \
		-L $TARGET \
		--min_base_quality_score 1 \
		--minReadsPerAlignmentStart 2 \
		--standard_min_confidence_threshold_for_calling 1.0 


	fi


	printf $"$WORKDIR/VARIANT_CALLING/GVCF/$2.g.vcf\n" >> $WORKDIR/gvcf.list
	printf $"\n =========>	Sample $2 => Variant Calling: Haplotype Caller: DONE\n\n"
	
}


GenotypeGVCFs () {
	
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_multi.txt

	printf $'\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling\n\n'

	java -jar -Xmx64g $GATK -T GenotypeGVCFs \
	-R $REF \
	-V:VCF $1 \
	-o $WORKDIR/VARIANT_CALLING/$2\_GATK.vcf

	python $SCRIPT_PIPELINE/header_fix.py -v G -f $WORKDIR/VARIANT_CALLING/$2\_GATK.vcf \
	> $WORKDIR/VARIANT_CALLING/$2\_GATK.fix.vcf

	$BCFTOOLS norm -m -both -D \
	-f $REF \
	$WORKDIR/VARIANT_CALLING/$2\_GATK.fix.vcf \
	> $WORKDIR/VARIANT_CALLING/$2\_GATK.norm.vcf

	mv $WORKDIR/VARIANT_CALLING/$2\_GATK.fix.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$2\_GATK.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$2\_GATK.vcf.idx $DELETE

	VCF_GATK=$WORKDIR/VARIANT_CALLING/$2\_GATK.norm.vcf

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
		-m 0 \
		-q 0 \
		-R 0 \
		-Y 0 \
		-Q 1 \
		-F 0.001 \
		-C 1 \
		> $WORKDIR/VARIANT_CALLING/$3\_FreeBayes.vcf
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
		-t $TARGETBED \
		-m 0 \
		-q 0 \
		-R 0 \
		-Y 0 \
		-Q 1 \
		-F 0.001 \
		-C 1 \
		> $WORKDIR/VARIANT_CALLING/$3\_FreeBayes.vcf
	fi

	python  $SCRIPT_PIPELINE/header_fix.py -v F -f $WORKDIR/VARIANT_CALLING/$3\_FreeBayes.vcf \
	> $WORKDIR/VARIANT_CALLING/$3\_FreeBayes.fix.vcf

 	$BCFTOOLS norm -m -both -D \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$3\_FreeBayes.fix.vcf \
 	> $WORKDIR/VARIANT_CALLING/$3\_FreeBayes.norm.vcf

 	mv $WORKDIR/VARIANT_CALLING/$3\_FreeBayes.fix.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$3\_FreeBayes.vcf $DELETE

	VCF_FREEBAYES=$WORKDIR/VARIANT_CALLING/$3\_FreeBayes.norm.vcf

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
		-m 0 \
		-q 0 \
		-R 0 \
		-Y 0 \
		-Q 1 \
		-F 0.01 \
		-C 1 > $WORKDIR/VARIANT_CALLING/$2\_FreeBayes.vcf
	else
		$FREEBAYES -f $REF \
		-b $1 \
		-K \
		-J \
		--genotype-qualities \
		--report-genotype-likelihood-max \
		--allele-balance-priors-off \
		-t $TARGETBED \
		-m 0 \
		-q 0 \
		-R 0 \
		-Y 0 \
		-Q 1 \
		-F 0.01 \
		-C 1 > $WORKDIR/VARIANT_CALLING/$2\_FreeBayes.vcf
	fi

	python  $SCRIPT_PIPELINE/header_fix.py -v F -f $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.vcf \
	> $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.fix.vcf

 	$BCFTOOLS norm -m -both -D \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.fix.vcf \
 	> $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.norm.vcf

 	mv $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.fix.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.vcf $DELETE

	VCF_FREEBAYES=$WORKDIR/VARIANT_CALLING/$SAMPLE_NAME\_FreeBayes.norm.vcf

	printf $'\n =========>	Variant Calling with FreeBayes: DONE\n\n'
	
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARSCAN2     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

VarScan2_germline_multisample () {	

	printf $'\n =========>	Variant Calling with VarScan2\n\n'
	if [ "$TARGETBED" == "" ]
	then
		samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
		-b $1 > $WORKDIR/VARIANT_CALLING/$3.mpileup
	else
		samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
		-l $TARGETBED -b $1 > $WORKDIR/VARIANT_CALLING/$3.mpileup
	fi

	java -jar -Xmx64g $VARSCAN mpileup2snp $WORKDIR/VARIANT_CALLING/$3.mpileup \
	--min-coverage 1 \
	--min-var-freq 0.01 \
	--pvalue 0.5 \
	--output-vcf 1 \
	--min-reads2 1 \
	--min-avg-qual 0 \
	--strand-filter 0 \
	--vcf-sample-list $2 > $WORKDIR/VARIANT_CALLING/$3\_VarScan_snp.vcf

	java -jar -Xmx64g $VARSCAN mpileup2indel $WORKDIR/VARIANT_CALLING/$3.mpileup \
	--min-coverage 1 \
	--min-var-freq 0.01 \
	--pvalue 0.5 \
	--output-vcf 1 \
	--min-reads2 1 \
	--min-avg-qual 0 \
	--strand-filter 0 \
	--vcf-sample-list $2 > $WORKDIR/VARIANT_CALLING/$3\_VarScan_Indel.vcf

	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $WORKDIR/VARIANT_CALLING/$3\_VarScan_Indel.vcf
	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $WORKDIR/VARIANT_CALLING/$3\_VarScan_snp.vcf

	bgzip $WORKDIR/VARIANT_CALLING/$3\_VarScan_snp.vcf
	bgzip $WORKDIR/VARIANT_CALLING/$3\_VarScan_Indel.vcf
	tabix $WORKDIR/VARIANT_CALLING/$3\_VarScan_snp.vcf.gz
	tabix $WORKDIR/VARIANT_CALLING/$3\_VarScan_Indel.vcf.gz

	vcf-concat $WORKDIR/VARIANT_CALLING/$3\_VarScan_snp.vcf.gz $WORKDIR/VARIANT_CALLING/$3\_VarScan_Indel.vcf.gz >$WORKDIR/VARIANT_CALLING/$3\_VarScan.Merge.vcf
	vcf-sort -c $WORKDIR/VARIANT_CALLING/$3\_VarScan.Merge.vcf > $WORKDIR/VARIANT_CALLING/$3\_VarScan.sort.vcf

	$BCFTOOLS norm -m -both -D \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$3\_VarScan.sort.vcf \
 	> $WORKDIR/VARIANT_CALLING/$3\_VarScan.norm.vcf

 	VCF_VARSCAN=$WORKDIR/VARIANT_CALLING/$3\_VarScan.norm.vcf

	mv $WORKDIR/VARIANT_CALLING/$3\_VarScan.Merge.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$3\_VarScan.sort.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$3\_VarScan_snp.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/$3\_VarScan_Indel.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/$3\_VarScan_snp.vcf.gz.tbi $DELETE
	mv $WORKDIR/VARIANT_CALLING/$3\_VarScan_Indel.vcf.gz.tbi $DELETE
	mv $WORKDIR/VARIANT_CALLING/$3.mpileup $DELETE

	printf $'\n =========>	Variant Calling with VarScan2: DONE\n'

}

VarScan2_germline_singlesample () {	

	printf $'\n =========>	Variant Calling with VarScan2\n\n'
	
	if [ "$TARGETBED" == "" ]
	then
		samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
		$1 -b $3 > $WORKDIR/VARIANT_CALLING/$2.mpileup
	else
		samtools mpileup -B -q 1 -d 50000 -L 50000 -f $REF \
		-l $TARGETBED $1 -b $3 > $WORKDIR/VARIANT_CALLING/$2.mpileup
	fi

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

	java -jar -Xmx64g $VARSCAN mpileup2snp $WORKDIR/VARIANT_CALLING/$2.mpileup \
	--min-coverage 0 \
	--min-var-freq 0.01 \
	--pvalue 0.5 \
	--output-vcf 1 \
	--min-reads2 0 \
	--min-avg-qual 0 \
	--strand-filter 0 > $WORKDIR/VARIANT_CALLING/$2\_VarScan_snp.vcf

	java -jar -Xmx64g $VARSCAN mpileup2indel $WORKDIR/VARIANT_CALLING/$2.mpileup \
	--min-coverage 0 \
	--min-var-freq 0.01 \
	--pvalue 0.5 \
	--output-vcf 1 \
	--min-reads2 0 \
	--min-avg-qual 0 \
	--strand-filter 0 > $WORKDIR/VARIANT_CALLING/$2\_VarScan_Indel.vcf

	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $WORKDIR/VARIANT_CALLING/$2\_VarScan_Indel.vcf
	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $WORKDIR/VARIANT_CALLING/$2\_VarScan_snp.vcf

	bgzip $WORKDIR/VARIANT_CALLING/$2\_VarScan_snp.vcf
	bgzip $WORKDIR/VARIANT_CALLING/$2\_VarScan_Indel.vcf
	tabix $WORKDIR/VARIANT_CALLING/$2\_VarScan_snp.vcf.gz
	tabix $WORKDIR/VARIANT_CALLING/$2\_VarScan_Indel.vcf.gz

	vcf-concat $WORKDIR/VARIANT_CALLING/$2\_VarScan_snp.vcf.gz $WORKDIR/VARIANT_CALLING/$2\_VarScan_Indel.vcf.gz >$WORKDIR/VARIANT_CALLING/$2\_VarScan.Merge.vcf
	vcf-sort -c $WORKDIR/VARIANT_CALLING/$2\_VarScan.Merge.vcf > $WORKDIR/VARIANT_CALLING/$2\_VarScan.sort.vcf

	$BCFTOOLS norm -m -both -D \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$2\_VarScan.sort.vcf \
 	> $WORKDIR/VARIANT_CALLING/$2\_VarScan.norm.vcf

 	VCF_VARSCAN=$WORKDIR/VARIANT_CALLING/$2\_VarScan.norm.vcf

	mv $WORKDIR/VARIANT_CALLING/$2\_VarScan.Merge.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$2\_VarScan.sort.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$2\_VarScan_snp.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/$2\_VarScan_Indel.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/$2\_VarScan_snp.vcf.gz.tbi $DELETE
	mv $WORKDIR/VARIANT_CALLING/$2\_VarScan_Indel.vcf.gz.tbi $DELETE
	mv $WORKDIR/VARIANT_CALLING/$2.mpileup $DELETE

	printf $'\n =========>	Variant Calling with VarScan2: DONE\n'

}


VarScan2_somatic () {
	
	printf $"\n =========>	Sample $3 => Variant Calling: VarScan2\n\n"
	
	samtools mpileup -f $REF -l $TARGETBED -d 50000 -L 50000 -q 1 -B $2 $1 > $WORKDIR/VARIANT_CALLING/$3\_$4\.mpileup 
 	#-Q 30	
	java -jar -Xmx64g $VARSCAN somatic $WORKDIR/VARIANT_CALLING/$3\_$4\.mpileup \
	$WORKDIR/VARIANT_CALLING/$3\_$4\_VarScan \
	--min-var-freq 0.005 \
	--output-vcf 1 \
	--mpileup 1

	python $SCRIPT_PIPELINE/header_fix.py -v V -f $WORKDIR/VARIANT_CALLING/$3\_$4\_VarScan.snp.vcf \
	> $WORKDIR/VARIANT_CALLING/$3\_$4\_VarScan.fix.snp.vcf

	$BCFTOOLS norm -m -both -D \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$3\_$4\_VarScan.fix.snp.vcf \
 	> $WORKDIR/VARIANT_CALLING/$3\_$4\_VarScan.norm.snp.vcf

 	python $SCRIPT_PIPELINE/header_fix.py -v V -f $WORKDIR/VARIANT_CALLING/$3\_$4\_VarScan.indel.vcf \
	> $WORKDIR/VARIANT_CALLING/$3\_$4\_VarScan.fix.indel.vcf

 	$BCFTOOLS norm -m -both -D \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$3\_$4\_VarScan.fix.indel.vcf \
 	> $WORKDIR/VARIANT_CALLING/$3\_$4\_VarScan.norm.indel.vcf

	VCF_VARSCAN_SNP=$WORKDIR/VARIANT_CALLING/$3\_$4\_VarScan.norm.snp.vcf
	VCF_VARSCAN_INDEL=$WORKDIR/VARIANT_CALLING/$3\_$4\_VarScan.norm.indel.vcf
	
	mv $WORKDIR/VARIANT_CALLING/$3\_$4\.mpileup $DELETE
	printf $"\n =========>	Sample $3 => Variant Calling: Varscan2: DONE\n\n"	
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARDICT     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

VarDict () {	

	printf $"\n =========> Variant Calling: VarDict\n\n"
	
	$VARDICT -G $REF -f 0.005 -N "$3|$4" -b "$1|$2" \
	-z -F 0 -c 1 -S 2 -E 3 -g 4 $TARGETBED | ~/NGS_TOOLS/VarDictJava-master/VarDict/testsomatic.R | ~/NGS_TOOLS/VarDictJava-master/VarDict/var2vcf_paired.pl \
	-N "$3|$4" -f 0.005 > $WORKDIR/VARIANT_CALLING/$3\_$4\_VarDict.vcf
	#-q 30
	VCF_VARDICT=$WORKDIR/VARIANT_CALLING/$3\_$4\_VarDict.vcf
	
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
	
	echo $1 
	echo $2

	sed "s/sample/$2/g" $WORKDIR/VARIANT_CALLING/SCALPEL/variants.indel.vcf > $WORKDIR/VARIANT_CALLING/$2\_Scalpel.vcf

	python  $SCRIPT_PIPELINE/header_fix.py -v L -f $WORKDIR/VARIANT_CALLING/$2\_Scalpel.vcf \
	> $WORKDIR/VARIANT_CALLING/$2\_Scalpel.fix.vcf

	$BCFTOOLS norm -m -both -D \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$2\_Scalpel.fix.vcf \
 	> $WORKDIR/VARIANT_CALLING/$2\_Scalpel.norm.vcf

	#rm $WORKDIR/VARIANT_CALLING/SCALPEL/variants.indel.vcf
	mv $WORKDIR/VARIANT_CALLING/$2\_Scalpel.fix.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$2\_Scalpel.vcf $DELETE
	
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

	python  $SCRIPT_PIPELINE/header_fix.py -v N -f $WORKDIR/VARIANT_CALLING/$2\_SNVer.sort.vcf \
	> $WORKDIR/VARIANT_CALLING/$2\_SNVer.fix.vcf

	$BCFTOOLS norm -m -both -D \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$2\_SNVer.fix.vcf \
 	> $WORKDIR/VARIANT_CALLING/$2\_SNVer.norm.vcf

 	mv $WORKDIR/VARIANT_CALLING/$2\_SNVer.fix.vcf $DELETE
 	mv $WORKDIR/VARIANT_CALLING/$2\_SNVer.sort.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/SNVER/$2\_SNVer.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/SNVER/$2.filter.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/SNVER/$2.indel.filter.vcf.gz $DELETE
	mv $WORKDIR/VARIANT_CALLING/SNVER/$2.filter.vcf.gz.tbi $DELETE
	mv $WORKDIR/VARIANT_CALLING/SNVER/$2.indel.filter.vcf.gz.tbi $DELETE

	VCF_SNVER=$WORKDIR/VARIANT_CALLING/$2\_SNVer.sort.vcf

	printf $"\n =========>	Sample $2 => Variant Calling: SNVer: DONE\n\n"

}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     PLATYPUS     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Platypus_multisample () {

	printf $'\n =========>	Variant Calling with Platypus\n\n'

	python $PLATYPUS callVariants \
	--refFile=$REF \
	--bamFiles=$1 \
	--regions=$TARGETBED \
	--nCPU=8 \
	--maxSize=5000 \
	--output=$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.dup.vcf

	awk '!a[$0]++' $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.dup.vcf> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.vcf

	python  $SCRIPT_PIPELINE/header_fix.py -v P -f $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.vcf \
	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.fix.vcf

	$BCFTOOLS norm -m -both -D \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.fix.vcf \
 	> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.norm.vcf 

 	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.dup.vcf $DELETE
 	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.fix.vcf $DELETE

	VCF_PLATYPUS=$WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.norm.vcf

 	printf $'\n =========>	Variant Calling with Platypus: DONE\n\n'


}


Platypus_somatic () {

	printf $'\n =========>	Variant Calling with Platypus\n\n'

	python $PLATYPUS callVariants \
	--refFile=$REF \
	--bamFiles=$1 \
	--bamFiles=$2 \
	--regions=$TARGETBED \
	--nCPU=8 \
	--maxSize=5000 \
	--output=$WORKDIR/VARIANT_CALLING/$3\_$4\_Platypus.nosom.vcf

	#awk '!a[$0]++' $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.dup.vcf> $WORKDIR/VARIANT_CALLING/$DATA\_$PANNELLO\_Platypus.vcf

	python ~/NGS_TOOLS/Platypus-master/extensions/Cancer/somaticMutationDetector.py \
	--inputVCF=$WORKDIR/VARIANT_CALLING/$3\_$4\_Platypus.nosom.vcf \
	--outputVCF=$WORKDIR/VARIANT_CALLING/$3\_$4\_Platypus.vcf \
	--tumorSample=$3 \
	--normalSample=$4 \
	--minPosterior=5

	python  $SCRIPT_PIPELINE/header_fix.py -v P -f $WORKDIR/VARIANT_CALLING/$3\_$4\_Platypus.vcf \
	> $WORKDIR/VARIANT_CALLING/$3\_$4\_Platypus.fix.vcf

	$BCFTOOLS norm -m -both -D \
 	-f $REF \
 	$WORKDIR/VARIANT_CALLING/$3\_$4\_Platypus.fix.vcf \
 	> $WORKDIR/VARIANT_CALLING/$3\_$4\_Platypus.norm.vcf

 	mv $WORKDIR/VARIANT_CALLING/$3\_$4\_Platypus.nosom.vcf $DELETE
 	mv $WORKDIR/VARIANT_CALLING/$3\_$4\_Platypus.vcf $DELETE
	mv $WORKDIR/VARIANT_CALLING/$3\_$4\_Platypus.fix.vcf $DELETE

	VCF_PLATYPUS=$WORKDIR/VARIANT_CALLING/$3\_$4\_Platypus.norm.vcf

 	printf $'\n =========>	Variant Calling with Platypus: DONE\n\n'


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
	#rm -f $WORKDIR/gvcf.list
	rm -f $WORKDIR/PostVariantCalling.cfg
	CFG=$WORKDIR/PostVariantCalling.cfg

	cat $1 | while read line
	do

		BAM=$(echo "$line" | cut -f1)
		SAMPLE_NAME=$(echo "$line" | cut -f2)
			
		#HaplotypeCaller $BAM $SAMPLE_NAME

		printf $"$BAM\n" >> $WORKDIR/Bam_list.txt
		printf $"$SAMPLE_NAME\n" >> $WORKDIR/Sample_list.txt

	done

	SAMPLE_NAME=$DATA\_$PANNELLO
	
	#GenotypeGVCFs $WORKDIR/gvcf.list $SAMPLE_NAME
	
	#FreeBayes_multisample $WORKDIR/Bam_list.txt $WORKDIR/Sample_list.txt $SAMPLE_NAME

	#VarScan2_germline_multisample $WORKDIR/Bam_list.txt $WORKDIR/Sample_list.txt $SAMPLE_NAME

	#Platypus_multisample $WORKDIR/Bam_list.txt

	#samtools_vc_multisample $WORKDIR/Bam_list.txt

	#save_in_storage $WORKDIR/Bam_list.txt
	VCF_GATK='/home/jarvis/NGS_ANALYSIS_TEMP/20190215_Run_Prova1_Germline_Cardio_1/VARIANT_CALLING/20190215_Cardio_GATK.norm.vcf'
	VCF_VARSCAN='/home/jarvis/NGS_ANALYSIS_TEMP/20190215_Run_Prova1_Germline_Cardio_1/VARIANT_CALLING/20190215_Cardio_VarScan.norm.vcf'
	VCF_FREEBAYES='/home/jarvis/NGS_ANALYSIS_TEMP/20190215_Run_Prova1_Germline_Cardio_1/VARIANT_CALLING/20190215_Cardio_FreeBayes.norm.vcf'

	#Vcf_merge $VCF_GATK $VCF_FREEBAYES $VCF_VARSCAN

	
	printf $"$VCF_GATK\t$VCF_FREEBAYES\t$VCF_VARSCAN\t$SAMPLE_NAME\n" >> $CFG
	#printf $"$VCF_MERGED\t$SAMPLE_NAME\n" >> $CFG
	
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
		printf $"$VCF_SAMTOOLS\t$VCF_FREEBAYES\t$VCF_VARSCAN\t$SAMPLE_NAME\n" >> $CFG

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
	mkdir -p $WORKDIR/VARIANT_CALLING/NORM
	
	rm -f $WORKDIR/Bam_list.txt
	rm -f $WORKDIR/Sample_list.txt
	rm -f $WORKDIR/gvcf.list
	rm -f $WORKDIR/PostVariantCalling.cfg
	CFG=$WORKDIR/PostVariantCalling.cfg
	

	cat $1 | while read line
	do
		mkdir -p $WORKDIR/VARIANT_CALLING/SCALPEL
		mkdir -p $WORKDIR/VARIANT_CALLING/SNVER
		BAM=$(echo "$line" | cut -f1)
		SAMPLE_NAME=$(echo "$line" | cut -f2)

		printf $"$SAMPLE_NAME" > $WORKDIR/Sample_list.txt

		#samtools_vc $BAM $SAMPLE_NAME
		#HaplotypeCaller $BAM $SAMPLE_NAME
		
		#GenotypeGVCFs $WORKDIR/gvcf.list $SAMPLE_NAME
		

		#FreeBayes_singlesample $BAM $SAMPLE_NAME
		#VarScan2_germline_singlesample $BAM $SAMPLE_NAME

		rm -f $WORKDIR/gvcf.list
		rm -f $WORKDIR/Sample_list.txt

		scalpel_germline_singlesample $BAM $SAMPLE_NAME
		 #rm -rf $WORKDIR/VARIANT_CALLING/SCALPEL

		SNVer_germline_singlesample $BAM $SAMPLE_NAME

		#rm -rf $WORKDIR/VARIANT_CALLING/SNVER
		printf $"$BAM\n" >> $WORKDIR/Bam_list.txt
		#printf $"$VCF_SAMTOOLS\t$VCF_SCALPEL\t$VCF_SNVER\n" >> $CFG
		printf $"$VCF_GATK\t$VCF_FREEBAYES\t$VCF_VARSCAN\t$SAMPLE_NAME\n" >> $CFG

	done
	
	#save_in_storage $WORKDIR/Bam_list.txt
	
	
}
