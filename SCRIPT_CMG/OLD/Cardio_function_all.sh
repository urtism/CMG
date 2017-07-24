#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::     FUNZIONE ALLINEAMENTO     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::     FUNZIONE VARIANT CALLING E ANNOTAZIONE     :::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


CARDIO_VARIANT_CALLING_AND_ANNOTATION () {

#	cd $PROCESSING/5_BQSR/
#	for filename in *.bam
#	do
#
#		cat ~/Scrivania/SCRIPT_PIPELINE/logo_variant.txt
#		printf $"\n =========>	Sample ${filename%.*} => Variant Calling: Haplotype Caller\n\n"
#
#		java -Xmx64g -jar $GATK -T HaplotypeCaller \
#		-R $REF \
#		-I $PROCESSING/5_BQSR/${filename%.*}.bam \
#		-o $PROCESSING/6_Variant/GATK/${filename%.*}.g.vcf \
#		-ERC GVCF \
#		--doNotRunPhysicalPhasing \
#		-bamout $PROCESSING/6_Variant/GATK/${filename%.*}.g.vcf.bam \
#		-L $TARGET
#
#		printf $"\n =========>	Sample ${filename%.*} => Variant Calling: Haplotype Caller: DONE\n\n"
#		cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
#		printf $"\n =========>	Sample ${filename%.*} => Move bam file to Google Genomic folder\n\n"
#
#		mv $PROCESSING/6_Variant/GATK/${filename%.*}.g.vcf.bam $STORAGE/$Name_Dir
#		mv $PROCESSING/6_Variant/GATK/${filename%.*}.g.vcf.bai $STORAGE/$Name_Dir
#		printf "\n${filename%.*}.g.vcf.bam : file bam del .g.vcf GATK per sample ${filename%.*}\n" >> $STORAGE/$Name_Dir/README.txt
#
#		printf $"\n =========>	Sample ${filename%.*} => Move bam file to Google Genomic folder: DONE\n\n"
#		cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
#		printf $"\n =========>	Sample ${filename%.*} => Copy g.vcf file in Google Genomic folder\n\n"
#
#		printf $"\n =========>	Sample ${filename%.*} => Copy g.vcf file in Google Genomic folder: DONE\n\n"
#
#	done




#::::::::::::::::::::::::::::::::::::::::::     INIZIO VARIANT CALLING MULTI-SAMPLE     :::::::::::::::::::::::::::::::::::::::::::::::




	cd $PROCESSING/6_Variant/GATK

	#Stampo tutti i file .g.vcf in una lista
	ls *.g.vcf > samples.list #stampo il nome di tutti  files in una lista .list
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_multi.txt
	echo $'\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling\n\n'

		java -jar -Xmx64g $GATK -T GenotypeGVCFs \
		-R $REF \
		-V:VCF samples.list \
		-o $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK.vcf	

	echo $'\n\n =========>	Variant Calling: Genotype GVCF & Multi-sample variant calling: DONE\n'


#	for file in *.g.vcf
#	do
#		mv $file $STORAGE/$Name_Dir
#		rm $file.idx
#	done



#::::::::::::::::::::::::::::::::::::::::::::::::::::::::     FILTRAGGIO GATK     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



	cat ~/Scrivania/SCRIPT_PIPELINE/logo_filtering.txt
	echo $'\n =========>	Hard Filtering \n\n'

		java -jar $GATK -T VariantFiltration \
		-R $REF \
		-V $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK.vcf \
		--filterExpression "QD < 2.0 || FS > 100.0 || ReadPosRankSum < -16.0 || DP < 15" \
		--filterName "FILTER" \
		-o $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK_Filter.vcf

	echo $'\n =========>	Hard Filtering: DONE\n'
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	GATK vcf manipulation for intersection\n\n'

	rm $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK.vcf

	# Salvo nella cartella storage
	cp $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK_Filter.vcf $STORAGE/$Name_Dir

	python $SCRIPTS_DIR/header_fix.py \
	-f $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK_Filter.vcf \
	-v G \
	> $PROCESSING/6_Variant/GATK/NORM_GATK/$DataRun\_Cardio_GATK_Filter.fixed.vcf

	$BCFTOOLS norm -m -both \
	-f $REF \
	$PROCESSING/6_Variant/GATK/NORM_GATK/$DataRun\_Cardio_GATK_Filter.fixed.vcf \
	> $PROCESSING/6_Variant/GATK/NORM_GATK/$DataRun\_Cardio_GATK_Filter.split.vcf	

	cp $PROCESSING/6_Variant/GATK/NORM_GATK/$DataRun\_Cardio_GATK_Filter.split.vcf $STORAGE/$Name_Dir
	cp $PROCESSING/6_Variant/GATK/NORM_GATK/$DataRun\_Cardio_GATK_Filter.split.vcf $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK_Filter.vcf
	mv $PROCESSING/6_Variant/GATK/NORM_GATK/$DataRun\_Cardio_GATK_Filter.split.vcf $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_GATK_Filter.vcf
	rm $PROCESSING/6_Variant/GATK/$DataRun\_Cardio_GATK.vcf.idx
	bgzip $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_GATK_Filter.vcf
	tabix $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_GATK_Filter.vcf.gz


	echo $'\n =========>	GATK vcf manipulation for intersection: DONE\n'
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	Variant Calling with FreeBayes\n\n'


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     FREEBAYES     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



	cd $PROCESSING/5_BQSR/
	ls *.bam > Bam_list.txt
	ls *.bam > Sample_list.txt
	sed -i -e "s/.bam//g" Sample_list.txt

	/home/jarvis/NGS_TOOLS/freebayes/bin/freebayes -f $REF \
	-L Bam_list.txt \
	-K \
	-J \
	-s Sample_list.txt \
	--genotype-qualities \
	--report-genotype-likelihood-max \
	--allele-balance-priors-off \
	-t $TARGETCARDIOBED > $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes.vcf

	echo $'\n =========>	Variant Calling with FreeBayes: DONE\n'
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	FreeBayes vcf manipulation for intersection\n\n'

	#Rimuovo le righe duplicate:
	cd $PROCESSING/6_Variant/FreeBayes/
	awk '!a[$0]++' $DataRun\_Cardio_FreeBayes.vcf > $DataRun\_Cardio_FreeBayes_RD.vcf
	
	#Converto il vcf alla versione 4.2
	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $DataRun\_Cardio_FreeBayes_RD.vcf

	java -Xmx64g -jar $PICARD SortVcf \
	I=$DataRun\_Cardio_FreeBayes_RD.vcf \
	O=$DataRun\_Cardio_FreeBayes_Sort.vcf

	cp $DataRun\_Cardio_FreeBayes_Sort.vcf $STORAGE/$Name_Dir/$DataRun\_Cardio_FreeBayes.vcf
	printf "\n${DataRun}_Cardio_FreeBayes.vcf : contiene il file vcf multi-sample di FreeBayes in data $DataRun\n" >> $STORAGE/$Name_Dir/README.txt

	python $SCRIPTS_DIR/header_fix.py \
	-f $PROCESSING/6_Variant/FreeBayes/$DataRun\_Cardio_FreeBayes_Sort.vcf \
	-v F \
	> $PROCESSING/6_Variant/FreeBayes/NORM_FREE/$DataRun\_Cardio_FreeBayes.fixed.vcf
		
	$BCFTOOLS norm -m -both \
	-f $REF \
	$PROCESSING/6_Variant/FreeBayes/NORM_FREE/$DataRun\_Cardio_FreeBayes.fixed.vcf \
	> $PROCESSING/6_Variant/FreeBayes/NORM_FREE/$DataRun\_Cardio_FreeBayes.split.vcf	


	#Copio il vcf nella cartella $PROCESSING/6_Variant/Intersect
	cp $PROCESSING/6_Variant/FreeBayes/NORM_FREE/$DataRun\_Cardio_FreeBayes.split.vcf $STORAGE/$Name_Dir
	mv $PROCESSING/6_Variant/FreeBayes/NORM_FREE/$DataRun\_Cardio_FreeBayes.split.vcf $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_FreeBayes_Sort.vcf
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



	cd $PROCESSING/5_BQSR
	samtools mpileup -B -q 1 -f $REF \
	-l $TARGETCARDIOBED -b $PROCESSING/5_BQSR/Bam_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup

	cp Bam_list.txt Sample_list.txt
	sed -i -e 's/.bam//g' Sample_list.txt

	java -jar -Xmx64g $VARSCAN mpileup2snp $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup \
	--min-coverage 10 \
	--min-var-freq 0.20 \
	--pvalue 0.05 \
	--output-vcf 1 \
	--vcf-sample-list Sample_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf

	java -jar -Xmx64g $VARSCAN mpileup2indel $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup \
	--min-coverage 10 \
	--min-var-freq 0.10 \
	--pvalue 0.1 \
	--output-vcf 1 \
	--vcf-sample-list Sample_list.txt > $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf

	echo $'\n =========>	Variant Calling with VarScan2: DONE\n'
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	VarScan2 vcf manipulation for intersection\n\n'
	
	cd $PROCESSING/6_Variant/VarScan/
	#Modifico la versione del vcf
	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $DataRun\_Cardio_VarScan_Indel.vcf
	sed -i -e 's/fileformat=VCFv4.1/fileformat=VCFv4.2/g' $DataRun\_Cardio_VarScan_snp.vcf

	#Copio in Intersect
	cp $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf $PROCESSING/6_Variant/Intersect/
	cp $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf $PROCESSING/6_Variant/Intersect/
	bgzip $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_snp.vcf
	bgzip $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Indel.vcf
	tabix $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_snp.vcf.gz
	tabix $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Indel.vcf.gz

	#Effettuo il merging dei file
	cd $PROCESSING/6_Variant/Intersect/
	vcf-concat $DataRun\_Cardio_VarScan_snp.vcf.gz $DataRun\_Cardio_VarScan_Indel.vcf.gz > $DataRun\_Cardio_VarScan_Merge.vcf
	vcf-sort -c $DataRun\_Cardio_VarScan_Merge.vcf > $DataRun\_Cardio_VarScan_Merge_Sort.vcf


	python $SCRIPTS_DIR/header_fix.py \
	-f $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge_Sort.vcf \
	-v V \
	> $PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$DataRun\_Cardio_VarScan.fixed.vcf
		
	$BCFTOOLS norm -m -both \
	-f $REF \
	$PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$DataRun\_Cardio_VarScan.fixed.vcf \
	> $PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$DataRun\_Cardio_VarScan.split.vcf	



	cp $PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$DataRun\_Cardio_VarScan.split.vcf $STORAGE/$Name_Dir/$DataRun\_Cardio_VarScan.vcf
	mv $PROCESSING/6_Variant/VarScan/NORM_VARSCAN/$DataRun\_Cardio_VarScan.split.vcf $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge_Sort.vcf

	bgzip $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge_Sort.vcf
	tabix $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge_Sort.vcf.gz

	#Rimuovo i file inutili:
	rm $DataRun\_Cardio_VarScan_snp.vcf.gz
	rm $DataRun\_Cardio_VarScan_Indel.vcf.gz
	rm $DataRun\_Cardio_VarScan_snp.vcf.gz.tbi
	rm $DataRun\_Cardio_VarScan_Indel.vcf.gz.tbi
	rm $DataRun\_Cardio_VarScan_Merge.vcf
	rm $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_Indel.vcf
	rm $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio_VarScan_snp.vcf
	rm $PROCESSING/6_Variant/VarScan/$DataRun\_Cardio.mpileup
	
	echo $'\n =========>	VarScan2 vcf manipulation for intersection: DONE\n'
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
	echo $'\n =========>	GATK + FreeBayes + VarScan2 INTERSECTION\n\n'


#------------------------------------------------SALVO I BAM DEL PRE-PROCESSING IN STORAGE---------------------------------------------------


#	cd $PROCESSING/5_BQSR/
#
#	for file in *.bam
#	do
#		mv $file $STORAGE/$Name_Dir
#	done
#
#	for file in *.bai
#	do
#		mv $file $STORAGE/$Name_Dir
#	done



#::::::::::::::::::::::::::::::::::::::::::::::::::::::::     INTERSECTION     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



	vcf-isec -p $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect. \
	$PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_GATK_Filter.vcf.gz \
	$PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge_Sort.vcf.gz \
	$PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_FreeBayes_Sort.vcf.gz	

	tabix $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.1.vcf.gz
	tabix $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.2.vcf.gz
	tabix $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.1_2.vcf.gz

	#Copio le intersezioni in google genomics:
	cp $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.0.vcf.gz $STORAGE/$Name_Dir
	printf "\n${DataRun}_Cardio_Intersect.0.vcf.gz : Intersezione contenente le varianti chiamate da solo GATK\n" >> $STORAGE/$Name_Dir/README.txt
	cp $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.1.vcf.gz $STORAGE/$Name_Dir
	printf "\n${DataRun}_Cardio_Intersect.1.vcf.gz : Intersezione contenente le varianti chiamate da solo VarScan\n" >> $STORAGE/$Name_Dir/README.txt
	cp $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.2.vcf.gz $STORAGE/$Name_Dir
	printf "\n${DataRun}_Cardio_Intersect.2.vcf.gz : Intersezione contenente le varianti chiamate da solo FreeBayes\n" >> $STORAGE/$Name_Dir/README.txt
	cp $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.0_1.vcf.gz $STORAGE/$Name_Dir
	printf "\n${DataRun}_Cardio_Intersect.0_1.vcf.gz : Intersezione contenente le varianti chiamate da GATK + VarScan\n" >> $STORAGE/$Name_Dir/README.txt
	cp $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.0_2.vcf.gz $STORAGE/$Name_Dir
	printf "\n${DataRun}_Cardio_Intersect.0_2.vcf.gz : Intersezione contenente le varianti chiamate GATK + FreeBayes\n" >> $STORAGE/$Name_Dir/README.txt
	cp $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.1_2.vcf.gz $STORAGE/$Name_Dir
	printf "\n${DataRun}_Cardio_Intersect.1_2.vcf.gz : Intersezione contenente le varianti chiamate VarScan + FreeBayes\n" >> $STORAGE/$Name_Dir/README.txt
	cp $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.0_1_2.vcf.gz $STORAGE/$Name_Dir
	printf "\n${DataRun}_Cardio_Intersect.0_1_2.vcf.gz : Intersezione contenente le varianti chiamate da GATK + VarScan + FreeBayes\n" >> $STORAGE/$Name_Dir/README.txt

	cd $PROCESSING/6_Variant/Intersect/

	#Controllo sull'esistenza dei file:
	if [ ! -f $DataRun\_Cardio_Intersect.1.vcf.gz ] && [ ! -f $DataRun\_Cardio_Intersect.2.vcf.gz ] && [ ! -f $DataRun\_Cardio_Intersect.1_2.vcf.gz ]; then
	printf "\nNON ESISTE NESSUN FILE DI INTERSEZIONE!\n"

	elif [ -f $DataRun\_Cardio_Intersect.1.vcf.gz ] && [ -f $DataRun\_Cardio_Intersect.2.vcf.gz ] && [ -f $DataRun\_Cardio_Intersect.1_2.vcf.gz ]; then
	vcf-concat $DataRun\_Cardio_Intersect.1.vcf.gz $DataRun\_Cardio_Intersect.2.vcf.gz \
	$DataRun\_Cardio_Intersect.1_2.vcf.gz > $DataRun\_Cardio_Total_Intersect.vcf
	vcf-sort -c $DataRun\_Cardio_Total_Intersect.vcf > $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf

	elif [ -f $DataRun\_Cardio_Intersect.1.vcf.gz ] && [ -f $DataRun\_Cardio_Intersect.1_2.vcf.gz ] && [ ! -f $DataRun\_Cardio_Intersect.2.vcf.gz ]; then
	vcf-concat $DataRun\_Cardio_Intersect.1.vcf.gz $DataRun\_Cardio_Intersect.1_2.vcf.gz > $DataRun\_Cardio_Total_Intersect.vcf
	vcf-sort -c $DataRun\_Cardio_Total_Intersect.vcf > $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf

	elif [ -f $DataRun\_Cardio_Intersect.2.vcf.gz ] && [ -f $DataRun\_Cardio_Intersect.1_2.vcf.gz ] && [ ! -f $DataRun\_Cardio_Intersect.1.vcf.gz ]; then
	vcf-concat $DataRun\_Cardio_Intersect.2.vcf.gz $DataRun\_Cardio_Intersect.1_2.vcf.gz > $DataRun\_Cardio_Total_Intersect.vcf
	vcf-sort -c $DataRun\_Cardio_Total_Intersect.vcf > $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf

	elif [ -f $DataRun\_Cardio_Intersect.1.vcf.gz ] && [ -f $DataRun\_Cardio_Intersect.2.vcf.gz ] && [ ! -f $DataRun\_Cardio_Intersect.1_2.vcf.gz ]; then
	vcf-concat $DataRun\_Cardio_Intersect.1.vcf.gz $DataRun\_Cardio_Intersect.2.vcf.gz > $DataRun\_Cardio_Total_Intersect.vcf
	vcf-sort -c $DataRun\_Cardio_Total_Intersect.vcf > $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf

	elif [ -f $DataRun\_Cardio_Intersect.1.vcf.gz ] && [ ! -f $DataRun\_Cardio_Intersect.2.vcf.gz ] && [ ! -f $DataRun\_Cardio_Intersect.1_2.vcf.gz ]; then
	gunzip $DataRun\_Cardio_Intersect.1.vcf.gz
	cp $DataRun\_Cardio_Intersect.1.vcf > $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf
	rm $DataRun\_Cardio_Intersect.1.vcf

	elif [ -f $DataRun\_Cardio_Intersect.1_2.vcf.gz ] && [ ! -f $DataRun\_Cardio_Intersect.1.vcf.gz ] && [ ! -f $DataRun\_Cardio_Intersect.2.vcf.gz ]; then
	gunzip $DataRun\_Cardio_Intersect.1_2.vcf.gz
	cp $DataRun\_Cardio_Intersect.1_2.vcf > $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf
	rm $DataRun\_Cardio_Intersect.1_2.vcf

	elif [ -f $DataRun\_Cardio_Intersect.2.vcf.gz ] && [ ! -f $DataRun\_Cardio_Intersect.1_2.vcf.gz ] && [ ! -f $DataRun\_Cardio_Intersect.1.vcf.gz ]; then
	gunzip $DataRun\_Cardio_Intersect.1.vcf.gz
	cp $DataRun\_Cardio_Intersect.2.vcf > $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf
	rm $DataRun\_Cardio_Intersect.2.vcf

	fi

	rm $DataRun\_Cardio_Total_Intersect.vcf

	# Salvo nella cartella storage
	cp $PROCESSING/7_Filter/$DataRun\_Cardio_Total_Intersect_Sort.vcf $STORAGE/$Name_Dir/$DataRun\_Cardio_Intersect.vcf
	printf "\n${DataRun}_Cardio_Intersect.vcf : contiene il vcf con le varianti chiamate da VarScan e/o FreeBayes ma non in GATK in data $DataRun\n" >> $STORAGE/$Name_Dir/README.txt

	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.0.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.1.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.1.vcf.gz.tbi
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.2.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.2.vcf.gz.tbi
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.0_1.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.0_2.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.1_2.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.1_2.vcf.gz.tbi
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect.0_1_2.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_FreeBayes_Sort.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_FreeBayes_Sort.vcf.gz.tbi
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_GATK_Filter.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_GATK_Filter.vcf.gz.tbi
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_Intersect._README
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge_Sort.vcf.gz
	rm $PROCESSING/6_Variant/Intersect/$DataRun\_Cardio_VarScan_Merge_Sort.vcf.gz.tbi



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



#:::::::::::::::::::::::::::::::::::::::::::::::::::    VCF MODIFICATION FOR EXCEL     ::::::::::::::::::::::::::::::::::::::::::::::::::::::



#Estraggo la HEADER a cui poi andrò a concatenare le varianti
sed -n -e '/#/p' $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_ANN.vcf > \
$PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN_SPLIT.vcf
sed -n -e '/#/p' $PROCESSING/8_Annotation/$DataRun\_Cardio_Total_Intersect_ANN.vcf > \
$PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_Total_Intersect_ANN_SPLIT.vcf

#Splitto le annotazioni per trascritti
perl -e'while(<>) { chomp; if(m/(.+?)(CSQ|ANN)\=([^;^\s]+)(.*)/) { foreach my $s(split ",", $3) { print "$1$2\=$s\;$4\n"}}}' \
$PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_ANN.vcf >> $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN_SPLIT.vcf
perl -e'while(<>) { chomp; if(m/(.+?)(CSQ|ANN)\=([^;^\s]+)(.*)/) { foreach my $s(split ",", $3) { print "$1$2\=$s\;$4\n"}}}' \
$PROCESSING/8_Annotation/$DataRun\_Cardio_Total_Intersect_ANN.vcf >> \
$PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_Total_Intersect_ANN_SPLIT.vcf

#Copio i file splittati nella cartella $VEP ed in quella dello STORAGE
cp $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN_SPLIT.vcf $STORAGE/$Name_Dir
printf "\n${DataRun}_Cardio_GATK_ANN_SPLIT.vcf : vcf GATK splittato per trascritti\n" >> $STORAGE/$Name_Dir/README.txt
cp $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_Total_Intersect_ANN_SPLIT.vcf $STORAGE/$Name_Dir/$DataRun\_Cardio_Intersect_ANN_SPLIT.vcf 
printf "\n${DataRun}_Cardio_Intersect_ANN_SPLIT.vcf : vcf Intersect splittato per trascritti\n" >> $STORAGE/$Name_Dir/README.txt
mv $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_ANN.vcf $PROCESSING/8_Annotation/FILE_PROCESSING/
rm $PROCESSING/8_Annotation/$DataRun\_Cardio_Total_Intersect_ANN.vcf

#Eseguo il filtraggio sui trascritti

cd $PROCESSING/8_Annotation/FILE_PROCESSING/

grep -e "#" -f ~/NGS_ANALYSIS/TARGET/Lista_trascritti_cardio.txt $DataRun\_Cardio_GATK_ANN_SPLIT.vcf > $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_Filter_Transcripts.vcf

grep -e "#" -f ~/NGS_ANALYSIS/TARGET/Lista_trascritti_cardio.txt $DataRun\_Cardio_Total_Intersect_ANN_SPLIT.vcf > $PROCESSING/8_Annotation/$DataRun\_Cardio_Intersect_Filter_Transcripts.vcf

#Copio il file gatk filtrato sui trascritti in file processing per fare il bgzip e poi intersecare i file!
cp $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_Filter_Transcripts.vcf $PROCESSING/8_Annotation/FILE_PROCESSING

#Aggiungi la copia in google genomics
cp $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_Filter_Transcripts.vcf $STORAGE/$Name_Dir
printf "\n${DataRun}_Cardio_GATK_Filter_Transcripts.vcf : vcf GATK filtrato per trascritti\n" >> $STORAGE/$Name_Dir/README.txt
cp $PROCESSING/8_Annotation/$DataRun\_Cardio_Intersect_Filter_Transcripts.vcf $STORAGE/$Name_Dir
printf "\n${DataRun}_Cardio_Intersect_Filter_Transcripts.vcf : vcf Intersect filtrato per trascritti\n" >> $STORAGE/$Name_Dir/README.txt
rm $DataRun\_Cardio_GATK_ANN_SPLIT.vcf
rm $DataRun\_Cardio_Total_Intersect_ANN_SPLIT.vcf



#::::::::::::::::::::::::::::::::::::::::::    VCF MODIFICATION FOR EXCEL OTHER TRANSCRIPTS     :::::::::::::::::::::::::::::::::::::::::::::



#Ripeto le operazioni di prima intersecando i file per avere gli altri trascritti. Prima interseco i file e poi li splitto.

bgzip $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN.vcf
tabix $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN.vcf.gz
bgzip $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_Filter_Transcripts.vcf
tabix $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_Filter_Transcripts.vcf.gz

vcf-isec -p $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN_Intersect. $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN.vcf.gz $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_Filter_Transcripts.vcf.gz

rm $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN.vcf.gz
rm $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN.vcf.gz.tbi
rm $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_Filter_Transcripts.vcf.gz
rm $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_Filter_Transcripts.vcf.gz.tbi
rm $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN_Intersect.0_1.vcf.gz
rm $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN_Intersect._README

gunzip $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN_Intersect.0.vcf.gz

sed -n -e '/#/p' $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN_Intersect.0.vcf > \
$PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_Other_Transcripts.vcf

perl -e'while(<>) { chomp; if(m/(.+?)(CSQ|ANN)\=([^;^\s]+)(.*)/) { foreach my $s(split ",", $3) { print "$1$2\=$s\;$4\n"}}}' \
$PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN_Intersect.0.vcf >> $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_Other_Transcripts.vcf

mv $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_ANN_Intersect.0.vcf  $STORAGE/$Name_Dir/$DataRun\_Cardio_GATK_ANN_Other_Transcripts.vcf
printf "\n${DataRun}\_Cardio_GATK_ANN_Other_Transcripts.vcf : vcf GATK annotato per altri trascritti non in studio\n" >> $STORAGE/$Name_Dir/README.txt

cp $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_Other_Transcripts.vcf $STORAGE/$Name_Dir/$DataRun\_Cardio_GATK_SPLIT_Other_Transcripts.vcf
printf "\n${DataRun}\_Cardio_GATK_SPLIT_Other_Transcripts.vcf : vcf GATK splittato per altri trascritti non in studio\n" >> $STORAGE/$Name_Dir/README.txt

mv $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_Other_Transcripts.vcf $PROCESSING/8_Annotation



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::   VCF MODIFICATION     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



sed -n -e '/#CHROM/p' $PROCESSING/8_Annotation/$DataRun\_Cardio_GATK_Filter_Transcripts.vcf > $PROCESSING/8_Annotation/FILE_PROCESSING/List_samples_for_split.txt
sed -i "s/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t//g" $PROCESSING/8_Annotation/FILE_PROCESSING/List_samples_for_split.txt
COUNT=$(awk '{print NF}' $PROCESSING/8_Annotation/FILE_PROCESSING/List_samples_for_split.txt | sort -nu | tail -n 1)

egrep "CIGAR|##|#" $PROCESSING/8_Annotation/$DataRun\_Cardio_Intersect_Filter_Transcripts.vcf \
> $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_Intersect_Transcripts_FREEBAYES.vcf
egrep "FREQ|##|#" $PROCESSING/8_Annotation/$DataRun\_Cardio_Intersect_Filter_Transcripts.vcf \
> $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_Intersect_Transcripts_VarScan.vcf
rm $PROCESSING/8_Annotation/$DataRun\_Cardio_Intersect_Filter_Transcripts.vcf

#Copio in Genomics
cp $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_Intersect_Transcripts_FREEBAYES.vcf $STORAGE/$Name_Dir/$DataRun\_Cardio_Transcripts_FREEBAYES.vcf
printf "\n${DataRun}_Cardio_Transcripts_FREEBAYES.vcf : vcf di sole varianti FreeBayes filtrato per trascritti\n" >> $STORAGE/$Name_Dir/README.txt
cp $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_Intersect_Transcripts_VarScan.vcf $STORAGE/$Name_Dir/$DataRun\_Cardio_Transcripts_VarScan.vcf
printf "\n${DataRun}_Cardio_Transcripts_VarScan.vcf : vcf di sole varianti VarScan filtrato per trascritti\n" >> $STORAGE/$Name_Dir/README.txt

cd $PROCESSING/8_Annotation/
#Modifico i file per splittare i campi in seguito. Aggiungo i ;
for vcf in *GATK_Filter_Transcripts.vcf
	do
	sed -i -e "s/\(AN=[[:digit:]]*\);\(DP=[[:digit:]]*\)/\1;;;\2/g" -e "s/FS=\([[:digit:]]*\|[[:digit:]]*.[[:digit:]]*\);MLEAC=\([[:digit:]]*\|[[:digit:]]*.[[:digit:]]\)/FS=\1;;MLEAC=\2/g" -e "s/MQ=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\);QD=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\)/MQ=\1;;QD=\2/g" -e "s/QD=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\);SOR=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\)/QD=\1;;SOR=\2/g" $PROCESSING/8_Annotation/$vcf
mv $PROCESSING/8_Annotation/$vcf $PROCESSING/8_Annotation/FILE_PROCESSING/
done

for vcf in *Other_Transcripts.vcf
	do
	sed -i -e "s/\(AN=[[:digit:]]*\);\(DP=[[:digit:]]*\)/\1;;;\2/g" -e "s/FS=\([[:digit:]]*\|[[:digit:]]*.[[:digit:]]*\);MLEAC=\([[:digit:]]*\|[[:digit:]]*.[[:digit:]]\)/FS=\1;;MLEAC=\2/g" -e "s/MQ=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\);QD=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\)/MQ=\1;;QD=\2/g" -e "s/QD=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\);SOR=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\)/QD=\1;;SOR=\2/g" $PROCESSING/8_Annotation/$vcf
mv $PROCESSING/8_Annotation/$vcf $PROCESSING/8_Annotation/FILE_PROCESSING/
done	



#---------------------------------------------------SPLITTO I VCF DI GATK IN FORMATO CSV-----------------------------------------------------



cd $PROCESSING/8_Annotation/FILE_PROCESSING/
for file in *GATK_Filter_Transcripts.vcf
do
	for (( a=1; a<$COUNT+1; a++ ))
	do
	NAME=$(cut -f $a List_samples_for_split.txt)
	b=$((9+$a))
	cut -f1-9,$b $file > $PROCESSING/8_Annotation/GATK/$NAME\_GATK.vcf
	sed -i -e '/0\/0/d' -e '/\.\/\./d' -e 's/|;/|-;/g' -e 's/;\t/\t/g' $PROCESSING/8_Annotation/GATK/$NAME\_GATK.vcf
	cp $PROCESSING/8_Annotation/GATK/$NAME\_GATK.vcf $STORAGE/$Name_Dir
	printf "\n${NAME}_GATK.vcf : vcf GATK splittato e annotato per il sample $NAME\n" >> $STORAGE/$Name_Dir/README.txt
	sed -i -e '/^##/d' $PROCESSING/8_Annotation/GATK/$NAME\_GATK.vcf
	# Aggiungo le intestazioni delle annotazioni:
	#Splitto in corrispondenza dei ; e di | e dei : poi aggiungo - in corrispondenza dei valori mancanti e rimuovo le stringhe in INFO, aggiungo il 	nome paziente nella colonna ID ecc...
	sed -i -e "s/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$NAME/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAC\tAF\tAN\tBaseQRankSum\tClippingRankSum\tDP\tExcessHet\tFS\tInbreedingCoeff\tMLEAC\tMLEAF\tMQ\tMQRankSum\tQD\tReadPosRankSum\tSOR\tAllele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tVARIANT_CLASS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tENSP\tSIFT\tPolyPhen\tHGVS_OFFSET\tGMAF\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tGT\tAD\tDP\tGQ\tPL/g" -e "s/;\tGT/\tGT/g" -e "s/GT:AD:DP:GQ:PL\t//g" $PROCESSING/8_Annotation/GATK/$NAME\_GATK.vcf
	sed -i -e "s/\(^chr[[:alnum:]]*\t[[:digit:]]*\t\)\.\t/\1$NAME\t/g" -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/:/\t/g' -e 's/\[\]/:/g' -e 's/;/\t/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/|/\t/g' -e "s/AC=\|AF=\|AN=\|BaseQRankSum=\|ClippingRankSum=\|DP=\|ExcessHet=\|FS=\|InbreedingCoeff=\|MLEAC=\|MLEAF=\|MQ=\|MQRankSum=\|QD=\|ReadPosRankSum=\|SOR=\|ANN=//g" -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' $PROCESSING/8_Annotation/GATK/$NAME\_GATK.vcf
	cut --complement -f26,28,29,31,42-44,46,47,52,55,56 $PROCESSING/8_Annotation/GATK/$NAME\_GATK.vcf > $STORAGE/$Name_Dir/$NAME\_GATK_Conferme.tsv
	mv $PROCESSING/8_Annotation/GATK/$NAME\_GATK.vcf $STORAGE/$Name_Dir/$NAME\_GATK_All_Tags.txt
	printf "\n${NAME}_GATK_All_Tags.txt : file txt GATK contenente sole varianti per il sample $NAME ma con tutte le tag del vcf\n" >> $STORAGE/$Name_Dir/README.txt
	printf "\n${NAME}_GATK_Conferme.tsv : file .tsv GATK definitivo contenente sole varianti per il sample $NAME\n" >> $STORAGE/$Name_Dir/README.txt
	cp $STORAGE/$Name_Dir/$NAME\_GATK_Conferme.tsv $OUTVCF/$Name_OUT
	done
done

rm $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_Filter_Transcripts.vcf



#--------------------------------------------------RIPETO OPERAZIONE SU FILE OTHER TRANSCRIPTS-----------------------------------------------



cd $PROCESSING/8_Annotation/FILE_PROCESSING/
for file in *Other_Transcripts.vcf
do
	for (( a=1; a<$COUNT+1; a++ ))
	do
	NAME=$(cut -f $a List_samples_for_split.txt)
	b=$((9+$a))
	cut -f1-9,$b $file > $PROCESSING/8_Annotation/GATK/$NAME\_GATK_Other_Transcripts.vcf
	sed -i -e '/0\/0/d' -e '/\.\/\./d' -e 's/|;/|-;/g' -e 's/;\t/\t/g' $PROCESSING/8_Annotation/GATK/$NAME\_GATK_Other_Transcripts.vcf
	cp $PROCESSING/8_Annotation/GATK/$NAME\_GATK_Other_Transcripts.vcf $STORAGE/$Name_Dir
	printf "\n${NAME}_GATK_Other_Transcripts.vcf : vcf GATK splittato e annotato per il sample $NAME contenente altri trascritti\n" >> $STORAGE/$Name_Dir/README.txt
	sed -i -e '/^##/d' $PROCESSING/8_Annotation/GATK/$NAME\_GATK_Other_Transcripts.vcf
	# Aggiungo le intestazioni delle annotazioni:
	#Splitto in corrispondenza dei ; e di | e dei : poi aggiungo - in corrispondenza dei valori mancanti e rimuovo le stringhe in INFO, aggiungo il 	nome paziente nella colonna ID ecc...
	sed -i -e "s/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$NAME/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAC\tAF\tAN\tBaseQRankSum\tClippingRankSum\tDP\tExcessHet\tFS\tInbreedingCoeff\tMLEAC\tMLEAF\tMQ\tMQRankSum\tQD\tReadPosRankSum\tSOR\tAllele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tVARIANT_CLASS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tENSP\tSIFT\tPolyPhen\tHGVS_OFFSET\tGMAF\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tGT\tAD\tDP\tGQ\tPL/g" -e "s/;\tGT/\tGT/g" -e "s/GT:AD:DP:GQ:PL\t//g" $PROCESSING/8_Annotation/GATK/$NAME\_GATK_Other_Transcripts.vcf
	sed -i -e "s/\(^chr[[:alnum:]]*\t[[:digit:]]*\t\)\.\t/\1$NAME\t/g" -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/:/\t/g' -e 's/\[\]/:/g' -e 's/;/\t/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/|/\t/g' -e "s/AC=\|AF=\|AN=\|BaseQRankSum=\|ClippingRankSum=\|DP=\|ExcessHet=\|FS=\|InbreedingCoeff=\|MLEAC=\|MLEAF=\|MQ=\|MQRankSum=\|QD=\|ReadPosRankSum=\|SOR=\|ANN=//g" -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' $PROCESSING/8_Annotation/GATK/$NAME\_GATK_Other_Transcripts.vcf
	cut --complement -f26,28,29,31,42-44,46,47,52,55,56 $PROCESSING/8_Annotation/GATK/$NAME\_GATK_Other_Transcripts.vcf > $STORAGE/$Name_Dir/$NAME\_GATK_Other_Transcripts_Conferme.tsv
	mv $PROCESSING/8_Annotation/GATK/$NAME\_GATK_Other_Transcripts.vcf $STORAGE/$Name_Dir/$NAME\_GATK_Other_Transcripts_All_Tags.txt
	printf "\n${NAME}_GATK_Other_Transcripts_All_Tags.txt : file txt GATK contenente sole varianti per il sample $NAME del file other transcripts ma con tutte le tag del vcf\n" >> $STORAGE/$Name_Dir/README.txt
	printf "\n${NAME}_GATK_Other_Transcripts_Conferme.tsv : file .tsv GATK definitivo contenente sole varianti per il sample $NAME\n" >> $STORAGE/$Name_Dir/README.txt
	cp $STORAGE/$Name_Dir/$NAME\_GATK_Other_Transcripts_Conferme.tsv $OUTVCF/$Name_OUT
	done
done

rm $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_GATK_Other_Transcripts.vcf



#-------------------------------------------------SPLITTO I VCF DI VARSCAN IN FORMATO CSV----------------------------------------------------



cd $PROCESSING/8_Annotation/FILE_PROCESSING
for file in *\_VarScan.vcf
do
	for (( a=1; a<$COUNT+1; a++ ))
	do
	NAME=$(cut -f $a List_samples_for_split.txt)
	b=$((9+$a))
	cut -f1-9,$b $file > $PROCESSING/8_Annotation/VARSCAN/$NAME\_VarScan.vcf
	sed -i -e '/0\/0/d' -e '/\.\/\./d' -e 's/|;/|-;/g' -e 's/;\t/\t/g' $PROCESSING/8_Annotation/VARSCAN/$NAME\_VarScan.vcf
	cp $PROCESSING/8_Annotation/VARSCAN/$NAME\_VarScan.vcf $STORAGE/$Name_Dir
	printf "\n${NAME}_VarScan.vcf : vcf VarScan splittato e annotato per il sample $NAME\n" >> $STORAGE/$Name_Dir/README.txt
	sed -i -e '/^##/d'  $PROCESSING/8_Annotation/VARSCAN/$NAME\_VarScan.vcf
	# Aggiungo le intestazioni delle annotazioni:
	#Splitto in corrispondenza dei ; e di | e dei : poi aggiungo - in corrispondenza dei valori mancanti e rimuovo le stringhe in INFO, aggiungo il 	nome paziente nella colonna ID ecc...
	sed -i -e "s/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$NAME/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tADP\tWT\tHET\tHOM\tNC\tAllele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tVARIANT_CLASS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tENSP\tSIFT\tPolyPhen\tHGVS_OFFSET\tGMAF\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tGT\tGQ\tSDP\tDP\tRD\tAD\tFREQ\tPVAL\tRBQ\tABQ\tRDF\tRDR\tADF\tADR/g" -e "s/;\tGT/\tGT/g" -e "s/GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t//g" $PROCESSING/8_Annotation/VARSCAN/$NAME\_VarScan.vcf
	sed -i -e "s/GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t//g" -e "s/\(^chr[[:alnum:]]*\t[[:digit:]]*\t\)\.\t/\1$NAME\t/g" -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/:/\t/g' -e 's/\[\]/:/g' -e 's/;/\t/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/|/\t/g' -e "s/ADP=\|WT=\|HET=\|HOM=\|NC=\|ANN=//g" -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' $PROCESSING/8_Annotation/VARSCAN/$NAME\_VarScan.vcf
	cut --complement -f15,17,18,20,31-33,35-36,41,44,45 $PROCESSING/8_Annotation/VARSCAN/$NAME\_VarScan.vcf > $STORAGE/$Name_Dir/$NAME\_VarScan_Conferme.tsv
	mv $PROCESSING/8_Annotation/VARSCAN/$NAME\_VarScan.vcf $STORAGE/$Name_Dir/$NAME\_VarScan_All_Tags.txt
	printf "\n${NAME}_VarScan_All_Tags.txt : file txt VarScan contenente sole varianti per il sample $NAME ma con tutte le tag del vcf\n" >> $STORAGE/$Name_Dir/README.txt
	printf "\n${NAME}_VarScan_Conferme.tsv : file .tsv VarScan definitivo contenente sole varianti per il sample $NAME\n" >> $STORAGE/$Name_Dir/README.txt
	cp $STORAGE/$Name_Dir/$NAME\_VarScan_Conferme.tsv $OUTVCF/$Name_OUT
	done
done

rm $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_Intersect_Transcripts_VarScan.vcf



#-----------------------------------------------SPLITTO I VCF DI FREEBAYES IN FORMATO CSV----------------------------------------------------



cd $PROCESSING/8_Annotation/FILE_PROCESSING
for file in *\_FREEBAYES.vcf
do
	for (( a=1; a<$COUNT+1; a++ ))
	do
	NAME=$(cut -f $a List_samples_for_split.txt)
	b=$((9+$a))
	cut -f1-9,$b $file > $PROCESSING/8_Annotation/FREEBAYES/$NAME\_FREEBAYES.vcf
	#Elimino i non varianti
	sed -i -e '/0\/0/d' -e '/\.\/\./d' -e '/GL\t\./d' -e 's/|;/|-;/g' -e 's/;\t/\t/g' $PROCESSING/8_Annotation/FREEBAYES/$NAME\_FREEBAYES.vcf
	cp $PROCESSING/8_Annotation/FREEBAYES/$NAME\_FREEBAYES.vcf $STORAGE/$Name_Dir
	printf "\n${NAME}_FREEBAYES.vcf : vcf FREEBAYES splittato e annotato per il sample $NAME\n" >> $STORAGE/$Name_Dir/README.txt
	# Elimino le righe con ##
	sed -i -e '/^##/d' $PROCESSING/8_Annotation/FREEBAYES/$NAME\_FREEBAYES.vcf
	# Aggiungo le intestazioni delle annotazioni:
	#Splitto in corrispondenza dei ; e di | e dei : poi aggiungo - in corrispondenza dei valori mancanti e rimuovo le stringhe in INFO, aggiungo il nome paziente nella colonna ID ecc...
	sed -i -e "s/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$NAME/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAB\tABP\tAC\tAF\tAN\tAO\tCIGAR\tDP\tDPB\tDPRA\tEPP\tEPPR\tGTI\tLEN\tMEANALT\tMQM\tMQMR\tNS\tNUMALT\tODDS\tPAIRED\tPAIREDR\tPAO\tPQA\tPQR\tPRO\tQA\tQR\tRO\tRPL\tRPP\tRPPR\tRPR\tRUN\tSAF\tSAP\tSAR\tSRF\tSRP\tSRR\tTYPE\ttechnology.ILLUMINA\tAllele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tVARIANT_CLASS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tENSP\tSIFT\tPolyPhen\tHGVS_OFFSET\tGMAF\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tGT\tGQ\tDP\tAD\tRO\tQR\tAO\tQA\tGL/g" -e "s/GT:GQ:DP:AD:RO:QR:AO:QA:GL\t//g" $PROCESSING/8_Annotation/FREEBAYES/$NAME\_FREEBAYES.vcf	
	sed -i -e "s/\(^chr[[:alnum:]]*\t[[:digit:]]*\t\)\.\t/\1$NAME\t/g" -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/:/\t/g' -e 's/\[\]/:/g' -e 's/;/\t/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/|/\t/g' -e "s/AB=\|ABP=\|AC=\|AF=\|AN=\|AO=\|CIGAR=\|DP=\|DPB=\|DPRA=\|EPP=\|EPPR=\|GTI=\|LEN=\|MEANALT=\|MQM=\|MQMR=\|NS=\|NUMALT=\|ODDS=\|PAIRED=\|PAIREDR=\|PAO=\|PQA=\|PQR=\|PRO=\|QA=\|QR=\|RO=\|RPL=\|RPP=\|RPPR=\|RPR=\|RUN=\|SAF=\|SAP=\|SAR=\|SRF=\|SRP=\|SRR=\|TYPE=\|technology.ILLUMINA=\|ANN=//g" -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' $PROCESSING/8_Annotation/FREEBAYES/$NAME\_FREEBAYES.vcf
	cut --complement -f9,14,17-22,27-33,37-41,43,46,49,52,54,55,57,68-70,72,73,78,81,82 $PROCESSING/8_Annotation/FREEBAYES/$NAME\_FREEBAYES.vcf > $STORAGE/$Name_Dir/$NAME\_FREEBAYES_Conferme.tsv
	mv $PROCESSING/8_Annotation/FREEBAYES/$NAME\_FREEBAYES.vcf $STORAGE/$Name_Dir/$NAME\_FREEBAYES_All_Tags.txt
	printf "\n${NAME}_FREEBAYES_All_Tags.txt : file txt FREEBAYES contenente sole varianti per il sample $NAME ma con tutte le tag del vcf\n" >> $STORAGE/$Name_Dir/README.txt
	printf "\n${NAME}_FREEBAYES_Conferme.tsv : file .tsv FREEBAYES definitivo contenente sole varianti per il sample $NAME\n" >> $STORAGE/$Name_Dir/README.txt
	cp $STORAGE/$Name_Dir/$NAME\_FREEBAYES_Conferme.tsv $OUTVCF/$Name_OUT
	done
done


	echo $'\n =========>	FILE MODIFICATION FOR CSV FORMAT: DONE \n\n'

cp $STORAGE/$Name_Dir/$DataRun\_Cardio_GATK_Filter.vcf $OUTVCF/$Name_OUT
cp $STORAGE/$Name_Dir/$DataRun\_Cardio_VarScan.vcf $OUTVCF/$Name_OUT
cp $STORAGE/$Name_Dir/$DataRun\_Cardio_FreeBayes.vcf $OUTVCF/$Name_OUT
cp $STORAGE/$Name_Dir/$DataRun\_Cardio_Intersect.vcf $OUTVCF/$Name_OUT


rm $PROCESSING/8_Annotation/FILE_PROCESSING/$DataRun\_Cardio_Intersect_Transcripts_FREEBAYES.vcf

#Rimuovo i file rimanenti inutili
rm $PROCESSING/5_BQSR/Bam_list.txt
rm $PROCESSING/5_BQSR/Sample_list.txt
rm $PROCESSING/6_Variant/GATK/samples.list
rm $PROCESSING/8_Annotation/FILE_PROCESSING/List_samples_for_split.txt

}
