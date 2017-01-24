PREPROCESSING () {

	cat ~/Scrivania/SCRIPT_PIPELINE/logo_processing.txt

	mkdir $WORKDIR/PREPROCESSING

	for bamsort in $WORKDIR/Alignment/*.sort.bam
	do

		sample_name="$(basename "${bamsort%.*}")"
		alig_path="$(dirname "${bamsort%.*}")"



		printf $"\n~~~>	Sample ${*.bamsort%.*.*} => Add Or Replace Read Groups\n\n"

		java -Xmx64g -jar $PICARD AddOrReplaceReadGroups \
		I=$bamsort \
		O=${bamsort%.*}.Add.bam \
		RGID=$SAMPLENAME \
		RGPL=ILLUMINA \
		RGSM=$SAMPLENAME \
		RGPU=ILLUMINA_$DATA \
		RGLB=$PANELLO \
		VALIDATION_STRINGENCY=LENIENT

		if [ "$DESIGN" == "ENRICHMENT" ]
		then
			printf $"\n~~~>	Sample $sample_name => Add Or Replace Read Groups: DONE\n\n"
			printf $"\n~~~>	Sample $sample_name => Mark Duplicates\n\n"

			java -Xmx64g -jar $PICARD MarkDuplicates \
			I=${bamsort%.*}.Add.bam \
			O=${bamsort%.*}.mark.bam \
			METRICS_FILE=${bamsort%.*}.Metrics.txt \
			READ_NAME_REGEX=null \
			VALIDATION_STRINGENCY=LENIENT \
			REMOVE_DUPLICATES=true ASSUME_SORTED=true

			printf $"\n~~~>	Sample $sample_name => Mark Duplicates: DONE\n\n"
			printf $"\n~~~>	Sample $sample_name => Build Bam Index\n\n"

			java -Xmx64g -jar $PICARD BuildBamIndex \
			I=${bamsort%.*}.mark.bam \
			O=${bamsort%.*}.mark.bai \
			VALIDATION_STRINGENCY=LENIENT

			printf $"\n~~~>	Sample $sample_name => Build Bam Index: DONE\n\n"
			printf $"\n~~~>	Sample $sample_name => Realigner Target Creator\n\n"

			java -Xmx64g -jar $GATK -T RealignerTargetCreator \
			-R $REF \
			-I ${bamsort%.*}.mark.bam \
			-o ${bamsort%.*}.IndelRealigner.intervals \
			--known $MILLS \
			-L $TARGET

			printf $"\n~~~>	Sample $sample_name => Realigner Target Creator: DONE\n\n"
			printf $"\n~~~>	Sample $sample_name => Indel Realigner\n\n"

			java -Xmx64g -jar $GATK -T IndelRealigner \
			-R $REF \
			-I ${bamsort%.*}.debup.bam \
			-targetIntervals ${bamsort%.*}.IndelRealigner.intervals \
			-o ${bamsort%.*}.Realigned.bam \
			-known $MILLS

			echo $'\n =========>	Indel Realigner: DONE\n'
			printf $"\n~~~>	Sample $sample_name => Base Quality Score Recalibration (BQSR): Base Recalibrator\n\n"

			java -Xmx64g -jar $GATK -T BaseRecalibrator \
			-R $REF \
			-I ${bamsort%.*}.Realigned.bam \
			-knownSites $DBSNP \
			-knownSites $MILLS \
			-o ${bamsort%.*}.Recal.table \
			-L $TARGET

			printf $"\n~~~>	Sample $sample_name => Base Quality Score Recalibration (BQSR): Base Recalibrator: DONE\n\n"
			printf $"\n~~~>	Sample $sample_name => Base Quality Score Recalibration (BQSR): Print reads\n\n"

			java -Xmx64g -jar $GATK -T PrintReads \
			-R $REF \
			-I ${bamsort%.*}.Realigned.bam \
			-BQSR ${bamsort%.*}.Recal.table  \
			-o $WORKDIR/PREPROCESSING/$sample_name.bam \
			-L $TARGET

			rm $bamsort
			rm ${bamsort%.*}.Add.bam
			rm ${bamsort%.*}.Metrics.txt
			rm ${bamsort%.*}.mark.bam
			rm ${bamsort%.*}.mark.bai
			rm ${bamsort%.*}.Recal.table
			rm ${bamsort%.*}.Realigned.bam
			rm ${bamsort%.*}.Realigned.bai
			rm ${bamsort%.*}.IndelRealigner.intervals

			printf $"\n~~~>	Sample $sample_name => Base Quality Score Recalibration (BQSR): Print reads: DONE\n\n"

		elif [ "$DESIGN" == "AMPLICON" ]
		then
			rm ${bamsort%.*}.Add.bam $WORKDIR/PREPROCESSING/$sample_name.bam

	done

}