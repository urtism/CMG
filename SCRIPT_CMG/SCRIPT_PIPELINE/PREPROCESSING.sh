
AddOrReplaceReadGroups () {

	printf $"\n~~~>	Sample $SAMPLE_NAME => Add Or Replace Read Groups\n\n"

	java -Xmx64g -jar $PICARD AddOrReplaceReadGroups \
	I=$1 \
	O=${1%.*.*}.Add.bam \
	RGID=$SAMPLE_NAME \
	RGPL=ILLUMINA \
	RGSM=$SAMPLE_NAME \
	RGPU=ILLUMINA_$DATA \
	RGLB=$PANNELLO \
	VALIDATION_STRINGENCY=LENIENT

	mv $1 $DELETE
	mv ${1%.*}.bai $DELETE
	INPUT=${1%.*.*}.Add.bam
	printf $"\n~~~>	Sample $SAMPLE_NAME => Add Or Replace Read Groups: DONE\n\n"
}

MarkDuplicates () {
	
	printf $"\n~~~>	Sample $SAMPLE_NAME => Mark Duplicates\n\n"

	# java -Xmx64g -jar $PICARD MarkDuplicates \
	# I=$1 \
	# O=${1%.*.*}.mark.bam \
	# METRICS_FILE=${1%.*.*}.Metrics.txt \
	# READ_NAME_REGEX=null \
	# VALIDATION_STRINGENCY=LENIENT \
	# REMOVE_DUPLICATES=true \
	# ASSUME_SORTED=true
	java -Xmx64g -jar $PICARD MarkDuplicates \
	I=$1 \
	O=${1%.*.*}.mark.bam \
	METRICS_FILE=${1%.*.*}.Metrics.txt \
	READ_NAME_REGEX=null \
	VALIDATION_STRINGENCY=LENIENT \
	ASSUME_SORTED=true

	mv ${1%.*.*}.Metrics.txt $DELETE
	mv $1 $DELETE
	INPUT=${1%.*.*}.mark.bam

	printf $"\n~~~>	Sample $SAMPLE_NAME => Mark Duplicates: DONE\n\n"

	BuildBamIndex $INPUT
}

BuildBamIndex () {

	printf $"\n~~~>	$1 => Build Bam Index\n\n"

	java -Xmx64g -jar $PICARD BuildBamIndex \
	I=$1 \
	O=${1%.*}.bai \
	VALIDATION_STRINGENCY=LENIENT

	printf $"\n~~~>	$1 => Build Bam Index: DONE\n\n"
}

IndelRealigner () {

	printf $"\n~~~>	Sample $SAMPLE_NAME => Realigner Target Creator\n\n"

	java -Xmx64g -jar $GATK -T RealignerTargetCreator \
	-R $REF \
	-I $1 \
	-o ${1%.*.*}.IndelRealigner.intervals \
	--known $MILLS \
	-L $TARGET

	printf $"\n~~~>	Sample $SAMPLE_NAME => Realigner Target Creator: DONE\n\n"

	printf $"\n~~~>	Sample $SAMPLE_NAME => Indel Realigner\n\n"

	java -Xmx64g -jar $GATK -T IndelRealigner \
	-R $REF \
	-I $1 \
	-targetIntervals ${1%.*.*}.IndelRealigner.intervals \
	-o ${1%.*.*}.Realigned.bam \
	-known $MILLS

	mv ${1%.*.*}.IndelRealigner.intervals $DELETE
	mv $1 $DELETE
	mv ${1%.*}.bai $DELETE
	INPUT=${1%.*.*}.Realigned.bam

	printf $"\n =========>	Indel Realigner: DONE\n"
}

BaseRecalibrator () {

	printf $"\n~~~>	Sample $SAMPLE_NAME => Base Quality Score Recalibration (BQSR): Base Recalibrator\n\n"

	java -Xmx64g -jar $GATK -T BaseRecalibrator \
	-R $REF \
	-I $1 \
	-knownSites $DBSNP \
	-knownSites $MILLS \
	-o ${1%.*.*}.Recal.table \
	-L $TARGET

	printf $"\n~~~>	Sample $SAMPLE_NAME => Base Quality Score Recalibration (BQSR): Base Recalibrator: DONE\n\n"

	printf $"\n~~~>	Sample $SAMPLE_NAME => Base Quality Score Recalibration (BQSR): Print reads\n\n"
	
	java -Xmx64g -jar $GATK -T PrintReads \
	-R $REF \
	-I $1 \
	-BQSR ${1%.*.*}.Recal.table \
	-o ${1%.*.*}.bam \
	-L $TARGET


	mv ${1%.*.*}.Recal.table $DELETE
	mv $1 $DELETE
	mv ${1%.*}.bai $DELETE
	INPUT=${1%.*.*}.bam
	printf $"\n~~~>	Sample $SAMPLE_NAME => Base Quality Score Recalibration (BQSR): Print reads: DONE\n\n"
	
}

PREPROCESSING () {

	cat $LOGHI/logo_processing.txt
	mkdir -p $WORKDIR/PREPROCESSING
	rm -f $WORKDIR/PostPreprocessing.cfg
	CFG=$WORKDIR/PostPreprocessing.cfg

	if [ "$ANALISI" == "Germline" ] ||  [ "$ANALISI" == "Sanger" ] 
		then
		cat $1 | while read line
		do
			INPUT=$(echo "$line" | cut -f1)
			SAMPLE_NAME=$(echo "$line" | cut -f2)
					
			AddOrReplaceReadGroups $INPUT

			if [ "$DESIGN" == "ENRICHMENT" ]
			then
				if [[ "$START" == *"M"* ]]
				then
					MarkDuplicates $INPUT
				fi
				if [[ "$START" == *"I"* ]]
				then
					IndelRealigner $INPUT
				fi
				if [[ "$START" == *"B"* ]]
				then
					BaseRecalibrator $INPUT
				fi

				mv $INPUT $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
				mv ${INPUT%.*}.bai $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bai
				BuildBamIndex $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam

			elif [ "$DESIGN" == "AMPLICON" ]
			then
				mv ${INPUT%.*.*}.Add.bam $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
				BuildBamIndex $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
			fi
			printf $"$WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam\t$SAMPLE_NAME\n" >> $CFG
		done
	elif [ "$ANALISI" == "Somatic" ]
		then
		cat $1 | while read line
		do
			INPUT=$(echo "$line" | cut -f1)
			SAMPLE_NAME=$(echo "$line" | cut -f2)
					
			AddOrReplaceReadGroups $INPUT

			if [ "$DESIGN" == "ENRICHMENT" ]
			then
				if [[ "$START" == *"M"* ]]
				then
					MarkDuplicates $INPUT
				fi
				if [[ "$START" == *"I"* ]]
				then
					IndelRealigner $INPUT
				fi
				if [[ "$START" == *"B"* ]]
				then
					BaseRecalibrator $INPUT
				fi

				mv $INPUT $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
				mv ${INPUT%.*}.bai $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bai

			elif [ "$DESIGN" == "AMPLICON" ]
			then
				mv ${INPUT%.*.*}.Add.bam $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
				BuildBamIndex $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
			fi

			INPUT_SOM=$WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
			SAMPLE_NAME_SOM=$SAMPLE_NAME
			
			INPUT=$(echo "$line" | cut -f3)
			SAMPLE_NAME=$(echo "$line" | cut -f4)

			AddOrReplaceReadGroups $INPUT

			if [ "$DESIGN" == "ENRICHMENT" ]
			then
				if [[ "$START" == *"M"* ]]
				then
					MarkDuplicates $INPUT
				fi
				if [[ "$START" == *"I"* ]]
				then
					IndelRealigner $INPUT
				fi
				if [[ "$START" == *"B"* ]]
				then
					BaseRecalibrator $INPUT
				fi

				mv $INPUT $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
				mv ${INPUT%.*}.bai $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bai

			elif [ "$DESIGN" == "AMPLICON" ]
			then
				mv ${INPUT%.*.*}.Add.bam $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
				BuildBamIndex $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
			fi

			INPUT_NORM=$WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
			SAMPLE_NAME_NORM=$SAMPLE_NAME
			
			printf $"$INPUT_SOM\t$SAMPLE_NAME_SOM\t$INPUT_NORM\t$SAMPLE_NAME_NORM\n" >> $CFG
		done

	elif [ "$ANALISI" == "CellFree" ]
		then
		cat $1 | while read line
		do
			
			IFS=$':' DIRS=(${line//$'\t'/:})
			dastampare=()

			for (( i=0 ; i<${#DIRS[@]} ; i++ ))
			do
				INPUT=${DIRS[i]}
				SAMPLE_NAME=${DIRS[i+1]}
				((i+=1))

				AddOrReplaceReadGroups $INPUT

				if [ "$DESIGN" == "ENRICHMENT" ]
				then
					if [[ "$START" == *"M"* ]]
					then
						MarkDuplicates $INPUT
					fi
					if [[ "$START" == *"I"* ]]
					then
						IndelRealigner $INPUT
					fi
					if [[ "$START" == *"B"* ]]
					then
						BaseRecalibrator $INPUT
					fi

					mv $INPUT $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
					mv ${INPUT%.*}.bai $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bai

				elif [ "$DESIGN" == "AMPLICON" ]
				then
					mv ${INPUT%.*.*}.Add.bam $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
					BuildBamIndex $WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam
				fi

				dastampare+=("$WORKDIR/PREPROCESSING/$SAMPLE_NAME.bam\t$SAMPLE_NAME")
			done
			bar=$(IFS=$'\t' ; echo "${dastampare[*]}")
			echo -e "$bar">>$CFG
		done

	fi


}