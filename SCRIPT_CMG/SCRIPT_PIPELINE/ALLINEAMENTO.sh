BWAMEM () {
	printf $"\n~~~>	Sample $SAMPLE_NAME => BWA MEM\n\n"
	$BWA/bwa mem $REF \
	-M $1 \
	-t 2 \
	$2 > $WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam
	
	INPUT=$WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam
	printf $"\n~~~>	Sample $SAMPLE_NAME => BWA MEM: DONE\n\n"
}

SamFormatConverter () {
	printf $"\n~~~>	Sample $SAMPLE_NAME => Sam Format Converter\n\n"
	java -Xmx64g -jar $PICARD SamFormatConverter \
	I=$1 \
	O=${1%.*}.converted.bam
	
	mv $1 $DELETE
	INPUT=${1%.*}.converted.bam
	printf $"\n~~~>	Sample $SAMPLE_NAME => Sam Format Converter: DONE\n\n"
}

SortSam () {
	printf $"\n~~~>	Sample $SAMPLE_NAME => Sort Sam\n\n"
	java -Xmx64g -jar $PICARD SortSam \
	I=$1 \
	O=${1%.*.*}.sort.bam \
	SORT_ORDER=coordinate

	mv $1 $DELETE
	INPUT=${1%.*.*}.sort.bam
	printf $"\n~~~>	Sample $SAMPLE_NAME => Sort Sam: DONE\n\n"
}

ALLINEAMENTO() {

	cat $LOGHI/logo_alignment.txt

	mkdir -p $WORKDIR/ALIGNMENT
	mkdir -p $WORKDIR/ALIGNMENT/

	rm -f $WORKDIR/PostAlignment.cfg
	CFG=$WORKDIR/PostAlignment.cfg

	if [ "$ANALISI" == "Germline" ]
		then
		cat $1 | while read line
		do 
			FASTQ1=$(echo "$line" | cut -f1)
			FASTQ2=$(echo "$line" | cut -f2)
			SAMPLE_NAME=$(echo "$line" | cut -f3)
			
			BWAMEM $FASTQ1 $FASTQ2
			SamFormatConverter $INPUT
			SortSam $INPUT
			
			cp $FASTQ1 $STORAGE
			cp $FASTQ2 $STORAGE
			
			printf $"$INPUT\t$SAMPLE_NAME\n" >> $CFG
		done

	elif [ "$ANALISI" == "Somatic" ]
		then
		cat $1 | while read line
		do
			FASTQ1=$(echo "$line" | cut -f1)
			FASTQ2=$(echo "$line" | cut -f2)
			SAMPLE_NAME=$(echo "$line" | cut -f3)

			BWAMEM $FASTQ1 $FASTQ2
			SamFormatConverter $INPUT
			SortSam $INPUT
			INPUT_SOM=$INPUT
			SAMPLE_NAME_SOM=$SAMPLE_NAME

			cp $FASTQ1 $STORAGE
			cp $FASTQ2 $STORAGE

			#mv $WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam $DELETE
			#mv $WORKDIR/ALIGNMENT/$SAMPLE_NAME.converted.bam $DELETE
			
			FASTQ1=$(echo "$line" | cut -f4)
			FASTQ2=$(echo "$line" | cut -f5)
			SAMPLE_NAME=$(echo "$line" | cut -f6)
			
			BWAMEM $FASTQ1 $FASTQ2
			SamFormatConverter $INPUT
			SortSam $INPUT
			INPUT_NORM=$INPUT
			SAMPLE_NAME_NORM=$SAMPLE_NAME
			
			cp $FASTQ1 $STORAGE
			cp $FASTQ2 $STORAGE
			#mv $WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam $DELETE
			#mv $WORKDIR/ALIGNMENT/$SAMPLE_NAME.converted.bam $DELETE
			
			printf $"$INPUT_SOM\t$SAMPLE_NAME_SOM\t$INPUT_NORM\t$SAMPLE_NAME_NORM\n" >> $CFG
		done
	fi
	cat $LOGHI/logo_cornice.txt
}