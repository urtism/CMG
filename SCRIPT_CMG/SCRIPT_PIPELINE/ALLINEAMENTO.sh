ALLINEAMENTO() {

	cat $LOGHI/logo_alignment.txt

	mkdir -p $WORKDIR/Alignment
	mkdir -p $WORKDIR/Alignment/Aligned_Sam/

	cat $FILEPATH | while read line
	do

		FASTQ1=$(echo "$line" | cut -f1)
		FASTQ2=$(echo "$line" | cut -f2)
		SAMPLE_NAME=$(echo "$line" | cut -f3)

		printf $"\n~~~>	Sample $SAMPLE_NAME => BWA MEM\n\n"
		$BWA/bwa mem $REF \
		-M $FASTQ1 \
		-t 12 \
		$FASTQ2 > $WORKDIR/Alignment/$SAMPLE_NAME.sam
		printf $"\n~~~>	Sample $SAMPLE_NAME => BWA MEM: DONE\n\n"

		#$FASTQC -f sam -o $WORKDIR/Alignment $WORKDIR/Alignment/$SAMPLE_NAME.sam
		printf $"\n~~~>	Sample $SAMPLE_NAME => Sam Format Converter\n\n"
		java -Xmx64g -jar $PICARD SamFormatConverter \
		I=$WORKDIR/Alignment/$SAMPLE_NAME.sam\
		O=$WORKDIR/Alignment/$SAMPLE_NAME.conv.bam
		printf $"\n~~~>	Sample $SAMPLE_NAME => Sam Format Converter: DONE\n\n"
		
		printf $"\n~~~>	Sample $SAMPLE_NAME => Sort Sam\n\n"
		java -Xmx64g -jar $PICARD SortSam \
		I=$WORKDIR/Alignment/$SAMPLE_NAME.conv.bam \
		O=$WORKDIR/Alignment/$SAMPLE_NAME.sort.bam \
		SORT_ORDER=coordinate
		printf $"\n~~~>	Sample $SAMPLE_NAME => Sort Sam: DONE\n\n"

		rm $WORKDIR/Alignment/$SAMPLE_NAME.conv.sam
		rm $WORKDIR/Alignment/$SAMPLE_NAME.sam
	
	done
	cat ~/Scrivania/SCRIPT_PIPELINE/logo_cornice.txt
}