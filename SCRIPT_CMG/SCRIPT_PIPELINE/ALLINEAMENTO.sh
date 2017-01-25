ALLINEAMENTO() {

	cat $LOGHI/logo_alignment.txt

	mkdir -p $WORKDIR/ALIGNMENT
	mkdir -p $WORKDIR/ALIGNMENT/

	rm -f $WORKDIR/PostAlignment.cfg
	CFG=$WORKDIR/PostAlignment.cfg

	cat $1 | while read line
	do

		FASTQ1=$(echo "$line" | cut -f1)
		FASTQ2=$(echo "$line" | cut -f2)
		SAMPLE_NAME=$(echo "$line" | cut -f3)

		printf $"\n~~~>	Sample $SAMPLE_NAME => BWA MEM\n\n"
		$BWA/bwa mem $REF \
		-M $FASTQ1 \
		-t  \
		$FASTQ2 > $WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam
		printf $"\n~~~>	Sample $SAMPLE_NAME => BWA MEM: DONE\n\n"

		#$FASTQC -f sam -o $WORKDIR/ALIGNMENT $WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam
		printf $"\n~~~>	Sample $SAMPLE_NAME => Sam Format Converter\n\n"
		java -Xmx64g -jar $PICARD SamFormatConverter \
		I=$WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam \
		O=$WORKDIR/ALIGNMENT/$SAMPLE_NAME.conv.bam
		printf $"\n~~~>	Sample $SAMPLE_NAME => Sam Format Converter: DONE\n\n"
		
		printf $"\n~~~>	Sample $SAMPLE_NAME => Sort Sam\n\n"
		java -Xmx64g -jar $PICARD SortSam \
		I=$WORKDIR/ALIGNMENT/$SAMPLE_NAME.conv.bam \
		O=$WORKDIR/ALIGNMENT/$SAMPLE_NAME.sort.bam \
		SORT_ORDER=coordinate
		printf $"\n~~~>	Sample $SAMPLE_NAME => Sort Sam: DONE\n\n"

		#rm $WORKDIR/ALIGNMENT/$SAMPLE_NAME.conv.sam
		rm $WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam
		
		printf $"$WORKDIR/ALIGNMENT/$SAMPLE_NAME.sort.bam\t$SAMPLE_NAME\n" >> $CFG
	done
	cat $LOGHI/logo_cornice.txt

}