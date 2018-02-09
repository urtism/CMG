BWAMEM () {
	printf $"\n~~~>	Sample $SAMPLE_NAME => BWA MEM\n\n"
	
	$BWA/bwa mem $REF \
	$1 $2 \
	-t 2 \
	 > $WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam
	
	INPUT=$WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam
	printf $"\n~~~>	Sample $SAMPLE_NAME => BWA MEM: DONE\n\n"
}

BWASW () {
	printf $"\n~~~>	Sample $SAMPLE_NAME => BWASW\n\n"
	
	$BWA/bwa bwasw $REF \
	$1 $2 \
	-t 2 \
	 > $WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam
	
	INPUT=$WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam
	printf $"\n~~~>	Sample $SAMPLE_NAME => BWASW: DONE\n\n"
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

SureCallTrimmer () {
	printf $"\n~~~>	Sample $SAMPLE_NAME => Trimming HaloPlex Adapters\n\n"
	
	java -jar $SURECALLTRIMMER -fq1 $1 -fq2 $2 $3 -out_loc $WORKDIR/ALIGNMENT/FASTQTRIMM/  

	trim1=$(basename $1 .gz)
	trim2=$(basename $2 .gz)

	FASTQ1="$(ls $WORKDIR/ALIGNMENT/FASTQTRIMM/$trim1*)"
	FASTQ2="$(ls $WORKDIR/ALIGNMENT/FASTQTRIMM/$trim2*)"
	printf $"\n~~~>	Sample $SAMPLE_NAME => Trimming HaloPlex Adapters: DONE\n\n"
}

LocatIt (){
	printf $"\n~~~>	Sample $SAMPLE_NAME => Merging reads by Molecular BarCodes\n\n"
	
	java -Xmx12G -jar $LOCATLT -U -IB -OB -i -o ${1%.*.*}.MBCmerged.bam $1 $2 

	mv $1 $DELETE
	INPUT=${1%.*.*}.MBCmerged.bam

	printf $"\n~~~>	Sample $SAMPLE_NAME => Merging reads by Molecular BarCodes: DONE\n\n"
}

CutAdapt (){
	printf $"\n~~~>	Sample $SAMPLE_NAME => Trimming HaloPlex Adapters\n\n"
	
	cutadapt -a GAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	$1 $2 \
	-o $WORKDIR/ALIGNMENT/FASTQTRIMM/$(basename $1) -p $WORKDIR/ALIGNMENT/FASTQTRIMM/$(basename $2) -m 1

	FASTQ1=$WORKDIR/ALIGNMENT/FASTQTRIMM/$(basename $1)
	FASTQ2=$WORKDIR/ALIGNMENT/FASTQTRIMM/$(basename $2)
	printf $"\n~~~>	Sample $SAMPLE_NAME => Trimming HaloPlex Adapters: DONE\n\n"
}

ALLINEAMENTO() {

	cat $LOGHI/logo_alignment.txt

	mkdir -p $WORKDIR/ALIGNMENT
	mkdir -p $WORKDIR/ALIGNMENT/

	rm -f $WORKDIR/PostAlignment.cfg
	CFG=$WORKDIR/PostAlignment.cfg

	if [ "$ANALISI" == "Germline" ] || [ "$ANALISI" == "SingleSample" ]
		then
		cat $1 | while read line
		do 	
			if [ "$PANNELLO" == "HaloPlex" ]
			then

				FASTQ1=$(echo "$line" | cut -f1)
				FASTQ2=$(echo "$line" | cut -f2)
				INDEX=$(echo "$line" | cut -f3)
				SAMPLE_NAME=$(echo "$line" | cut -f4)

				SureCallTrimmer $FASTQ1 $FASTQ2 -hs
				BWAMEM $FASTQ1 $FASTQ2
				SamFormatConverter $INPUT
				LocatIt $INPUT $INDEX
				SortSam $INPUT

			elif [ "$PANNELLO" == "SureSelect" ]
			then

				FASTQ1=$(echo "$line" | cut -f1)
				FASTQ2=$(echo "$line" | cut -f2)
				SAMPLE_NAME=$(echo "$line" | cut -f3)

				#SureCallTrimmer $FASTQ1 $FASTQ2 -qxt
				BWAMEM $FASTQ1 $FASTQ2
				SamFormatConverter $INPUT
				SortSam $INPUT

			else
				FASTQ1=$(echo "$line" | cut -f1)
				FASTQ2=$(echo "$line" | cut -f2)
				SAMPLE_NAME=$(echo "$line" | cut -f3)
				
				BWAMEM $FASTQ1 $FASTQ2
				SamFormatConverter $INPUT
				SortSam $INPUT
			fi
			
			printf $"$INPUT\t$SAMPLE_NAME\n" >> $CFG
		done

	elif [ "$ANALISI" == "Somatic" ]
		then
		cat $1 | while read line
		do
			if [ "$PANNELLO" == "HaloPlex" ]
				then
				FASTQ1=$(echo "$line" | cut -f1)
				FASTQ2=$(echo "$line" | cut -f2)
				INDEX=$(echo "$line" | cut -f3)
				SAMPLE_NAME=$(echo "$line" | cut -f4)

				SureCallTrimmer $FASTQ1 $FASTQ2 -hs
				BWAMEM $FASTQ1 $FASTQ2
				SamFormatConverter $INPUT
				LocatIt $INPUT $INDEX
				SortSam $INPUT
				INPUT_SOM=$INPUT
				SAMPLE_NAME_SOM=$SAMPLE_NAME

				cp $FASTQ1 $STORAGE
				cp $FASTQ2 $STORAGE

				#mv $WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam $DELETE
				#mv $WORKDIR/ALIGNMENT/$SAMPLE_NAME.converted.bam $DELETE
				
				FASTQ1=$(echo "$line" | cut -f5)
				FASTQ2=$(echo "$line" | cut -f6)
				INDEX=$(echo "$line" | cut -f7)
				SAMPLE_NAME=$(echo "$line" | cut -f8)

				SureCallTrimmer $FASTQ1 $FASTQ2 -hs
				BWAMEM $FASTQ1 $FASTQ2
				SamFormatConverter $INPUT
				LocatIt $INPUT $INDEX
				SortSam $INPUT
				INPUT_NORM=$INPUT
				SAMPLE_NAME_NORM=$SAMPLE_NAME
				
				cp $FASTQ1 $STORAGE
				cp $FASTQ2 $STORAGE
				#mv $WORKDIR/ALIGNMENT/$SAMPLE_NAME.sam $DELETE
				#mv $WORKDIR/ALIGNMENT/$SAMPLE_NAME.converted.bam $DELETE
				
				printf $"$INPUT_SOM\t$SAMPLE_NAME_SOM\t$INPUT_NORM\t$SAMPLE_NAME_NORM\n" >> $CFG
			else
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
			fi
		done

	elif [ "$ANALISI" == "CellFree" ]
		then
		
		cat $1 | while read line
		do
			dastampare=()
			IFS=$':' DIRS=(${line//$'\t'/:})

			for (( i=0 ; i<${#DIRS[@]} ; i++ ))
			do
				FASTQ1=${DIRS[i]}
				FASTQ2=${DIRS[i+1]}
				SAMPLE_NAME=${DIRS[i+2]}


				BWAMEM $FASTQ1 $FASTQ2
				SamFormatConverter $INPUT
				SortSam $INPUT
				INPUT_CF=$INPUT
				SAMPLE_NAME_CF=$SAMPLE_NAME

				dastampare+=("$INPUT_CF\t$SAMPLE_NAME_CF")
				((i+=2))
			done
			bar=$(IFS=$'\t' ; echo "${dastampare[*]}")
			echo -e "$bar">>$CFG
		done

	elif [ "$ANALISI" == "Sanger" ]
		then
		cat $1 | while read line
		do 
			FASTQ1=$(echo "$line" | cut -f1)
			SAMPLE_NAME=$(echo "$line" | cut -f2)
			
			BWASW $FASTQ1
			SamFormatConverter $INPUT
			SortSam $INPUT
			
			cp $FASTQ1 $STORAGE
			cp $FASTQ2 $STORAGE
			
			printf $"$INPUT\t$SAMPLE_NAME\n" >> $CFG
		done

	fi
	cp $FASTQ1 $STORAGE
	cp $FASTQ2 $STORAGE

	cat $LOGHI/logo_cornice.txt
}