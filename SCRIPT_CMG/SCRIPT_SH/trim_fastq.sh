CFG=$2/config.cfg

cat $1 | while read line
do
	R1=$(echo "$line" | cut -f1)
	R2=$(echo "$line" | cut -f2)
	sample_name=$(echo "$line" | cut -f3)

	sample_name_R1=${R1##*/}
	sample_name_R2=${R2##*/}

	cutadapt -a GAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	$R1 $R2 \
	-o $2/$sample_name_R1 -p $2/$sample_name_R2 -m 1

	printf $"$2/$sample_name_R1\t$2/$sample_name_R2\t$sample_name\t" >> $CFG


	R1=$(echo "$line" | cut -f4)
	R2=$(echo "$line" | cut -f5)
	sample_name=$(echo "$line" | cut -f6)

	sample_name_R1=${R1##*/}
	sample_name_R2=${R2##*/}

	cutadapt -a GAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	$R1 $R2 \
	-o $2/$sample_name_R1 -p $2/$sample_name_R2 -m 1


	printf $"$2/$sample_name_R1\t$2/$sample_name_R2\t$sample_name\n" >> $CFG
done