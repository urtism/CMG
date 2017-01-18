#!bin/bash

DATA=20161223

TIPO=BRCA

NOMEDIR=BRCA

Pannello=Cardio

Upload=2.1

NumRun=70

cd ~/Scrivania/RUN_ILLUMINA/$NOMEDIR/$DATA\_Run_$NumRun\_$Pannello\_$Upload
for (( a=1; a<10; a++ ))
do
mv $DATA-0$a-$TIPO\_S$a\_L001_R1_001.fastq.gz $DATA\_0$a\_$TIPO\_L001_R1_001.fastq.gz
mv $DATA-0$a-$TIPO\_S$a\_L001_R2_001.fastq.gz $DATA\_0$a\_$TIPO\_L001_R2_001.fastq.gz
done

cd ~/Scrivania/RUN_ILLUMINA/$NOMEDIR/$DATA\_Run_$NumRun\_$Pannello\_$Upload
for (( a=10; a<16; a++ ))
do
mv $DATA-$a-$TIPO\_S$a\_L001_R1_001.fastq.gz $DATA\_$a\_$TIPO\_L001_R1_001.fastq.gz
mv $DATA-$a-$TIPO\_S$a\_L001_R2_001.fastq.gz $DATA\_$a\_$TIPO\_L001_R2_001.fastq.gz
done
