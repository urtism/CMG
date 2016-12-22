#!bin/bash

cd /home/minime/Scrivania/STORICO_CONFERME_RUN/CARDIO
#cd /home/minime/Scrivania/a

for dir in *	
do
echo $dir
cd /home/minime/Scrivania/STORICO_CONFERME_RUN/CARDIO/$dir
	for file in *.ods
	do
	mv "$file" "${file%.ods}.tsv"
	done
done
