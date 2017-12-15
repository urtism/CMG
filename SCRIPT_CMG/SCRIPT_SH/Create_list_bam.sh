#!bin/bash

cd /media/jarvis/HD2/Storage_Run_cmgcv/CANCER/

#Questo script serve per processare tutti i file dentro ogni cartella nel path /media/jarvis/HD2/Storage_Run_cmgcv/CARDIO/
#Estraggo poi la data prendendo solo il primo campo del nome ed infine scrivo tutti i path del bam nel file .list
#Per estrarre il path assoluto del file uso il comando ----> ls -d -1 $PWD/*.bam

for d in */
	do	mkdir /media/jarvis/HD1/VCF_Ann_iEVA_Cancer/$d
		cd $d
		extract=$(echo $d| cut -d'_' -f 1)
		#rm $extract.list
		ls -d -1 $PWD/*.bam >/media/jarvis/HD1/VCF_Ann_iEVA_Cancer/$d/$extract.list
		echo "$d"
		cd ..
	done
