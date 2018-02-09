#!bin/bash

cd /media/matteo/CHEWBWCCA/JARVIS-STORAGE/ANNOTATED_VCF_iEVA/BRCA/

#Questo script serve per processare tutti i file dentro ogni cartella nel path /media/jarvis/HD2/Storage_Run_cmgcv/CARDIO/
#Estraggo poi la data prendendo solo il primo campo del nome ed infine scrivo tutti i path del bam nel file .list
#Per estrarre il path assoluto del file uso il comando ----> ls -d -1 $PWD/*.bam

for d in */
	do	
		newD=$d
		replace='2.1_VCF_IEVA'
		result_string="${newD/2.1/$replace}"
		#echo $result_string
		mv /media/matteo/CHEWBWCCA/JARVIS-STORAGE/ANNOTATED_VCF_iEVA/BRCA/$d/*.list /media/matteo/CHEWBWCCA/JARVIS-STORAGE/OUT_VCF_IEVA_BRCA/$result_string/
		#cd $d
		#extract=$(echo $d| cut -d'_' -f 1)
		#rm $extract.list
		#ls -d -1 $PWD/*.bam > /media/matteo/CHEWBWCCA/JARVIS-STORAGE/ANNOTATED_VCF_iEVA/BRCA/$d/$extract.list
		#echo "$d"
		#cd ..
	done
