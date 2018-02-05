#!bin/bash

REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta

cd /media/minime/IANSOLO/JARVIS-STORAGE/OUT_VCF_IEVA/

#Questo script serve per lanciare iEVA su tutte le run contenute in /media/jarvis/HD2/Storage_Run_cmgcv/CARDIO/
#Poi vado a prendere il vcf in /home/jarvis/Scrivania/VCF-IEVA/NORM ed infine calcolo il tempo computazionale.

# for file in *.vcf
# 	do

# 		DATA=$(echo $file| cut -d'_' -f 1)
# 		echo "$DATA"

# 		mv $file $DATA\_*\_IEVA

# 	done

for d in */
	do
		cd $d

		rename 's/FreeBayes.vcf/FreeBayes.norm.vcf/' *.vcf
		rename 's/VarScan.vcf/VarScan.norm.vcf/' *.vcf

		ls *VarScan.norm.vcf
		cd ..

	done
