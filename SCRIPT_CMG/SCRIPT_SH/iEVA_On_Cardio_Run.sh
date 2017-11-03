#!bin/bash

REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta

cd /media/jarvis/HD2/Storage_Run_cmgcv/CARDIO

#Questo script serve per lanciare iEVA su tutte le run contenute in /media/jarvis/HD2/Storage_Run_cmgcv/CARDIO/
#Poi vado a prendere il vcf in /home/jarvis/Scrivania/VCF-IEVA/NORM ed infine calcolo il tempo computazionale.

for d in */
	do

		#DATE=$(date +"%Y%m%d");
		STARTTIME=$(date +%s)
		extract=$(echo $d| cut -d'_' -f 1)
		run=$(echo $d| cut -d'_' -f 3)
		printf '\n'
		echo "Annoto con iEVA i file della run $d"

		#LANCIO iEVA SU FREEBAYES
		STARTTIME_FREE=$(date +%s)

		python ~/iEVA/iEVA/__main__.py \
		-Ref $REF \
		-I /home/jarvis/Scrivania/VCF-IEVA/NORM/$extract\_*FreeBayes*.vcf \
		-L /media/jarvis/HD2/Storage_Run_cmgcv/CARDIO/$d/$extract.list \
		-O /media/jarvis/HD2/Storage_Run_cmgcv/VCF_Ann_iEVA_Default/$extract\_Cardio_iEVA_FreeBayes_Default.vcf \
		-v \
		-WS 40 \
		-SR \
		-SRL \
		-SRU \
		-PNC \
		-RM \
		-GC \
		-VC \
		-SBR \
		-UnMap \
		-MMQ \
		-MQ0 \
		-NPA \
		-SA \
		-NP \
		-NPP \
		-AS \
		-XS \
		-DUP \
		-TDP \
		-iDP \
		-iAD \
		-iADup \
		-iDDup \
		-iQual \
		-iMMQ \
		-iAS \
		-iXS \
		-iXS0 \
		-iUnMap \
		-iMQ0 \
		-iNPA \
		-iSA \
		-iNP \
		-iNPP \
		-iCR

		ENDTIME_FREE=$(date +%s)
		PROCESSING_FREE=$(($ENDTIME_FREE - $STARTTIME_FREE))


		#LANCIO iEVA SU VARSCAN
		STARTTIME_VARSCAN=$(date +%s)

		python ~/iEVA/iEVA/__main__.py \
		-Ref $REF \
		-I /home/jarvis/Scrivania/VCF-IEVA/NORM/$extract\_*VarScan*.vcf \
		-L /media/jarvis/HD2/Storage_Run_cmgcv/CARDIO/$d/$extract.list \
		-O /media/jarvis/HD2/Storage_Run_cmgcv/VCF_Ann_iEVA_Default/$extract\_Cardio_iEVA_VarScan_Default.vcf \
		-v \
		-WS 40 \
		-SR \
		-SRL \
		-SRU \
		-PNC \
		-RM \
		-GC \
		-VC \
		-SBR \
		-UnMap \
		-MMQ \
		-MQ0 \
		-NPA \
		-SA \
		-NP \
		-NPP \
		-AS \
		-XS \
		-DUP \
		-TDP \
		-iDP \
		-iAD \
		-iADup \
		-iDDup \
		-iQual \
		-iMMQ \
		-iAS \
		-iXS \
		-iXS0 \
		-iUnMap \
		-iMQ0 \
		-iNPA \
		-iSA \
		-iNP \
		-iNPP \
		-iCR

		ENDTIME_VARSCAN=$(date +%s)
		PROCESSING_VARSCAN=$(($ENDTIME_VARSCAN - $STARTTIME_VARSCAN))


		#LANCIO iEVA SU GATK
		STARTTIME_GATK=$(date +%s)

		python ~/iEVA/iEVA/__main__.py \
		-Ref $REF \
		-I /home/jarvis/Scrivania/VCF-IEVA/NORM/$extract\_*GATK*.vcf \
		-L /media/jarvis/HD2/Storage_Run_cmgcv/CARDIO/$d/$extract.list \
		-O /media/jarvis/HD2/Storage_Run_cmgcv/VCF_Ann_iEVA_Default/$extract\_Cardio_iEVA_GATK_Default.vcf \
		-v \
		-WS 40 \
		-SR \
		-SRL \
		-SRU \
		-PNC \
		-RM \
		-GC \
		-VC \
		-SBR \
		-UnMap \
		-MMQ \
		-MQ0 \
		-NPA \
		-SA \
		-NP \
		-NPP \
		-AS \
		-XS \
		-DUP \
		-TDP \
		-iDP \
		-iAD \
		-iADup \
		-iDDup \
		-iQual \
		-iMMQ \
		-iAS \
		-iXS \
		-iXS0 \
		-iUnMap \
		-iMQ0 \
		-iNPA \
		-iSA \
		-iNP \
		-iNPP \
		-iCR

		ENDTIME_GATK=$(date +%s)
		PROCESSING_GATK=$(($ENDTIME_GATK - $STARTTIME_GATK))

		cd ..
		ENDTIME=$(date +%s)
		PROCESSING=$(($ENDTIME - $STARTTIME))
		echo "$PROCESSING"
		printf '\n'
		printf "\n${extract}_Run_${run}\t${PROCESSING_GATK}\t${PROCESSING_FREE}\t${PROCESSING_VARSCAN}\t${PROCESSING}\t\n" >> /media/jarvis/HD2/Storage_Run_cmgcv/iEVA_log_DEFAULT.txt

	done