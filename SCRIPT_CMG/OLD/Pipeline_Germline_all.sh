#!/bin/bash

cat ~/Scrivania/SCRIPT_PIPELINE/logo.txt 

FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.15
PICARD=~/NGS_TOOLS/picard-tools-2.7.1/picard.jar
GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
BCFTOOLS=~/NGS_TOOLS/bcftools-1.3.1/bcftools
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_138.hg19.vcf
INPUT=~/NGS_ANALYSIS/INPUT_DATA/CARDIO
INPUTBRCA=~/NGS_ANALYSIS/INPUT_DATA/BRCA
INPUTCANCER=~/NGS_ANALYSIS/INPUT_DATA/CANCER
INPUTEXOME=~/NGS_ANALYSIS/INPUT_DATA/EXOME
PROCESSING=~/NGS_ANALYSIS/PROCESSING
TARGET=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.list
TARGETCARDIOBED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.bed
TARGETMETRICS=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a.list
TARGETBRCA=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.list
TARGETBRCABED=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.bed
TARGBRCAFREE=~/NGS_ANALYSIS/TARGET/BRCA_FreeBayes_amplicon.bed
TARGETEXOME=~/NGS_ANALYSIS/TARGET/TruSight_One_v1.1_ESTESO+-1000.list
TARGETEXOMEBED=~/NGS_ANALYSIS/TARGET/TruSight_One_v1.1_ESTESO+-1000.bed
TARGETCANCER=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.list
TARGETCANCERBED=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.bed
OUTVCF=~/NGS_ANALYSIS/OUTPUT_DATA
STORAGE=~/NGS_ANALYSIS/STORAGE
SCRIPTS_DIR=~/git/CMG/scripts
VEP=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/filter_vep.pl

DATE=$(date +"%Y%m%d");
STARTTIME=$(date +%s)

cd $PROCESSING

echo $'\nPipeline NGS sviluppata da: Ing. Matteo Di Giovannantonio\n\nCentro Malattie Genetiche Cardiovascolari. Fondazione IRCCS Policlinico San Matteo. \n\n'

printf "Prima di iniziare, ti spiegherò come far partire l'analisi.\n\n 1) Premessa: quando parlo di analisi CARDIO intendo sia connettivi che cardio, mentre un'analisi a parte viene fatta per i BRCA.\nPer prima cosa ti verrà chiesto di inserire la data della run. N.B per data run si intende la data corrispondente alla run e non la data attuale! Se ad esempio la run è stata eseguita il 10/08/2015 (20150810) ma stai eseguendo l'analisi nella data odierna, inserisci 20150810 come data e non $DATE.\n\n 2) Prima di ogni cosa verifica che i file siano nominati nel seguente formato: Data_paziente_analisi.estensione Ad esempio, un file deve essere nominato come ${DATE}_01_Cardio o BRCA o CONN quindi ad esempio i file .bam dovranno essere così nominati:\n\n ${DATE}_01_Cardio.bam \n ${DATE}_02_Conn.bam \n ${DATE}_03_BRCA.bam \n ${DATE}_04_Cardio.bam \n ${DATE}_05_Cardio.bam \n ${DATE}_06_Conn.bam \n ${DATE}_07_Cardio.bam \n ${DATE}_08_BRCA.bam \n ${DATE}_09_Cardio.bam \n ${DATE}_10_Cardio.bam \n ${DATE}_11_Cardio.bam \n ${DATE}_12_Cardio.bam \n\nNota bene che i pazienti da 1 a 9 hanno uno 0 prima.\n\n 3) Successivamente ti verrà chiesto quale tipo di analisi vuoi svolgere, se Cardio o BRCA. Ricorda che non puoi effettuare analisi cardio e BRCA contemporaneamente, gli step di processing sono diversi! \n\n 4) Per concludere ti verrà chiesto da quale step della pipeline preferisci iniziare, se effettuare l'analisi completa a partire dall'allineamento oppure dal pre-processing o infine dal solo variant-calline + annotazione. A seconda dello step dipenderà il path in cui dovrai caricare i file. Di seguito le istruzioni:\n\n -------------------CARDIO------------------- \n\n --> Analisi completa (Partire dall'allineamento sino all'annotazione):\n	Caricare i file fastq con il seguente nome: ${DATE}_n°paziente_Cardio_L001_R1_001.fastq.gz oppure ${DATE}_n°paziente_Conn_L001_R1_001.fastq.gz (se si sta analizzando un connettivo) all'interno del path ~/NGS_ANALYSIS/INPUT_DATA/CARDIO\n\n --> Analisi a partire dal pre-processing sino all'annotazione:\n	Caricare i file nominati come ${DATE}_n°paziente_Cardio.bam (o Conn) nel path $PROCESSING/1_Align/Sort_bam\n\n --> Analisi a partire dal variant-calling sino all'annotazione:\n	Caricare i file nominati come ${DATE}_n°paziente_Cardio.bam (o Conn) nel path $PROCESSING/5_BQSR con i rispettivi file .bai es:\n\n ${DATE}_01_Cardio.bam \n ${DATE}_01_Cardio.bai \n\n\nN.B -----> I file bam dal pre-processing e quelli dal variant calling NON SONO GLI STESSI! \n\n --------------------BRCA-------------------- \n\n --> Analisi completa (Partire dall'allineamento sino all'annotazione):\n	Caricare i file fastq con il seguente nome: ${DATE}_n°paziente_BRCA_L001_R1_001.fastq.gz all'interno del path ~/NGS_ANALYSIS/INPUT_DATA/BRCA\n\n --> Analisi a partire dal pre-processing sino all'annotazione:\n	Caricare i file nominati come ${DATE}_n°paziente_BRCA.bam nel path $PROCESSING/1_Align/Sort_bam\n\n --> Analisi a partire dal variant-calling sino all'annotazione:\n	Caricare i file nominati come ${DATE}_n°paziente_BRCA.bam nel path $PROCESSING/4_Indel con i rispettivi file .bai es:\n\n ${DATE}_01_BRCA.bam \n ${DATE}_01_BRCA.bai \n\n\nN.B -----> I file bam dal pre-processing e quelli dal variant calling NON SONO GLI STESSI!\n\n\nBUON DIVERTIMENTO! La pipeline ha inizio!\n\n\n"

echo -n "Inserisci la data della run corrispondente a quella presente nel nome file (es: $DATE). La data inserita verrà utilizzata per nominare il file finale Multi-sample (es: ${DATE}_Cardio_GATK.vcf) e per creare la directory da caricare in Google Genomics [Enter/Invio]:"

read DataRun

if [[ $DataRun =~ ^([0-9]{4})(0[1-9]|1[0-2])(0[1-9]|[1-2][0-9]|3[0-1])$  ]]; 
then
   	printf "\n\nHai inserito una data valida.\n\n"
else
	printf "\n\n"
	echo "Formato data non corretto. Inserire la data nel formato yyyymmdd (AnnoMeseGiorno) (es: $DATE). Arrivederci!"
	printf "\n\n"
exit
fi

echo -n "Inserisci il numero della run corrispondente. Il numero inserito verrà utilizzato per nominare la cartella storage da caricare in Google Genomics [Enter/Invio]:"

read NumRun

PipeVersion=2.1

#-------------------------------------------------SELEZIONE DEL TIPO DI ANALISI--------------------------------------------------------


printf "\n"
echo -n "Selezionare il tipo di analisi che si sta svolgendo: Cardio-Connettivo (C) [Enter/Invio]:"

read TIPO


#---------------------------------------------------------------CARDIO-----------------------------------------------------------------


if [ "$TIPO" == "C" ] || [ "$TIPO" == "c" ]
then

cd $OUTVCF
Name_OUT=$DataRun\_Run_$NumRun\_Cardio_Conferme_$PipeVersion
mkdir "$Name_OUT"
cd $STORAGE
Name_Dir=$DataRun\_Run_$NumRun\_Cardio_$PipeVersion
mkdir "$Name_Dir"

	printf "\n\n"
	echo -n "Potrai usare questa analisi soltanto a partire dal variant calling. Continuare? (y/n) "

	read Analisi

	until [ "$Analisi" == "N" ] || [ "$Analisi" == "n" ] || [ "$Analisi" == "Y" ] || [ "$Analisi" == "y" ]
	do
	printf "\nInserire solo Y/y/N/n [Enter/Invio]: "
	read Analisi
	done
	printf "$Analisi"


		if [ "$Analisi" == "N" ] || [ "$Analisi" == "n" ]
		then
		printf "\n\n"
		echo "Arrivederci!"
		printf "\n\n"
		exit
		fi


		if [ "$Analisi" == "Y" ] || [ "$Analisi" == "y" ]
		then
		printf "\n\nL'analisi avrà inizio a partire dal Variant calling.\n\n"
		cd $PROCESSING/5_BQSR/
		#Controllo che esistano file bam da cui partire con l'analisi
		count=$(ls -1 *.bam 2>/dev/null | wc -l)
		echo -n "Quanti sample stai analizzando? "
		read Sample
			if [[ $count != "$Sample" ]];
			then 
			printf "\n\nIl numero di bam non corrisponde al numero di pazienti che hai inserito. I file sono i seguenti:\n\n"
			ls *.bam
			printf "\n\nArrivederci!\n\n"
			exit 1;
			else
			printf "\n\nIl numero di bam è uguale al numero di pazienti. I file sono i seguenti:\n\n"
			ls *.bam
			printf "\n"
			ls *.bai
			printf "\n\n"
			echo -n "Confermi che i file sono corretti (S/N)? Ricorda che devono essere presenti anche i file .bai per ciascun .bam : "
			read Answer
				until [ "$Answer" == "S" ] || [ "$Answer" == "s" ] || [ "$Answer" == "N" ] || [ "$Answer" == "n" ]
				do
				printf "\nInserire Sì (S/s) o No (N/n): "
				read Answer
				done
					if [ "$Answer" == "S" ] || [ "$Answer" == "s" ]
					then
					printf "\nL'analisi ha inizio! Mettiti comodo, JARVIS farà il resto!\n\n"
					else
					printf "\nCarica i file corretti nel path $PROCESSING/5_BQSR e poi torna da me!\n\n"
					exit 1;
					fi
			fi

		source ~/Scrivania/SCRIPT_PIPELINE/Pipeline_ALL_CARDIO/Cardio_function_all.sh
		CARDIO_VARIANT_CALLING_AND_ANNOTATION
		fi


else
printf "\n\n"
echo "Error! Non esiste alcuna analisi del tipo $TIPO. Arrivederci!"
printf "\n\n"
exit
fi



ENDTIME=$(date +%s)
echo "TEMPO TOTALE DI PROCESSING: $(($ENDTIME - $STARTTIME)) sec"
printf "\n\n"

