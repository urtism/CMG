#!/bin/bash

#cat ~/Scrivania/SCRIPT_PIPELINE/logo.txt 


#echo $'\nPipeline NGS sviluppata da: Ing. Matteo Di Giovannantonio\n\nCentro Malattie Genetiche Cardiovascolari. Fondazione IRCCS Policlinico San Matteo. \n\n'

PIPELINE_GERMLINE() {

	ALLINEAMENTO

	if [ "$PANNELLO" == "Cardio" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CARDIO_1000
		echo "cardio"

	elif [ "$PANNELLO" == "Cancer" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CANCER_1000
		echo "cancer"

	elif [ "$PANNELLO" == "Exome" ] 
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_EXOME_1000
		echo "esoma"

	elif [ "$PANNELLO" == "BRCA" ]
	then
		DESIGN="AMPLICON"
		TARGET=$TARGET_BRCA
		echo "brca"

	fi

	PREPROCESSING

}



# 	#---------------------------------------------------------------CARDIO-----------------------------------------------------------------
# 		printf "\n\n"
# 		echo -n "Vuoi partire dall'allineamento (A/a), dal pre-processing (P/p) o direttamente dal variant-calling (V/v) ? ==> N.B Se si decide di saltare l'allineamento, verificare prima che i files .bam allineati siano contenuti nella cartella $PROCESSING/1_Align/Sort_bam con estensione DataRun_n°paziente_tipo-analisi (es: ${DataRun}_01_Cardio). I file da analizzare dovranno essere già ordinati con il comando SortSam del Picardtool. Qualora invece si parta dal variant Calling, verificare che i file siano presenti nella cartella $PROCESSING/5_BQSR. In questo caso i file utilizzati devono già essere stati sottoposti a pre-processing! [Enter/Invio]: "

# 		read Analisi

# 		until [ "$Analisi" == "A" ] || [ "$Analisi" == "a" ] || [ "$Analisi" == "P" ] || [ "$Analisi" == "p" ] || [ "$Analisi" == "V" ] || [ "$Analisi" == "v" ]
# 		do
# 		printf "\nInserire (A/a) per partire dall'allineamento, (P/p) per partire dal pre-processing oppure (V/v) per iniziare con il variant calling [Enter/Invio]: "
# 		read Analisi
# 		done
# 		printf "$Analisi"


# 			if [ "$Analisi" == "A" ] || [ "$Analisi" == "a" ]
# 			then
# 			printf "\n\nL'analisi avrà inizio a partire dall'allineamento. Mettiti comodo, JARVIS farà il resto!\n\n\n"
# 			source ~/Scrivania/SCRIPT_PIPELINE/Cardio_function.sh
# 			CARDIO_ALL
# 			CARDIO_PREPROCESSING
# 			CARDIO_VARIANT_CALLING_AND_ANNOTATION
# 			fi
		
		
# 			if [ "$Analisi" == "P" ] || [ "$Analisi" == "p" ]
# 			then
# 			printf "\n\nNon verrà effettuato l'allineamento. L'analisi avrà inizio a partire dal pre-processing.\n\n"
# 			#Controllo la presenza di file bam
# 			cd $PROCESSING/1_Align/Sort_bam
# 			count=$(ls -1 *.bam 2>/dev/null | wc -l)
# 			echo -n "Quanti sample stai analizzando? "
# 			read Sample
# 				if [[ $count != "$Sample" ]];
# 				then 
# 				printf "\n\nIl numero di bam non corrisponde al numero di pazienti che hai inserito. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n\nArrivederci!\n\n"
# 				exit 1;
# 				else
# 				printf "\n\nIl numero di bam è uguale al numero di pazienti. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n"
# 				ls *.bai
# 				printf "\n\n"
# 				echo -n "Confermi che i file sono corretti (S/N)? "
# 				read Answer
# 					until [ "$Answer" == "S" ] || [ "$Answer" == "s" ] || [ "$Answer" == "N" ] || [ "$Answer" == "n" ]
# 					do
# 					printf "\nInserire Sì (S/s) o No (N/n): "
# 					read Answer
# 					done
# 						if [ "$Answer" == "S" ] || [ "$Answer" == "s" ]
# 						then
# 						printf "\nL'analisi ha inizio! Mettiti comodo, JARVIS farà il resto!\n\n"
# 						else
# 						printf "\nCarica i file corretti nel path $PROCESSING/1_Align/Sort_bam e poi torna da me!\n\n"
# 						exit 1;
# 						fi
# 				fi
# 			source ~/Scrivania/SCRIPT_PIPELINE/Cardio_function.sh
# 			CARDIO_PREPROCESSING
# 			CARDIO_VARIANT_CALLING_AND_ANNOTATION
# 			fi


# 			if [ "$Analisi" == "V" ] || [ "$Analisi" == "v" ]
# 			then
# 			printf "\n\nL'analisi avrà inizio a partire dal Variant calling.\n\n"
# 			cd $PROCESSING/5_BQSR/
# 			#Controllo che esistano file bam da cui partire con l'analisi
# 			count=$(ls -1 *.bam 2>/dev/null | wc -l)
# 			echo -n "Quanti sample stai analizzando? "
# 			read Sample
# 				if [[ $count != "$Sample" ]];
# 				then 
# 				printf "\n\nIl numero di bam non corrisponde al numero di pazienti che hai inserito. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n\nArrivederci!\n\n"
# 				exit 1;
# 				else
# 				printf "\n\nIl numero di bam è uguale al numero di pazienti. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n"
# 				ls *.bai
# 				printf "\n\n"
# 				echo -n "Confermi che i file sono corretti (S/N)? Ricorda che devono essere presenti anche i file .bai per ciascun .bam : "
# 				read Answer
# 					until [ "$Answer" == "S" ] || [ "$Answer" == "s" ] || [ "$Answer" == "N" ] || [ "$Answer" == "n" ]
# 					do
# 					printf "\nInserire Sì (S/s) o No (N/n): "
# 					read Answer
# 					done
# 						if [ "$Answer" == "S" ] || [ "$Answer" == "s" ]
# 						then
# 						printf "\nL'analisi ha inizio! Mettiti comodo, JARVIS farà il resto!\n\n"
# 						else
# 						printf "\nCarica i file corretti nel path $PROCESSING/5_BQSR e poi torna da me!\n\n"
# 						exit 1;
# 						fi
# 				fi

# 			source ~/Scrivania/SCRIPT_PIPELINE/Cardio_function.sh
# 			CARDIO_VARIANT_CALLING_AND_ANNOTATION
# 			fi


# 	#-----------------------------------------------------------------BRCA-----------------------------------------------------------------


# 	elif [ "$PANNELLO" == "B" ] || [ "$PANNELLO" == "b" ]
# 	then

# 	cd $OUTVCF
# 	Name_OUT=$DataRun\_Run_$NumRun\_BRCA_Conferme_$PipeVersion
# 	mkdir "$Name_OUT"
# 	cd $STORAGE
# 	Name_Dir=$DataRun\_Run_$NumRun\_BRCA_$PipeVersion
# 	mkdir "$Name_Dir"

# 		printf "\n\n"
# 		echo -n "Vuoi partire dall'allineamento (A/a), dal pre-processing (P/p) o direttamente dal variant-calling (V/v) ? ==> N.B Se si decide di saltare l'allineamento, verificare prima che i files .bam allineati siano contenuti nella cartella $PROCESSING/1_Align/Sort_bam con estensione DataRun_n°paziente_tipo-analisi (es: ${DataRun}_01_Cardio). I file da analizzare dovranno essere già ordinati con il comando SortSam del Picardtool. Qualora invece si parta dal variant Calling, verificare che i file siano presenti nella cartella $PROCESSING/5_BQSR. In questo caso i file utilizzati devono già essere stati sottoposti a pre-processing! [Enter/Invio]: "

# 		read Analisi

# 		until [ "$Analisi" == "A" ] || [ "$Analisi" == "a" ] || [ "$Analisi" == "P" ] || [ "$Analisi" == "p" ] || [ "$Analisi" == "V" ] || [ "$Analisi" == "v" ]
# 		do
# 		printf "\nInserire (A/a) per partire dall'allineamento, (P/p) per partire dal pre-processing oppure (V/v) per iniziare con il variant calling [Enter/Invio]: "
# 		read Analisi
# 		done


# 			if [ "$Analisi" == "A" ] || [ "$Analisi" == "a" ]
# 			then
# 			printf "\n\nL'analisi avrà inizio a partire dall'allineamento.  Mettiti comodo, JARVIS farà il resto!\n\n\n"
# 			source ~/Scrivania/SCRIPT_PIPELINE/BRCA_function.sh
# 			BRCA_ALL
# 			BRCA_PREPROCESSING
# 			BRCA_VARIANT_CALLING_AND_ANNOTATION
# 			fi
			

# 			if [ "$Analisi" == "P" ] || [ "$Analisi" == "p" ]
# 			then
# 			printf "\n\nNon verrà effettuato l'allineamento. L'analisi avrà inizio a partire dal pre-processing.\n\n"
# 			#Controllo la presenza di file bam
# 			cd $PROCESSING/1_Align/Sort_bam
# 			count=$(ls -1 *.bam 2>/dev/null | wc -l)
# 			echo -n "Quanti sample stai analizzando? "
# 			read Sample
# 				if [[ $count != "$Sample" ]];
# 				then 
# 				printf "\n\nIl numero di bam non corrisponde al numero di pazienti che hai inserito. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n\nArrivederci!\n\n"
# 				exit 1;
# 				else
# 				printf "\n\nIl numero di bam è uguale al numero di pazienti. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n"
# 				ls *.bai
# 				printf "\n\n"
# 				echo -n "Confermi che i file sono corretti (S/N)? "
# 				read Answer
# 					until [ "$Answer" == "S" ] || [ "$Answer" == "s" ] || [ "$Answer" == "N" ] || [ "$Answer" == "n" ]
# 					do
# 					printf "\nInserire Sì (S/s) o No (N/n): "
# 					read Answer
# 					done
# 						if [ "$Answer" == "S" ] || [ "$Answer" == "s" ]
# 						then
# 						printf "\nL'analisi ha inizio! Mettiti comodo, JARVIS farà il resto!\n\n"
# 						else
# 						printf "\nCarica i file corretti nel path $PROCESSING/1_Align/Sort_bam e poi torna da me!\n\n"
# 						exit 1;
# 						fi
# 				fi
# 			source ~/Scrivania/SCRIPT_PIPELINE/BRCA_function.sh
# 			BRCA_PREPROCESSING
# 			BRCA_VARIANT_CALLING_AND_ANNOTATION
# 			fi


# 			if [ "$Analisi" == "V" ] || [ "$Analisi" == "v" ]
# 			then
# 			printf "\n\nL'analisi avrà inizio a partire dal variant calling.\n\n"
# 			#Controllo la presenza di file bam
# 			cd $PROCESSING/4_Indel/
# 			count=$(ls -1 *.bam 2>/dev/null | wc -l)
# 			echo -n "Quanti sample stai analizzando? "
# 			read Sample
# 				if [[ $count != "$Sample" ]];
# 				then 
# 				printf "\n\nIl numero di bam non corrisponde al numero di pazienti che hai inserito. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n\nArrivederci!\n\n"
# 				exit 1;
# 				else
# 				printf "\n\nIl numero di bam è uguale al numero di pazienti. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n"
# 				ls *.bai
# 				printf "\n\n"
# 				echo -n "Confermi che i file sono corretti (S/N)? Ricorda che devono essere presenti anche i file .bai per ciascun .bam : "
# 				read Answer
# 					until [ "$Answer" == "S" ] || [ "$Answer" == "s" ] || [ "$Answer" == "N" ] || [ "$Answer" == "n" ]
# 					do
# 					printf "\nInserire Sì (S/s) o No (N/n): "
# 					read Answer
# 					done
# 						if [ "$Answer" == "S" ] || [ "$Answer" == "s" ]
# 						then
# 						printf "\nL'analisi ha inizio! Mettiti comodo, JARVIS farà il resto!\n\n"
# 						else
# 						printf "\nCarica i file corretti nel path $PROCESSING/4_Indel e poi torna da me!\n\n"
# 						exit 1;
# 						fi
# 				fi
# 			source ~/Scrivania/SCRIPT_PIPELINE/BRCA_function.sh
# 			BRCA_VARIANT_CALLING_AND_ANNOTATION
# 			fi


# 	#---------------------------------------------------------------CANCER-----------------------------------------------------------------


# 	elif [ "$PANNELLO" == "Cancer" ] || [ "$PANNELLO" == "cancer" ]
# 	then

# 	cd $OUTVCF
# 	Name_OUT=$DataRun\_Run_$NumRun\_Cancer_Conferme_$PipeVersion
# 	mkdir "$Name_OUT"
# 	cd $STORAGE
# 	Name_Dir=$DataRun\_Run_$NumRun\_Cancer_$PipeVersion
# 	mkdir "$Name_Dir"

# 		printf "\n\n"
# 		echo -n "Vuoi partire dall'allineamento (A/a), dal pre-processing (P/p) o direttamente dal variant-calling (V/v) ? ==> N.B Se si decide di saltare l'allineamento, verificare prima che i files .bam allineati siano contenuti nella cartella $PROCESSING/1_Align/Sort_bam con estensione DataRun_n°paziente_tipo-analisi (es: ${DataRun}_01_Cardio). I file da analizzare dovranno essere già ordinati con il comando SortSam del Picardtool. Qualora invece si parta dal variant Calling, verificare che i file siano presenti nella cartella $PROCESSING/5_BQSR. In questo caso i file utilizzati devono già essere stati sottoposti a pre-processing! [Enter/Invio]: "

# 		read Analisi

# 		until [ "$Analisi" == "A" ] || [ "$Analisi" == "a" ] || [ "$Analisi" == "P" ] || [ "$Analisi" == "p" ] || [ "$Analisi" == "V" ] || [ "$Analisi" == "v" ]
# 		do
# 		printf "\nInserire (A/a) per partire dall'allineamento, (P/p) per partire dal pre-processing oppure (V/v) per iniziare con il variant calling [Enter/Invio]: "
# 		read Analisi
# 		done
# 		printf "$Analisi"


# 			if [ "$Analisi" == "A" ] || [ "$Analisi" == "a" ]
# 			then
# 			printf "\n\nL'analisi avrà inizio a partire dall'allineamento. Mettiti comodo, JARVIS farà il resto!\n\n\n"
# 			source ~/Scrivania/SCRIPT_PIPELINE/Cancer_function.sh
# 			CANCER_ALL
# 			CANCER_PREPROCESSING
# 			CANCER_VARIANT_CALLING_AND_ANNOTATION
# 			fi
		
		
# 			if [ "$Analisi" == "P" ] || [ "$Analisi" == "p" ]
# 			then
# 			printf "\n\nNon verrà effettuato l'allineamento. L'analisi avrà inizio a partire dal pre-processing.\n\n"
# 			#Controllo la presenza di file bam
# 			cd $PROCESSING/1_Align/Sort_bam
# 			count=$(ls -1 *.bam 2>/dev/null | wc -l)
# 			echo -n "Quanti sample stai analizzando? "
# 			read Sample
# 				if [[ $count != "$Sample" ]];
# 				then 
# 				printf "\n\nIl numero di bam non corrisponde al numero di pazienti che hai inserito. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n\nArrivederci!\n\n"
# 				exit 1;
# 				else
# 				printf "\n\nIl numero di bam è uguale al numero di pazienti. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n"
# 				ls *.bai
# 				printf "\n\n"
# 				echo -n "Confermi che i file sono corretti (S/N)? "
# 				read Answer
# 					until [ "$Answer" == "S" ] || [ "$Answer" == "s" ] || [ "$Answer" == "N" ] || [ "$Answer" == "n" ]
# 					do
# 					printf "\nInserire Sì (S/s) o No (N/n): "
# 					read Answer
# 					done
# 						if [ "$Answer" == "S" ] || [ "$Answer" == "s" ]
# 						then
# 						printf "\nL'analisi ha inizio! Mettiti comodo, JARVIS farà il resto!\n\n"
# 						else
# 						printf "\nCarica i file corretti nel path $PROCESSING/1_Align/Sort_bam e poi torna da me!\n\n"
# 						exit 1;
# 						fi
# 				fi
# 			source ~/Scrivania/SCRIPT_PIPELINE/Cancer_function.sh
# 			CANCER_PREPROCESSING
# 			CANCER_VARIANT_CALLING_AND_ANNOTATION
# 			fi


# 			if [ "$Analisi" == "V" ] || [ "$Analisi" == "v" ]
# 			then
# 			printf "\n\nL'analisi avrà inizio a partire dal Variant calling.\n\n"
# 			cd $PROCESSING/5_BQSR/
# 			#Controllo che esistano file bam da cui partire con l'analisi
# 			count=$(ls -1 *.bam 2>/dev/null | wc -l)
# 			echo -n "Quanti sample stai analizzando? "
# 			read Sample
# 				if [[ $count != "$Sample" ]];
# 				then 
# 				printf "\n\nIl numero di bam non corrisponde al numero di pazienti che hai inserito. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n\nArrivederci!\n\n"
# 				exit 1;
# 				else
# 				printf "\n\nIl numero di bam è uguale al numero di pazienti. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n"
# 				ls *.bai
# 				printf "\n\n"
# 				echo -n "Confermi che i file sono corretti (S/N)? Ricorda che devono essere presenti anche i file .bai per ciascun .bam : "
# 				read Answer
# 					until [ "$Answer" == "S" ] || [ "$Answer" == "s" ] || [ "$Answer" == "N" ] || [ "$Answer" == "n" ]
# 					do
# 					printf "\nInserire Sì (S/s) o No (N/n): "
# 					read Answer
# 					done
# 						if [ "$Answer" == "S" ] || [ "$Answer" == "s" ]
# 						then
# 						printf "\nL'analisi ha inizio! Mettiti comodo, JARVIS farà il resto!\n\n"
# 						else
# 						printf "\nCarica i file corretti nel path $PROCESSING/5_BQSR e poi torna da me!\n\n"
# 						exit 1;
# 						fi
# 				fi

# 			source ~/Scrivania/SCRIPT_PIPELINE/Cancer_function.sh
# 			CANCER_VARIANT_CALLING_AND_ANNOTATION
# 			fi


# 	#---------------------------------------------------------------EXOME-----------------------------------------------------------------


# 	elif [ "$PANNELLO" == "E" ] || [ "$PANNELLO" == "e" ]
# 	then

# 	cd $OUTVCF
# 	Name_OUT=$DataRun\_Run_$NumRun\_Exome_Conferme_$PipeVersion
# 	mkdir "$Name_OUT"
# 	cd $STORAGE
# 	Name_Dir=$DataRun\_Run_$NumRun\_Exome_$PipeVersion
# 	mkdir "$Name_Dir"

# 		printf "\n\n"
# 		echo -n "Vuoi partire dall'allineamento (A/a), dal pre-processing (P/p) o direttamente dal variant-calling (V/v) ? ==> N.B Se si decide di saltare l'allineamento, verificare prima che i files .bam allineati siano contenuti nella cartella $PROCESSING/1_Align/Sort_bam con estensione DataRun_n°paziente_tipo-analisi (es: ${DataRun}_01_Cardio). I file da analizzare dovranno essere già ordinati con il comando SortSam del Picardtool. Qualora invece si parta dal variant Calling, verificare che i file siano presenti nella cartella $PROCESSING/5_BQSR. In questo caso i file utilizzati devono già essere stati sottoposti a pre-processing! [Enter/Invio]: "

# 		read Analisi

# 		until [ "$Analisi" == "A" ] || [ "$Analisi" == "a" ] || [ "$Analisi" == "P" ] || [ "$Analisi" == "p" ] || [ "$Analisi" == "V" ] || [ "$Analisi" == "v" ]
# 		do
# 		printf "\nInserire (A/a) per partire dall'allineamento, (P/p) per partire dal pre-processing oppure (V/v) per iniziare con il variant calling [Enter/Invio]: "
# 		read Analisi
# 		done
# 		printf "$Analisi"


# 			if [ "$Analisi" == "A" ] || [ "$Analisi" == "a" ]
# 			then
# 			printf "\n\nL'analisi avrà inizio a partire dall'allineamento. Mettiti comodo, JARVIS farà il resto!\n\n\n"
# 			source ~/Scrivania/SCRIPT_PIPELINE/Exome_function.sh
# 			EXOME_ALL
# 			EXOME_PREPROCESSING
# 			EXOME_VARIANT_CALLING_AND_ANNOTATION
# 			fi
		
		
# 			if [ "$Analisi" == "P" ] || [ "$Analisi" == "p" ]
# 			then
# 			printf "\n\nNon verrà effettuato l'allineamento. L'analisi avrà inizio a partire dal pre-processing.\n\n"
# 			#Controllo la presenza di file bam
# 			cd $PROCESSING/1_Align/Sort_bam
# 			count=$(ls -1 *.bam 2>/dev/null | wc -l)
# 			echo -n "Quanti sample stai analizzando? "
# 			read Sample
# 				if [[ $count != "$Sample" ]];
# 				then 
# 				printf "\n\nIl numero di bam non corrisponde al numero di pazienti che hai inserito. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n\nArrivederci!\n\n"
# 				exit 1;
# 				else
# 				printf "\n\nIl numero di bam è uguale al numero di pazienti. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n"
# 				ls *.bai
# 				printf "\n\n"
# 				echo -n "Confermi che i file sono corretti (S/N)? "
# 				read Answer
# 					until [ "$Answer" == "S" ] || [ "$Answer" == "s" ] || [ "$Answer" == "N" ] || [ "$Answer" == "n" ]
# 					do
# 					printf "\nInserire Sì (S/s) o No (N/n): "
# 					read Answer
# 					done
# 						if [ "$Answer" == "S" ] || [ "$Answer" == "s" ]
# 						then
# 						printf "\nL'analisi ha inizio! Mettiti comodo, JARVIS farà il resto!\n\n"
# 						else
# 						printf "\nCarica i file corretti nel path $PROCESSING/1_Align/Sort_bam e poi torna da me!\n\n"
# 						exit 1;
# 						fi
# 				fi
# 			source ~/Scrivania/SCRIPT_PIPELINE/Exome_function.sh
# 			EXOME_PREPROCESSING
# 			EXOME_VARIANT_CALLING_AND_ANNOTATION
# 			fi


# 			if [ "$Analisi" == "V" ] || [ "$Analisi" == "v" ]
# 			then
# 			printf "\n\nL'analisi avrà inizio a partire dal Variant calling.\n\n"
# 			cd $PROCESSING/5_BQSR/
# 			#Controllo che esistano file bam da cui partire con l'analisi
# 			count=$(ls -1 *.bam 2>/dev/null | wc -l)
# 			echo -n "Quanti sample stai analizzando? "
# 			read Sample
# 				if [[ $count != "$Sample" ]];
# 				then 
# 				printf "\n\nIl numero di bam non corrisponde al numero di pazienti che hai inserito. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n\nArrivederci!\n\n"
# 				exit 1;
# 				else
# 				printf "\n\nIl numero di bam è uguale al numero di pazienti. I file sono i seguenti:\n\n"
# 				ls *.bam
# 				printf "\n"
# 				ls *.bai
# 				printf "\n\n"
# 				echo -n "Confermi che i file sono corretti (S/N)? Ricorda che devono essere presenti anche i file .bai per ciascun .bam : "
# 				read Answer
# 					until [ "$Answer" == "S" ] || [ "$Answer" == "s" ] || [ "$Answer" == "N" ] || [ "$Answer" == "n" ]
# 					do
# 					printf "\nInserire Sì (S/s) o No (N/n): "
# 					read Answer
# 					done
# 						if [ "$Answer" == "S" ] || [ "$Answer" == "s" ]
# 						then
# 						printf "\nL'analisi ha inizio! Mettiti comodo, JARVIS farà il resto!\n\n"
# 						else
# 						printf "\nCarica i file corretti nel path $PROCESSING/5_BQSR e poi torna da me!\n\n"
# 						exit 1;
# 						fi
# 				fi

# 			source ~/Scrivania/SCRIPT_PIPELINE/Exome_function.sh
# 			EXOME_VARIANT_CALLING_AND_ANNOTATION
# 			fi


# 	else
# 	printf "\n\n"
# 	echo "Error! Non esiste alcuna analisi del tipo $PANNELLO. Arrivederci!"
# 	printf "\n\n"
# 	exit
# 	fi



# 	ENDTIME=$(date +%s)
# 	echo "TEMPO TOTALE DI PROCESSING: $(($ENDTIME - $STARTTIME)) sec"
# 	printf "\n\n"


