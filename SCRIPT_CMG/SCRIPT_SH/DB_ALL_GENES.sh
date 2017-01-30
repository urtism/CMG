#!/bin/bash


while [[ $# -gt 0 ]]
do
key="$1"
case $key in
	-p|--path)
	PATH="$2"
	shift # past argument
	;;
	-g|--genelist)
	GENELIST="$2"
	shift # past argument
	;;
	-t|--tipo)
	TIPO="$2"
	shift # past argument
	;;
	-l|--paz_list)
	PAZ_LIST="$2"
	shift # past argument
	;;
	-o|--out)
	OUT="$2"
	shift # past argument
	;;
	--default)
	DEFAULT=1
	;;
	*)
	echo "ERROR: Wrong command(s). To read the HELP type [-h] [--help] option."
	exit 1;
	;;
esac
shift 
done

# PATH=/home/jarvis/Scrivania/DB_GENI_VAR/CANCER_20170127/20170127_Run_ALL_Cancer_Conferme_2.1
# GENELIST=/home/jarvis/NGS_ANALYSIS/TARGET/gene_list_trusight_cancer.txt
# TIPO=Cancer
# PAZ_LIST=/home/jarvis/Scrivania/Lista_pazienti_20170130.txt
# OUT=/home/jarvis/Scrivania/DB_GENI_VAR/CANCER_20170127

while read -r line
do

GENE=$line

echo -e "/usr/bin/python ~/git/CMG/SCRIPT_CMG/SCRIPT_PYTHON/DB_per_Gene.py --path $PATH --gene $GENE --out $OUT/$GENE --paz_list $PAZ_LIST --tipo $TIPO"
/usr/bin/python ~/git/CMG/SCRIPT_CMG/SCRIPT_PYTHON/DB_per_Gene.py --path $PATH --gene $GENE --out $OUT/$GENE --paz_list $PAZ_LIST --tipo $TIPO

done < "$GENELIST"

