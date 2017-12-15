#!bin/bash


if [ "$#" == "0" ]
then
	echo "\nERROR: No INPUT command(s). To read the HELP type [-h] [--help] option.\n"
	exit 1;
else
	while [[ $# -gt 0 ]]
	do
	key="$1"
	case $key in
		-d|--db_path)
		DB_PATH="$2"
		shift # past argument
		;;
		-g|--gene_list)
		GENELIST="$2"
		shift # past argument
		;;
		--default)
		DEFAULT=1
		;;
		*)
		exit 1;
		;;
	esac
	shift 
	done
fi

mkdir $DB_PATH/MERGE/
mkdir $DB_PATH/MERGE/GENI

while read -r line
do

GENE=${line%$'\r'}
echo $GENE

mkdir $DB_PATH/MERGE/GENI/$GENE/

python ~/git/CMG/SCRIPT_CMG/SCRIPT_PYTHON/unisci_Cardio_Conn_per_gene.py --cardio $DB_PATH/CARDIO/GENI/$GENE/lista_varianti_$GENE\_STATS.tsv --conn $DB_PATH/CONN/GENI/$GENE/lista_varianti_$GENE\_STATS.tsv --gene $GENE --out $DB_PATH/MERGE/GENI/$GENE/lista_varianti_$GENE\_STATS.tsv


done < "$GENELIST"
