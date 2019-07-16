#!bin/bash

###FUNZIONI###


HaploTypeCaller() {
	echo "java -Xmx64g -jar $GATK -T HaplotypeCaller \
	-R $REF \
	-I $1 \
	-o $2 \
	-ERC GVCF \
	--doNotRunPhysicalPhasing \
	-L $TARGET"

	java -Xmx64g -jar $GATK -T HaplotypeCaller \
	-R $REF \
	-I $1 \
	-o $2 \
	-ERC GVCF \
	--doNotRunPhysicalPhasing \
	-L $TARGET
}

eliminaBRCA_Cancer() {
	cd $GVCF_PATH
	while read -r line
	do

	PAZ=$line

	rm $PAZ
	rm $PAZ.idx

	done < "/home/jarvis/Scrivania/DB_GENI_VAR/lista_BRCA_che_hanno_fatto_Cancer.list"
	cd 
}

scarica_gvcf() {
	TIPO=$2

	cd $GVCF_PATH

	cat $PAZ_LIST | while read line
		do
			NAME=$(echo "$line" | cut -f1)
			NUM_RUN=$(echo "$line" | cut -f3)
			DATA_RUN=$(echo "$line" | cut -f4)
			SAMPLE_ID=$(echo "$line" | cut -f5)
			CHECK=$(echo "$line" | cut -f7)

			gvcf=$SAMPLE_ID.g.vcf
			if [ -e $gvcf ] 
			then
				#echo "$gvcf presente: check $CHECK"
				if [ "$CHECK" == "0" ] 
				then
					rm $gvcf
				else
					printf $"$GVCF_PATH/$gvcf\n" >> $gvcflist
				fi;
			else 
				if [[ "$SAMPLE_ID" == *"$TIPO"* ]] && [ "$CHECK" == "1" ] 
				then
					
					if [[ $((gsutil -m cp gs://storage_run_cmgcv/$1/$DATA_RUN\_Run_$NUM_RUN\_$1\_2.1/$gvcf .) 2>&1 | grep "Exception" ) ]]; then
						

						if  [[ "$TIPO" == "Conn" ]] || [[ "$TIPO" == "Cardio" ]]
						then

							dir=/media/jarvis/HD2/Storage_Run_cmgcv/CARDIO/$DATA_RUN\_Run_$NUM_RUN\_$1\_2.1
							bam=$dir/$SAMPLE_ID.bam

							if test -d $dir
								then
									bam=$dir/$SAMPLE_ID.bam
									if [ -f $bam ]
									then
										HaploTypeCaller $bam $GVCF_PATH/$SAMPLE_ID.g.vcf
										printf $"$GVCF_PATH/$gvcf\n" >> $gvcflist
									fi
							else
								dir=/media/jarvis/IANSOLO/JARVIS-STORAGE/Storage_Run_Hard_Disk_Esterno/CARDIO/$DATA_RUN\_Run_$NUM_RUN\_$1\_2.1
								bam=$dir/$SAMPLE_ID.bam
									if [ -f $bam ]
										then
										HaploTypeCaller $bam $GVCF_PATH/$SAMPLE_ID.g.vcf
										printf $"$GVCF_PATH/$gvcf\n" >> $gvcflist
									fi
							fi

						elif [[ "$TIPO" == "Cancer" ]]
							then
							dir=/media/jarvis/IANSOLO/JARVIS-STORAGE/Storage_Run_Hard_Disk_Esterno/CANCER/$DATA_RUN\_Run_$NUM_RUN\_$1\_2.1
							bam=$dir/$SAMPLE_ID.bam
							if [ -f $bam ]
							then
								HaploTypeCaller $bam $GVCF_PATH/$SAMPLE_ID.g.vcf
								printf $"$GVCF_PATH/$gvcf\n" >> $gvcflist
							fi


						fi


					fi
				fi
			fi
		done

	cd
}




### TOOLS ###
SCRIPT_PIPELINE=~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE
FASTQC=~/NGS_TOOLS/FastQC/fastqc
BWA=~/NGS_TOOLS/bwa-0.7.15
PICARD=~/NGS_TOOLS/picard-tools-2.7.1/picard.jar

GATK=~/NGS_TOOLS/GATK/GenomeAnalysisTK.jar
VARSCAN=~/NGS_TOOLS/VarScan/VarScan.v2.3.9.jar
FREEBAYES=~/NGS_TOOLS/freebayes/bin/freebayes
VARDICT=~/NGS_TOOLS/VarDictJava-master/build/install/VarDict/bin/VarDict
BCFTOOLS=bcftools

VEP=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/
VEPANN=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEPFILTER=~/NGS_TOOLS/ensembl-tools-release-86/scripts/variant_effect_predictor/filter_vep.pl

### DATABASES & FILES ###
LISTAFEATURES_GERMLINE=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_germline.list
LISTAFEATURES_SOMATIC=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_somatic.list
#LISTAFEATURES_SOMATIC=/home/jarvis/NGS_ANALYSIS/TARGET/Features_lists/lista_features_somatic_CF_20170315.list
ANN_LIST_GERMLINE=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_annotazione.list
ANN_LIST_SOMATIC=~/NGS_ANALYSIS/TARGET/Features_lists/lista_features_annotazione.list
REF=~/NGS_TOOLS/hg19/ucsc.hg19.fasta
MILLS=~/NGS_TOOLS/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DBSNP=~/NGS_TOOLS/hg19/dbsnp_138.hg19.vcf
LOGHI=~/git/CMG/LOGHI
TRASCR_CARDIO=~/NGS_ANALYSIS/TARGET/Lista_trascritti_Cardio.txt
TRASCR_BRCA=~/NGS_ANALYSIS/TARGET/Lista_trascritti_BRCA.txt
TRASCR_CANCER=~/NGS_ANALYSIS/TARGET/Lista_trascritti_Cancer.txt

### TARGET ###
#TARGET_CARDIO_1000=~/Scrivania/TEST/col3a1/COL3A1.list
#TARGET_CARDIO_1000_BED=~/Scrivania/TEST/col3a1/COL3A1.bed
TARGET_CARDIO_1000=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.list
TARGET_CARDIO_1000_BED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-1000.bed
TARGET_CARDIO_100=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-100.list
TARGET_CARDIO_100_BED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-100.bed
TARGET_CARDIO_25=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-25.list
TARGET_CARDIO_25_BED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a_ESTESO+-25.bed
TARGET_CARDIO=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a.list
TARGET_CARDIO_BED=~/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a.bed
TARGET_BRCA=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.list
TARGET_BRCA_BED=~/NGS_ANALYSIS/TARGET/AFP2_manifest_v1.bed
TARGET_BRCA_FREEBAYES=~/NGS_ANALYSIS/TARGET/BRCA_FreeBayes_amplicon.bed
TARGET_EXOME_1000=~/NGS_ANALYSIS/TARGET/TruSight_One_v1.1_ESTESO+-1000.list
TARGET_EXOME_1000_BED=~/NGS_ANALYSIS/TARGET/TruSight_One_v1.1_ESTESO+-1000.bed
TARGET_CF=~/NGS_ANALYSIS/TARGET/ctDNA_2_113416_AmpliconsExport.list
TARGET_CF_BED=~/NGS_ANALYSIS/TARGET/ctDNA_2_113416_AmpliconsExport.bed
TARGET_CANCER_1000=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.list
TARGET_CANCER_1000_BED=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-1000.bed
TARGET_CANCER_25=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-25.list
TARGET_CANCER_25_BED=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-25.bed
TARGET_CANCER_100=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-100.list
TARGET_CANCER_100_BED=~/NGS_ANALYSIS/TARGET/trusight_cancer_manifest_a_ESTESO+-100.bed

if [ "$#" == "0" ]
then
	echo "\nERROR: No INPUT command(s). To read the HELP type [-h] [--help] option.\n"
	exit 1;
else
	while [[ $# -gt 0 ]]
	do
	key="$1"
	case $key in
		--gvcf_path)
		GVCF_PATH="$2"
		shift # past argument
		;;
		--gvcf_merged)
		MGVCF="$2"
		shift # past argument
		;;
		-o|--out_path)
		OUT_PATH="$2"
		shift # past argument
		;;
		-o|--data)
		DATA="$2"
		shift # past argument
		;;
		-p|--pannello)
		PANNELLO="$2"
		shift # past argument
		;;
		-g|--gene_list)
		GENELIST="$2"
		shift # past argument
		;;
		-l|--paz_list)
		PAZ_LIST="$2"
		shift # past argument
		;;
		-h|--help)
		HELP
		exit 1;
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
fi

mkdir $OUT_PATH
mkdir $OUT_PATH/VCF
mkdir $OUT_PATH/TSV
mkdir $OUT_PATH/GENI

gvcflist=$OUT_PATH/gvcf.list

if [ "$PANNELLO" == "Cardio" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CARDIO_100
		TARGETBED=$TARGET_CARDIO_100_BED
		TRANSCR_LIST=$TRASCR_CARDIO
		tipo=Cardio-Conn
		echo "Aggiungo i gvcf alla lista"
		scarica_gvcf $PANNELLO Cardio
		scarica_gvcf $PANNELLO Conn

		
	elif [ "$PANNELLO" == "Cancer" ]
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_CANCER_100
		TARGETBED=$TARGET_CANCER_100_BED
		TRANSCR_LIST=$TRASCR_CANCER
		tipo=Cancer
		echo "Aggiungo i gvcf alla lista"
		scarica_gvcf $PANNELLO $tipo 
		
	elif [ "$PANNELLO" == "Exome" ] 
	then
		DESIGN="ENRICHMENT"
		TARGET=$TARGET_EXOME_1000
		TARGETBED=$TARGET_EXOME_1000_BED
		
	elif [ "$PANNELLO" == "BRCA" ]
	then
		
		eliminaBRCA_Cancer
		DESIGN="AMPLICON"
		TARGET=$TARGET_BRCA
		TARGETBED=$TARGET_BRCA_BED
		TRANSCR_LIST=$TRASCR_BRCA
fi

##ls $GVCF_PATH/*.g.vcf > $gvcflist

mkdir $GVCF_PATH/TMP
echo "Faccio il merge dei gvcf"

python /home/jarvis/git/CMG/SCRIPT_CMG/SCRIPT_PYTHON/DATABASE_SCRIPTING/merge_gvcf.py \
-R /home/jarvis/NGS_TOOLS/hg19/ucsc.hg19.fasta \
-V $gvcflist \
-L $TARGET \
-g 10 \
-O $GVCF_PATH/TMP/$tipo.merge.g.vcf

MGVCF=$GVCF_PATH/TMP/$tipo.merge.g.vcf
#MGVCF="/media/jarvis/HD1/DB.gvcf/CARDIO/TMP/Cardio-Conn.merge.g.vcf"
#MGVCF="/media/jarvis/HD1/DB.gvcf/CANCER/TMP/Cancer.merge.g.vcf"

echo "Faccio il Variant Calling"

java -jar -Xmx45g $GATK -T GenotypeGVCFs \
-R $REF \
-V:VCF $MGVCF \
-o $OUT_PATH/$DATA\_$PANNELLO.vcf
#--disable_auto_index_creation_and_locking_when_reading_rods

java -jar $GATK -T VariantFiltration \
-R $REF \
-V $OUT_PATH/$DATA\_$PANNELLO.vcf \
--filterExpression "QD < 2.0 || FS > 100.0 || ReadPosRankSum < -16.0 || DP < 15" \
--filterName "FILTER" \
-o $OUT_PATH/$DATA\_$PANNELLO.FILTER.vcf

python $SCRIPT_PIPELINE/header_fix.py -v G -f $OUT_PATH/$DATA\_$PANNELLO.FILTER.vcf \
> $OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.vcf

$BCFTOOLS norm -m -both \
-f $REF \
$OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.vcf \
> $OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.vcf


echo "Annoto le varianti"

/home/jarvis/NGS_TOOLS/ensembl-vep-release-95/vep -i $OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.vcf \
-o $OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.ANN.vcf \
--stats_file ${1%.*}.ANN.html \
--fasta $REF \
--cache \
--format vcf \
--assembly GRCh37 \
--merged \
--offline \
--force_overwrite \
-v \
--fork 4 \
--variant_class \
--sift b \
--poly b \
--vcf_info_field ANN \
--hgvs \
--protein \
--check_existing \
--af \
--species homo_sapiens \
--failed 1 \
--vcf

sed -n -e '/#/p' $OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.ANN.vcf > \
$OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.ANN.SPLIT.vcf

perl -e 'while(<>) { chomp; if(m/(.+?)(CSQ|ANN)\=([^;^\s]+)(.*)/) { foreach my $s(split ",", $3) { print "$1$2\=$s\;$4\n"}}}' \
$OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.ANN.vcf >> $OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.ANN.SPLIT.vcf


grep -e "#" -f $TRANSCR_LIST $OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.ANN.SPLIT.vcf > $OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.ANN.SPLIT.TRANSCR.vcf

sed -n -e '/#CHROM/p' $OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.ANN.SPLIT.TRANSCR.vcf > $OUT_PATH/List_samples_for_split.txt
sed -i "s/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t//g" $OUT_PATH/List_samples_for_split.txt

COUNT=$(awk '{print NF}' $OUT_PATH/List_samples_for_split.txt | sort -nu | tail -n 1)

#COUNT=$(awk '{print NF}' /home/jarvis/Scrivania/DB_GENI_VAR/CARDIO_20190319/List_samples_for_split.txt | sort -nu | tail -n 1)

sed -i -e "s/\(AN=[[:digit:]]*\);\(DP=[[:digit:]]*\)/\1;;;\2/g" \
-e "s/FS=\([[:digit:]]*\|[[:digit:]]*.[[:digit:]]*\);MLEAC=\([[:digit:]]*\|[[:digit:]]*.[[:digit:]]\)/FS=\1;;MLEAC=\2/g" \
-e "s/MQ=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\);QD=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\)/MQ=\1;;QD=\2/g" \
-e "s/QD=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\);SOR=\([[:digit:]]*\|[0-9]*\|[0-9]*.[0-9]*\)/QD=\1;;SOR=\2/g" $OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.ANN.SPLIT.TRANSCR.vcf

echo "Estraggo le varianti per ciascun gene"

for (( a=1; a<$COUNT+1; a++ ))
	do
	NAME=$(cut -f $a $OUT_PATH/List_samples_for_split.txt)
	b=$((9+$a))
	cut -f1-9,$b $OUT_PATH/$DATA\_$PANNELLO.FILTER.FIX.NORM.ANN.SPLIT.TRANSCR.vcf > $OUT_PATH/VCF/$NAME.vcf
	sed -i -e '/0\/0/d' -e '/\.\/\./d' -e 's/|;/|-;/g' -e 's/;\t/\t/g' $OUT_PATH/VCF/$NAME.vcf
	sed -i -e '/^##/d' $OUT_PATH/VCF/$NAME.vcf

sed -i -e "s/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$NAME/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAC\tAF\tAN\tBaseQRankSum\tClippingRankSum\tDP\tExcessHet\tFS\tInbreedingCoeff\tMLEAC\tMLEAF\tMQ\tMQRankSum\tQD\tReadPosRankSum\tSOR\tAllele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tVARIANT_CLASS\tSYMBOL_SOURCE\tHGNC_ID\tENSP\tREFSEQ_MATCH\tSOURCE\tGIVEN_REF\tUSED_REF\tBAM_EDIT\tSIFT\tPolyPhen\tHGVS_OFFSET\tAF\tCLIN_SIG\tSOMATIC\tPHENO\tGT\tAD\tDP\tGQ\tPL/g" -e "s/;\tGT/\tGT/g" -e "s/GT:AD:DP:GQ:PL\t//g" $OUT_PATH/VCF/$NAME.vcf
sed -i -e "s/\(^chr[[:alnum:]]*\t[[:digit:]]*\t\)\.\t/\1$NAME\t/g" -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/\(|.*\):\(.*|\)/\1[]\2/g' -e 's/:/\t/g' -e 's/\[\]/:/g' -e 's/;/\t/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/||/|-|/g' -e 's/|/\t/g' -e "s/AC=\|AF=\|AN=\|BaseQRankSum=\|ClippingRankSum=\|DP=\|ExcessHet=\|FS=\|InbreedingCoeff=\|MLEAC=\|MLEAF=\|MQ=\|MQRankSum=\|QD=\|ReadPosRankSum=\|SOR=\|ANN=//g" -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' $OUT_PATH/VCF/$NAME.vcf
cut -f1,2,3,4,5,6,7,8,10,25,27,32,33,34,35,58,61,62,63,64,65,66 $OUT_PATH/VCF/$NAME.vcf > $OUT_PATH/TSV/$NAME.tsv
done


while read -r line
do

GENE=$line
echo $GENE

#echo -e "/usr/bin/python ~/git/CMG/SCRIPT_CMG/SCRIPT_PYTHON/DB_per_Gene.py --path $PATH --gene $GENE --out $OUT/$GENE --paz_list $PAZ_LIST --tipo $TIPO"
/usr/bin/python /home/jarvis/git/CMG/SCRIPT_CMG/SCRIPT_PYTHON/DATABASE_SCRIPTING/DB_per_Gene.py --path $OUT_PATH/TSV --gene $GENE --out $OUT_PATH/GENI/$GENE --paz_list $PAZ_LIST --tipo $PANNELLO

done < "$GENELIST"

