#!bin/bash

cd /home/minime/Scrivania/TEST_ANNOTAZIONE_DB

DATE=$(date +"%Y%m%d");
STARTTIME=$(date +%s)


python ~/git/CMG/SCRIPT_CMG/SCRIPT_PIPELINE/vcf_Annotation.py \
-i ~/Scrivania/SCRIPT_PYTHON/INPUT/20170118_Cancer_TOTAL.vcf \
-clnv ~/DB_ANNOTAZIONE/ClinVar/clinvar.target_cancer.vcf,CLNORIGIN,CLNDBN,CLNSIG \
-ESP ~/DB_ANNOTAZIONE/ExomeVariantServer/ESP6500SI-V2-SSA137.GRCh38-liftover.allchr.target_cancer+-1000.snps_indels.vcf,ESP_EA,ESP_AA,ESP_All \
-hvar ~/DB_ANNOTAZIONE/Humsavar/humsavar.cancer.txt,Variant_type \
-gerp ~/DB_ANNOTAZIONE/GERP++/,RS_Score \
-P Cancer \
--phastCons ~/DB_ANNOTAZIONE/phastCons/,placental,vertebrate,primate \
--phyloP ~/DB_ANNOTAZIONE/phyloP/,placental,vertebrate,primate \
-o ~/Scrivania/SCRIPT_PYTHON/OUTPUT/20170118_TOTAL_ANN_DB.vcf


ENDTIME=$(date +%s)

ENDTIME=$(date +%s)
echo "TEMPO TOTALE DI PROCESSING: $(($ENDTIME - $STARTTIME)) sec"
echo "TEMPO TOTALE DI PROCESSING: $(($ENDTIME - $STARTTIME)) sec" > /home/minime/Scrivania/TEST_ANNOTAZIONE_DB/Time_log.txt
printf "\n\n"


#5336 secondi x 12 Cancer con tutte le annotazioni