#!bin/bash

cd /home/minime/DB_ANNOTAZIONE/phastCons46way/phastCons46way_placental_Cancer

for (( a=1; a<23; a++ ))
do
grep "chr$a	" phastCons46way_placental_Cancer.tsv > phastCons46way_placental_Cancer_chr$a.tsv
done

grep "chrX	" phastCons46way_placental_Cancer.tsv > phastCons46way_placental_Cancer_chrX.tsv
grep "chrY	" phastCons46way_placental_Cancer.tsv > phastCons46way_placental_Cancer_chrY.tsv



cd /home/minime/DB_ANNOTAZIONE/phastCons46way/phastCons46way_placental_Cardio

for (( a=1; a<23; a++ ))
do
grep "chr$a	" phastCons46way_placental_Cardio.tsv > phastCons46way_placental_Cardio_chr$a.tsv
done

grep "chrX	" phastCons46way_placental_Cardio.tsv > phastCons46way_placental_Cardio_chrX.tsv
grep "chrY	" phastCons46way_placental_Cardio.tsv > phastCons46way_placental_Cardio_chrY.tsv



cd /home/minime/DB_ANNOTAZIONE/phastCons46way/phastCons46way_primate_Cancer

for (( a=1; a<23; a++ ))
do
grep "chr$a	" phastCons46way_primate_Cancer.tsv > phastCons46way_primate_Cancer_chr$a.tsv
done

grep "chrX	" phastCons46way_primate_Cancer.tsv > phastCons46way_primate_Cancer_chrX.tsv
grep "chrY	" phastCons46way_primate_Cancer.tsv > phastCons46way_primate_Cancer_chrY.tsv



cd /home/minime/DB_ANNOTAZIONE/phastCons46way/phastCons46way_primate_Cardio

for (( a=1; a<23; a++ ))
do
grep "chr$a	" phastCons46way_primate_Cardio.tsv > phastCons46way_primate_Cardio_chr$a.tsv
done

grep "chrX	" phastCons46way_primate_Cardio.tsv > phastCons46way_primate_Cardio_chrX.tsv
grep "chrY	" phastCons46way_primate_Cardio.tsv > phastCons46way_primate_Cardio_chrY.tsv


