#!bin/bash

cd /home/minime/DB_ANNOTAZIONE/phastCons46way/phastCons46way_placental

for file in *.chr.tsv
do
  mv "$file" "${file//phastCons46way/phastCons}"
done



cd /home/minime/DB_ANNOTAZIONE/phastCons46way/phastCons46way_primate

for file in *.chr.tsv
do
  mv "$file" "${file//phastCons46way/phastCons}"
done




cd /home/minime/DB_ANNOTAZIONE/phastCons100way/phastCons100way_vertebrate

for file in *.chr.tsv
do
  mv "$file" "${file//phastCons100way/phastCons}"
done



cd /home/minime/DB_ANNOTAZIONE/phyloP46way/phyloP46way_placental

for file in *.chr.tsv
do
  mv "$file" "${file//phyloP46way/phyloP}"
done


cd /home/minime/DB_ANNOTAZIONE/phyloP46way/phyloP46way_primate

for file in *.chr.tsv
do
  mv "$file" "${file//phyloP46way/phyloP}"
done



cd /home/minime/DB_ANNOTAZIONE/phyloP100way/phyloP100way_vertebrate

for file in *.chr.tsv
do
  mv "$file" "${file//phyloP100way/phyloP}"
done