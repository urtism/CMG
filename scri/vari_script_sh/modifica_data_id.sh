#!bin/bash


cd /home/minime/Scrivania/STORICO_CONFERME_RUN/CARDIO/20151113_Run_14_Cardio_Conferme_2.1
#cd /home/minime/Scrivania/STORICO_CONFERME_RUN/CARDIO/20151119_Run_16_Cardio_Conferme_2.1
#cd /home/minime/Scrivania/a/20151113_Run_14_Cardio_Conferme_2.1

for file in *.tsv
do
echo "$file"
	sed -i -e "s/20151119/20151113/g" "/home/minime/Scrivania/STORICO_CONFERME_RUN/CARDIO/20151113_Run_14_Cardio_Conferme_2.1/$file"
done




