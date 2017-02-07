vcf = open('/home/minime/Scrivania/SCRIPT_PYTHON/INPUT/20161215_Cardio_GATK_Filter_Transcripts.vcf','r')

out = open('/home/minime/Scrivania/SCRIPT_PYTHON/INPUT/20161215_Cardio_GATK_Filter_Transcripts_No_Duplicati.vcf','w')

chrom = vcf.readlines()
chrom2 = chrom
stampa = []

for line in chrom:
	line.rstrip()
	line = line.split('\t')
	print line

	for ann in line:
		print ann

		for riga in chrom2:
			riga.rstrip()
			riga = riga.split('\t')

			#for elem in riga:

				#if elem[0] != ann[0] and elem[1] != ann[1] and elem[3] != ann[3] and elem[4] != ann[4]:
					#stampa += [elem]
					#print stampa [0]
					#out.write('\t'.join(stampa) + '\n')



