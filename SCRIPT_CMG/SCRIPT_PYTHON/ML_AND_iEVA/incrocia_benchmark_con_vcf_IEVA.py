import argparse




def main():
	
	parser = argparse.ArgumentParser('incrocia il vcf con la var_list mantenendo soltanto le varianti a + o meno una window e ristampando in formato vcf.')
	parser.add_argument('-l', '--var_list', help="GATK vcf output file name")
	parser.add_argument('-v', '--vcf', help="Freebayes vcf output file name")
	parser.add_argument('-w', '--window', help="finestra entro il quale una variante viene mantenuta [[+-w]")
	parser.add_argument('-o', '--out', help="file name in output. It returns file_name.vcf ")
	
	global opts 
	opts = parser.parse_args()

	vars=dict()
	samples_list=dict()
	with open(opts.var_list,'r') as vlist:
		for line in vlist:
			if line.startswith('Paziente'):
				continue
			#print line.rstrip().split('\t')
			else:
				Paziente,ID,CHR,POS,REF,ALT,HGVSc,HGVSp,Gene,Class,Numero_Esone = line.rstrip().split('\t')
				id_var='\t'.join([CHR,POS])
				vars[id_var]=[CHR,POS,ID]
				samples_list[ID]=''

	outvcf = open(opts.out,'w')
	non_printare = 1
	with open(opts.vcf,'r') as invcf:
		for line in invcf:
			line = line.rstrip()
			if line.startswith('##'):
				outvcf.write(line+'\n')
			elif line.startswith('#CHROM'):
				outvcf.write(line+'\n')
				samples = line.split('\t')[9:]
				for s in samples:
					data,num,pannello = s.split('_')
					for sample in samples_list.keys():
						datasample,numsample,pannellosample = sample.split('_')
						if datasample == data and numsample == num:
							non_printare = 0
							if pannellosample != pannello:
								print sample,'errato, non presente nel vcf'
			else:
				if non_printare == 0:
					line_split = line.split('\t')
					CHR,POS = line_split[:2]
					for pos in range(int(POS)-int(opts.window),int(POS)+int(opts.window),1):
						if '\t'.join([CHR,str(pos)]) in vars.keys():
							outvcf.write(line+'\n')

main()

