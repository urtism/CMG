import argparse

if __name__ == '__main__':

	parser = argparse.ArgumentParser('Genera un file per ogni paziente che ha piU di due varianti nel gene in entrata')
	parser.add_argument('-g','--gene',help="file tsv di varianti nel gene db")
	parser.add_argument('-n','--num_var',help="soglia di varianti minimo")
	parser.add_argument('-o','--out',help="path di out")

	global opts
	opts = parser.parse_args()
	gene=open(opts.gene,'r')
	pazienti=dict()
	excluded=['synonymous_variant','3_prime_UTR_variant','intron_variant','splice_region_variant&intron_variant','downstream_gene_variant']
	for line in gene:
		if line.startswith('CHROM'):
			header=line.rstrip().split('\t')
		elif not line.split('\t')[header.index('CONSEQUENCE')] in excluded:
			arr_paz= line.rstrip().split('\t')[-2].split(';')+line.rstrip().split('\t')[-1].split(';')
			for paz in arr_paz:
				if paz in pazienti.keys():
					pazienti[paz]+=[line.rstrip().split('\t')]
				else:
					pazienti[paz]=[line.rstrip().split('\t')]

	header[2:2]=['GT']
	for paz in pazienti.keys():
		if len(pazienti[paz])>=int(opts.num_var):
			nome_paz='_'.join(paz.split(' '))
			out=open(opts.out+'/'+nome_paz+'.tsv','w')
			out.write('\t'.join(header) +'\n')
			for line in pazienti[paz]:
				if paz in line[-1]:
					line[2:2]=['OMO']
				elif paz in line[-2]:
					line[2:2]=['HET']
				out.write('\t'.join(line) +'\n')
			out.close()





