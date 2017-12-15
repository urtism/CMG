import argparse


def main():
	parser = argparse.ArgumentParser('Adds number of exon for each interval')
	
	parser.add_argument('-f','--file',help="path della lista dei geni")
	parser.add_argument('-o','--omim',help="path di omim")

	global opts
	opts = parser.parse_args()
	lista_geni=open(opts.file,'r')
	print '\t'.join(['SYMBOL','PANEL','GENE_NAME','PHENOTYPE'])
	
	for gene in lista_geni:
		if gene.startswith('GENE'):
			continue
		else:
			gene=gene.rstrip()
			find=0

			omim=open(opts.omim,'r')
			for line in omim:
				line=line.rstrip()
				if line.startswith('GeneSymbols'):
					continue
				else:
					
					line_split=line.split('\t')
					symbol=line_split[0]
					try:
						nome=line_split[1]
					except:
						nome='.'
					try:
						fenotipo=line_split[2]
					except:
						fenotipo='.'
					
					if gene.split('\t')[0] in symbol.split(','):
						print gene + '\t' + '\t'.join([nome,fenotipo])
						find=1
						break
			omim.close()
			if find==0:
				print gene + '\t' + '\t'.join(['.','.'])

main()