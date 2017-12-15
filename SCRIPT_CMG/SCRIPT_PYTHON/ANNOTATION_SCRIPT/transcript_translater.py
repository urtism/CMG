
import argparse

def search_genes(file,geni):
	lista_geni_tr = open(opts.out,'w')
	
	lista_geni = []
	for gene in geni:
		gene = gene.rstrip()
		lista_geni = lista_geni + [gene]
	
	for line in file:
		line = line.rstrip()
		if line.startswith('Ensembl'):
			continue
		gene_id = line.split(',')[0]
		transcript_id = line.split(',')[1]
		refseq_id = line.split(',')[2]
		gene_name = line.split(',')[3]
		if gene_name in lista_geni:
			lista_geni_tr.write('\t'.join([gene_name,gene_id,transcript_id,refseq_id]) +'\n')

		




def main():
	parser = argparse.ArgumentParser('Adds number of exon for each interval')
	
	parser.add_argument('-f','--file',help="variants file")
	parser.add_argument('-g','--geni',help="targets file",default= None)
	parser.add_argument('-o','--out',help="new variants file path",default= None)

	global opts
	opts = parser.parse_args()
	search_genes(open(opts.file,'r'),open(opts.geni,'r'))

main()