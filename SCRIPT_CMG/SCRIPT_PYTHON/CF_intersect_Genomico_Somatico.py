import argparse 
def main():

	parser = argparse.ArgumentParser('Esegue l intersezione tra le varianti dell analisi germinale e quelle dell analisi somatica utilizzando solo chrom e pos.  Output is to stdout.')
	parser.add_argument('-g', '--genomic', help="tsv germinale")
	parser.add_argument('-s', '--somatic', help="tsv somatico")
	
	global opts 
	opts = parser.parse_args()
	
	varianti = dict()
	genomic=open(opts.genomic,'r')
	for posiz in genomic:
		posiz=posiz.rstrip()
		chr=posiz.split('\t')[0]
		pos=posiz.split('\t')[1]
		varianti[chr+pos]=chr+'\t'+pos

	somatic=open(opts.somatic,'r')
	for posiz in somatic:
		posiz=posiz.rstrip()
		chr=posiz.split('\t')[0]
		pos=posiz.split('\t')[1]
		if varianti.has_key(chr+pos):
			varianti[chr+pos]=posiz

	for var in varianti.keys():
		print varianti.get(var)
			
	
main()
