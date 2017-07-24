import argparse

def clinsig_map(CLNSIG):
	map={'0' : 'Uncertain significance', '1' :'not provided', '2':'Benign', '3' : 'Likely benign', '4' : 'Likely pathogenic', '5' : 'Pathogenic', '6' : 'drug response', '7' : 'histocompatibility', '255' : 'other'}
	gne=[]
	for el in CLNSIG.split('|'):
		i=[]
		for il in el.split(','):
			i+=[map[il]]
		gne+=[','.join(i)]
	return '|'.join(gne)



def main():
	parser = argparse.ArgumentParser('estrae informazioni dal vcf di clinvar filtrando per gene. Out in stdout')
	
	parser.add_argument('-f','--vcf',help="path del vcf di clinvar")
	parser.add_argument('-g','--gene',help="gene ")
	global opts
	opts = parser.parse_args()
	vcf=open(opts.vcf,'r')

	print '\t'.join(['CHROM','POS','ID','REF','ALT','GENE','MAF','CLNSIG','CLNDBN'])

	for line in vcf:
		line=line.rstrip()
		if line.startswith('#'):
			continue
		else:
			line_split=line.split('\t')
			chr=line_split[0]
			pos=line_split[1]
			id=line_split[2]
			ref=line_split[3]
			alt=line_split[4]
			qual=line_split[5]
			filter=line_split[6]
			info=line_split[7]
			CLNSIG='.'
			CLNDBN='.'
			MAF='.'
			GENE='.'
			for elem in info.split(';'):
				if elem.startswith('CLNSIG='):
					CLNSIG=clinsig_map(elem.split('=')[1])

				elif elem.startswith('CLNDBN='):
					CLNDBN=elem.split('=')[1]
				elif elem.startswith('CAF='):
					MAF=elem.split('=')[1]
				elif elem.startswith('GENEINFO='):
					GENE=(elem.split('=')[1]).split(':')[0]

			if GENE == (opts.gene).rstrip():
				print '\t'.join([chr,pos,id,ref,alt,GENE,MAF,CLNSIG,CLNDBN])


main()