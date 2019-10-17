import argparse
if __name__ == '__main__':

	parser = argparse.ArgumentParser('estrae da un vcf annotato le varianti.\n')
	parser.add_argument('-i', '--vcf', default=None, help="vcf input")
	parser.add_argument('-o', '--out', default=None, help="output file")

	global opts
	opts = parser.parse_args()

	vcf = open(opts.vcf,'r')
	out = open(opts.out,'w')
	out.write('\t'.join(['CHROM','POS','ID','REF','ALT','Consequence','GENE','HGVSc','HGVSp','clinvar'])+'\n')

	for line in vcf:
		line=line.rstrip()
		if line.startswith('##INFO=<ID=ANN,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: '):
			annotation_header = line.split('Format: ')[1].split('|')
		elif line.startswith('#CHROM'):
			header = line.split('\t')
		elif line.startswith('chr'):
			CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE=line.split('\t')
			ANN=INFO.split('ANN=')[1].split('|')
			cons = ANN[annotation_header.index('Consequence')]
			gene = ANN[annotation_header.index('SYMBOL')]
			hgvsc = ANN[annotation_header.index('HGVSc')]
			hgvsp = ANN[annotation_header.index('HGVSp')]
			clinvar = ANN[annotation_header.index('CLIN_SIG')]

			out.write('\t'.join([CHROM,POS,ID,REF,ALT,cons,gene,hgvsc,hgvsp,clinvar])+'\n')