import argparse




if __name__ == '__main__':

	parser = argparse.ArgumentParser('Controlal le varianti in comune tra due tsv')
	parser.add_argument('-i1','--input1',help="File tab delimited 1 contenente le varianti da controllare")
	parser.add_argument('-i2','--input2',help="File tab delimited 2 contenente le varianti da controllare")
	parser.add_argument('-gl','--gene_list',help="lista di geni da controllare",default=None)
	parser.add_argument('-o','--out',help="lista di varianti in comune tra input1 e input2 presenti nei geni della lista")

	global opts
	opts = parser.parse_args()

	tsv1 = open(opts.input1,'r')
	tsv2 = open(opts.input2,'r')

	out= open(opts.out,'w')
	out.write('\t'.join(["CHROM","POS","REF","ALT","FILTER","AC","Consequence","SYMBOL","HGVSc","HGVSp","Existing_variation","SIFT","PolyPhen","GMAF","CLIN_SIG","GT","AD"]) +'\n')
	gene_l=[]

	tsv1_l = dict()
	tsv2_l = dict()

	if opts.gene_list != None:
	
		gene_list = open(opts.gene_list,'r')

		for line in gene_list:
			gene_l += [line.rstrip()]


	for line in tsv1:
		line=line.rstrip()
		if line.startswith('CHROM'):
			header = line.split('\t')
		else:
			CHROM,POS,ID,REF,ALT,QUAL,FILTER,AC,AF,AN,BaseQRankSum,ClippingRankSum,DP,ExcessHet,FS,InbreedingCoeff,MLEAC,MLEAF,MQ,MQRankSum,QD,ReadPosRankSum,SOR,Allele,Consequence,SYMBOL,Feature,HGVSc,HGVSp,Amino_acids,Codons,Existing_variation,VARIANT_CLASS,SIFT,PolyPhen,GMAF,CLIN_SIG,PUBMED,GT,AD,DP,GQ,PL=line.split('\t')
			id = CHROM,POS,REF,ALT
			if gene_l != []:
				if SYMBOL in gene_l:
					tsv1_l[id] = '\t'.join([CHROM,POS,REF,ALT,FILTER,AC,Consequence,SYMBOL,HGVSc,HGVSp,Existing_variation,SIFT,PolyPhen,GMAF,CLIN_SIG,GT,AD])
			else:
				tsv1_l[id] = '\t'.join([CHROM,POS,REF,ALT,FILTER,AC,Consequence,SYMBOL,HGVSc,HGVSp,Existing_variation,SIFT,PolyPhen,GMAF,CLIN_SIG,GT,AD])

	

	for line in tsv2:
		line=line.rstrip()
		if line.startswith('CHROM'):
			header = line.split('\t')
		else:
			CHROM,POS,ID,REF,ALT,QUAL,FILTER,AC,AF,AN,BaseQRankSum,ClippingRankSum,DP,ExcessHet,FS,InbreedingCoeff,MLEAC,MLEAF,MQ,MQRankSum,QD,ReadPosRankSum,SOR,Allele,Consequence,SYMBOL,Feature,HGVSc,HGVSp,Amino_acids,Codons,Existing_variation,VARIANT_CLASS,SIFT,PolyPhen,GMAF,CLIN_SIG,PUBMED,GT,AD,DP,GQ,PL=line.split('\t')
			id = CHROM,POS,REF,ALT
			if id in tsv1_l.keys():
				out.write('\t'.join([CHROM,POS,REF,ALT,FILTER,AC,Consequence,SYMBOL,HGVSc,HGVSp,Existing_variation,SIFT,PolyPhen,GMAF,CLIN_SIG,GT,AD]) +'\n')
