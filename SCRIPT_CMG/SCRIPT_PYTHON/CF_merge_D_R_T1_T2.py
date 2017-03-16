import argparse


def read_snp(varianti):
	lista_snp = open(opts.snp_list,'r')
	for line in lista_snp:
		line=line.rstrip()
		if line.startswith('#'):
			continue
		else:
			chrom = line.split('\t')[0]
			pos = line.split('\t')[1]
			var=['','','','','']
			varianti[chrom + '\t'+ pos]=var
	lista_snp.close()

def read_donor(varianti,header):
	lista_snp = open(opts.donor,'r')
	for line in lista_snp:
		line=line.rstrip()
		if line.startswith('CHROM'):
			HEADER=line.split('\t')
			header = header + '\t' +'\t'.join(['GT_GATK_D','GT_VARSCAN_D','GT_FREEBAYES_D','DP_D','AF_D','QB_D'])
		else:
			chrom = line.split('\t')[0]
			pos = line.split('\t')[1]
			ref = line.split('\t')[3]
			alt = line.split('\t')[4]
			gt = '\t'.join(line.split('\t')[HEADER.index('GT_GATK'):HEADER.index('GT_Freebayes') +1])
			#print gt
			#dp = '\t'.join(line.split('\t')[11:14])
			dp = line.split('\t')[14]
			qb = line.split('\t')[18]
			AF = line.split('\t')[-1]
			if varianti.has_key(chrom + '\t'+ pos):

				donor = '\t'.join([gt,dp,AF,qb])
				var = varianti.get(chrom + '\t'+ pos)
				if var[0]=='':
					var[0]='.\t'+ref+'\t'+alt
				var[1] = donor
				varianti[chrom + '\t'+ pos]=var

	lista_snp.close()
	return header

def read_rec(varianti,header):
	lista_snp = open(opts.receiver,'r')
	for line in lista_snp:
		line=line.rstrip()
		if line.startswith('CHROM'):
			HEADER=line.split('\t')
			header = header + '\t' + '\t'.join(['GT_GATK_R','GT_VARSCAN_R','GT_FREEBAYES_R','DP_R','AF_R','QB_R'])
		else:
			chrom = line.split('\t')[0]
			pos = line.split('\t')[1]
			ref = line.split('\t')[3]
			alt = line.split('\t')[4]
			gt = '\t'.join(line.split('\t')[HEADER.index('GT_GATK'):HEADER.index('GT_Freebayes') +1])
			#print gt
			#dp = '\t'.join(line.split('\t')[11:14])
			dp = line.split('\t')[14]
			qb = line.split('\t')[18]
			AF = line.split('\t')[-1]
			if varianti.has_key(chrom + '\t'+ pos):
				rec = '\t'.join([gt,dp,AF,qb])
				var = varianti.get(chrom + '\t'+ pos)
				if var[0]=='':
					var[0]='.\t'+ref+'\t'+alt
				var[2] = rec
				varianti[chrom + '\t'+ pos]=var

	lista_snp.close()
	return header

def read_sT1(varianti,header):
	lista_snp = open(opts.sT1,'r')
	for line in lista_snp:
		line=line.rstrip()
		if line.startswith('CHROM'):
			HEADER=line.split('\t')
			header = header + '\t' + '\t'.join(['SOMATIC_VARSCAN_T1','FILTER_MUTECT_T1','STATUS_VARDICT_T1','GT_TUM_MUTECT_T1','GT_TUM_VARSCAN_T1','GT_TUM_VARDICT_T1',
										'GT_NORM_MUTECT_T1','GT_NORM_VARSCAN_T1','GT_NORM_VARDICT_T1',
										'DP_TUM_T1','AF_TUM_T1','QB_TUM_T1','DP_NORM_T1','AF_NORM_T1','QB_NORM_T1','DELTA_MEDIA_T1'])
		else:
			
			chrom = line.split('\t')[0]
			pos = line.split('\t')[1]
			ref = line.split('\t')[3]
			alt = line.split('\t')[4]
			somatic_varscan =line.split('\t')[8]
			filter_mutect = line.split('\t')[10]
			status_vardict = line.split('\t')[11]

			gt_t = '\t'.join(line.split('\t')[HEADER.index('GT_Mutect'):HEADER.index('GT_Vardict') +1])
			gt_n = '\t'.join(line.split('\t')[15:18])
			dp_t = line.split('\t')[18]
			af_t = line.split('\t')[22]
			qb_t = line.split('\t')[25]
			dp_n = line.split('\t')[26]
			af_n = line.split('\t')[30]
			qb_n = line.split('\t')[33]
			delta = line.split('\t')[-1]
			if varianti.has_key(chrom + '\t'+ pos):
				t1 = '\t'.join([somatic_varscan,filter_mutect,status_vardict,gt_t,gt_n,dp_t,af_t,qb_t,dp_n,af_n,qb_n,delta])
				var = varianti.get(chrom + '\t'+ pos)
				if var[0]=='':
					var[0]='.\t'+ref+'\t'+alt
				var[3] = t1
				varianti[chrom + '\t'+ pos]=var
			

	lista_snp.close()
	return header

def read_sT2(varianti,header):
	lista_snp = open(opts.sT2,'r')
	for line in lista_snp:
		line=line.rstrip()

		if line.startswith('CHROM'):
			HEADER=line.split('\t')
			header = header +  '\t' +'\t'.join(['SOMATIC_VARSCAN_T2','FILTER_MUTECT_T2','STATUS_VARDICT_T2','GT_TUM_MUTECT_T2','GT_TUM_VARSCAN_T2','GT_TUM_VARDICT_T2',
										'GT_NORM_MUTECT_T2','GT_NORM_VARSCAN_T2','GT_NORM_VARDICT_T2',
										'DP_TUM_T2','AF_TUM_T2','QB_TUM_T2','DP_NORM_T2','AF_NORM_T2','QB_NORM_T2','DELTA_MEDIA_T2'])
		else:
			chrom = line.split('\t')[0]
			pos = line.split('\t')[1]
			ref = line.split('\t')[3]
			alt = line.split('\t')[4]
			somatic_varscan =line.split('\t')[8]
			filter_mutect = line.split('\t')[10]
			status_vardict = line.split('\t')[11]

			gt_t = '\t'.join(line.split('\t')[12:15])
			gt_n = '\t'.join(line.split('\t')[15:18])
			dp_t = line.split('\t')[18]
			af_t = line.split('\t')[22]
			qb_t = line.split('\t')[25]
			dp_n = line.split('\t')[26]
			af_n = line.split('\t')[30]
			qb_n = line.split('\t')[33]
			delta = line.split('\t')[-1]
			if varianti.has_key(chrom + '\t'+ pos):
				t2 = '\t'.join([somatic_varscan,filter_mutect,status_vardict,gt_t,gt_n,dp_t,af_t,qb_t,dp_n,af_n,qb_n,delta])
				var = varianti.get(chrom + '\t'+ pos)
				if var[0]=='':
					var[0]='.\t'+ref+'\t'+alt
				var[4] = t2
				varianti[chrom + '\t'+ pos]=var

	lista_snp.close()
	return header

def check(varianti):
	for variante in varianti.keys():
		var = varianti.get(variante)
		if var[0]=='':
			var[0] = '\t'.join(['.','.','.'])
		if var[1]=='':
			var[1] = '\t'.join(['.','.','.','.','.','.'])
		if var[2]=='':
			var[2] = '\t'.join(['.','.','.','.','.','.'])
		if var[3]=='':
			var[3] = '\t'.join(['.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.'])
		if var[4]=='':
			var[4] = '\t'.join(['.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.'])


def print_varianti(varianti,out,header):
	out.write(header+'\n')
	for variante in varianti.keys():
		var = varianti.get(variante)
		out.write('\t'.join([variante] + var) + '\n')





def main():

	parser = argparse.ArgumentParser('Parse VCF output from Variant callers to output a variant_dataset.txt.  Output is to stdout.')
	parser.add_argument('-d', '--donor', help="Freebayes vcf output file name")
	parser.add_argument('-r', '--receiver', help="gatk vcf output file name")
	parser.add_argument('--sT1', help="Varscan vcf output file name")
	parser.add_argument('--sT2', help="Lista di features da stampare")
	parser.add_argument('-o', '--out',help="path di output")
	parser.add_argument('-l', '--snp_list',help="gvcf path")

	global opts 
	opts = parser.parse_args()
	
	varianti = dict()
	header='CHROM\tPOS\tID\tREF\tALT'
	read_snp(varianti)
	header=read_donor(varianti,header)
	header=read_rec(varianti,header)
	header=read_sT1(varianti,header)
	header=read_sT2(varianti,header)
	check(varianti)
	out=open(opts.out,'w')
	print_varianti(varianti,out,header)
	
	
main()