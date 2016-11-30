import argparse

def extract_info(vcf):
	out = open (opts.out,'w')

	for line in vcf:
		line = line.rstrip()
		if line.startswith('##'):
			continue
		elif line.startswith('#CHROM'):
			header = ['CHROM','POS','REF','ALT','QUAL','FILTER','INFO','GT_TUM','GT_NORM','DP_TUM','DP_NORM','AD_TUM','AD_NORM','FREQ_TUM','FREQ_NORM','HIAF_TUM','HIAF_NORM','BQUAL_TUM','BQUAL_NORM','MQ_TUM','MQ_NORM']
			out.write('\t'.join(header) + '\n')
			continue
		else:
			line_split = line.split('\t')
			chrom = line_split[0]
			pos = line_split[1]
			ref = line_split[3]
			alt = line_split[4]
			qual = line_split[5]
			filter = line_split[6]
			info  = line_split[7]
			format = line_split[8]
			tumor = line_split[9]
			normal = line_split[10]
			tumor_split = tumor.split(':')
			normal_split = normal.split(':')
			#GT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM
			variant = [chrom,pos,ref,alt,qual,filter,info,tumor_split[0],normal_split[0],tumor_split[1],normal_split[1],tumor_split[5],normal_split[5],tumor_split[6],normal_split[6],tumor_split[16],normal_split[16],tumor_split[10],normal_split[10],tumor_split[14],normal_split[14]]
			out.write('\t'.join(variant) + '\n')

def main():
	parser = argparse.ArgumentParser('find info in vcf')
	parser.add_argument('-f','--infile',help="vcf file")
	parser.add_argument('-o','--out',help="variants file")

	global opts
	opts = parser.parse_args()

	extract_info(open(opts.infile,'r'))

main()