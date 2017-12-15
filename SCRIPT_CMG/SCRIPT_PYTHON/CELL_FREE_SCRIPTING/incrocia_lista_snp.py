import os
import glob2
import argparse
import shutil


def main():
	parser = argparse.ArgumentParser('Adds number of exon for each interval')
	
	parser.add_argument('-i','--snp1',help="path degli snp del panello 1")
	parser.add_argument('-t','--snp2',help="path degli snp del panello 2")
	parser.add_argument('-o','--out',help="path di output")

	global opts
	opts = parser.parse_args()
	lista_snp=dict()

	pan1=open(opts.snp1,'r')
	pan2=open(opts.snp2,'r')
	out=open(opts.out,'w')
	
	for snp in pan1:
		snp=snp.rstrip()
		if snp.startswith('#'):
			continue
		else:
			SNP=snp.split('\t')[2]
			lista_snp[SNP]=['SI','NO']


	# for snp in pan1:
	# 	snp=snp.rstrip()
	# 	if snp.startswith('#'):
	# 		continue
	# 	else:
	# 		SNP=snp.split('\t')[0]+'\t'+snp.split('\t')[1]
	# 		lista_snp[SNP]=['SI','NO']

	# for snp in pan2:
	# 	snp=snp.rstrip()
	# 	if snp.startswith('track'):
	# 		continue
	# 	else:
	# 		SNP=snp.split('\t')[0]+'\t'+snp.split('\t')[2]
	# 		#print snp,SNP in lista_snp.keys()
	# 		if SNP in lista_snp.keys():
	# 			lista_snp[SNP]=['SI','SI']
	# 		else:
	# 			lista_snp[SNP]=['NO','SI']

	for snp in pan2:
		snp=snp.rstrip()
		if snp.startswith('Locus'):
			continue
		else:
			SNP=snp.split('\t')[0]
			#print snp,SNP in lista_snp.keys()
			if SNP in lista_snp.keys():
				lista_snp[SNP]=['SI','SI']
			else:
				lista_snp[SNP]=['NO','SI']


	#out.write('\t'.join(['CHROM','POS','SNP_TRAP','HID-Ion'])+'\n')
	out.write('\t'.join(['LOCUS','SNP_TRAP','IISNPs'])+'\n')
	
	#print lista_snp.keys()
	for snp in lista_snp.keys():
		out.write('\t'.join([snp]+lista_snp.get(snp))+'\n')


main()
