import os
import glob2
import argparse
import shutil


def main():
	parser = argparse.ArgumentParser('incrocia due file di varianti e vede in quali file sono presenti')
	
	parser.add_argument('-i','--file1',help="path del file di varianti da illumina")
	parser.add_argument('-s','--file2',help="path del file di varianti da sanger")
	parser.add_argument('-o','--out',help="file di output")

	global opts
	opts = parser.parse_args()
	lista_snp=dict()

	file1=open(opts.file1,'r')
	file2=open(opts.file2,'r')
	out=open(opts.out,'w')
	
	for snp in file1:
		snp=snp.rstrip()
		if snp.startswith('CHROM'):
			header=snp.split('\t')
		else:
			chrom=snp.split('\t')[header.index('CHROM')]
			pos=snp.split('\t')[header.index('POS')]
			ref=snp.split('\t')[header.index('REF')]
			alt=snp.split('\t')[header.index('ALT')]
			id_var='\t'.join([chrom,pos,ref,alt])
			lista_snp[id_var]=['file1','.']


	for snp in file2:
		snp=snp.rstrip()
		if snp.startswith('CHROM'):
			header=snp.split('\t')
		else:
			chrom=snp.split('\t')[header.index('CHROM')]
			pos=snp.split('\t')[header.index('POS')]
			ref=snp.split('\t')[header.index('REF')]
			alt=snp.split('\t')[header.index('ALT')]
			id_var='\t'.join([chrom,pos,ref,alt])
			if id_var in lista_snp.keys():
				lista_snp[id_var][1]='file2'
			else:
				lista_snp[id_var]=['.','file2']

	for id_var in lista_snp.keys():
		out.write(id_var +'\t'+'\t'.join(lista_snp[id_var]) +'\n')
main()
