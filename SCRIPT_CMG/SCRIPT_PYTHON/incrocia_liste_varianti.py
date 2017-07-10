import os
import glob2
import argparse
import shutil


def main():
	parser = argparse.ArgumentParser('Adds number of exon for each interval')
	
	parser.add_argument('-i','--file1',help="path del file di varianti da illumina")
	parser.add_argument('-s','--file2',help="path del file di varianti da sanger")
	parser.add_argument('-o','--out',help="path di output")

	global opts
	opts = parser.parse_args()
	lista_snp=dict()

	file1=open(opts.file1,'r')
	file2=open(opts.file2,'r')
	inters=open(opts.out+'.intersect.tsv','w')
	inters_all=open(opts.out+'.intersect.ALL.tsv','a')
	sanger=open(opts.out+'.sanger.tsv','w')

	#inters.write('\t'.join(['CHROM','POS','ID','REF','ALT','GENE','HGVSc','HGVSp'])+'\n')
	sanger.write('\t'.join(['CHROM','POS','REF','ALT'])+'\n')
	
	for snp in file1:
		snp=snp.rstrip()
		if snp.startswith('CHROM'):
			header=snp.split('\t')
			inters.write(snp +'\n')
		else:
			chrom=snp.split('\t')[header.index('CHROM')]
			pos=snp.split('\t')[header.index('POS')]
			ref=snp.split('\t')[header.index('REF')]
			id=snp.split('\t')[header.index('ID')]
			alt=snp.split('\t')[header.index('ALT')]
			HGVSc=snp.split('\t')[header.index('HGVSc')]
			HGVSp=snp.split('\t')[header.index('HGVSp')]
			gene=snp.split('\t')[header.index('SYMBOL')]
			id_var='\t'.join([chrom,pos,ref,alt])
			lista_snp[id_var]=snp

    if not any(x.startswith('CHROM') for x in inters_all):
        f.write('\t'.join(header) + '\n')

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
				inters.write(lista_snp[id_var]+'\n')
				inters_all.write(lista_snp[id_var]+'\n')
			else:
				sanger.write('\t'.join([chrom,pos,ref,alt])+'\n')


main()
