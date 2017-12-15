import argparse
import os
import glob2
import sys


def main():
	parser = argparse.ArgumentParser('filtra il file --file per la lista di geni --gene_list')
	
	parser.add_argument('-f','--file',help="database dei trascritti")
	parser.add_argument('-g','--gene_list',help="lista di gene su cui filtrare il database dei strascritti")

	global opts
	opts = parser.parse_args()
	geni=[]
	gene_list=open(opts.gene_list,'r')
	for gene in gene_list:
		#print gene
		gene=gene.rstrip()
		if gene.startswith('CHROM'):
			continue
		else:
			geni+=[gene]

	#print geni
	prova=geni
	gene_list.close()

	trdb=open(opts.file,'r')
	
	for line in trdb:
		line=line.rstrip()
		if line.startswith('Gene'):
			header=line.split('\t')
			ind=header.index('Associated Gene Name')
			print line
			#print ind
		else:
			line_split=line.split('\t')
			#print line_split[ind]
			if line_split[ind] in geni:
					prova.remove(line_split[ind])
					print line
			

	trdb.close()
	if prova != []:
		sys.stderr.write("Sono rimasti fuori i seguenti geni: \n"+'\n'.join(prova))

main()