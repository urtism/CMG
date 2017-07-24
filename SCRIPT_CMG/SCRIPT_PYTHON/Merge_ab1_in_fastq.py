import argparse
import os
from Bio import SeqIO

def convert(file_a,tipo_a,file_b,tipo_b):
	SeqIO.convert(file_a, tipo_a,file_b, tipo_b)

def merge_ab1_in_fastq(list,fastq):
	for ab in list:
		nome_file=ab.split('/')[-1]
		if 'fastq' not in nome_file and 'phd' not in nome_file:
			convert(ab,"abi",ab+'.fastq',"fastq-illumina")
			ab_file=open(ab+'.fastq','r')
			for line in ab_file:
				fastq.write(line)


def merge_phd_in_fastq(list,fastq):
	for phd in list:
		nome_file=phd.split('/')[-1]
		if 'fastq' not in nome_file and 'ab' not in nome_file:
			convert(phd,"phd",phd+'.fastq',"fastq-illumina")
			phd_file=open(phd+'.fastq','r')
			for line in phd_file:
				fastq.write(line)


def main():
	
	parser = argparse.ArgumentParser('Genera un fastq mergiando i file ab1 presenti in una lista o in una dir')
	parser.add_argument('-d', '--dir', help="directory dove sono contenuti i file ba1 da mergiare in fastq",default=None)
	parser.add_argument('-l', '--list', help="list dei file ab1 da mergiare in fastq",default=None)
	parser.add_argument('-t', '--tipo', help="tipo di file in entrata phd o abi")
	parser.add_argument('-o', '--out', help="file di output")
	

	global opts 
	opts = parser.parse_args()
	fastq=open(opts.out,'w')
	if opts.dir:
		ab_list=(opts.dir+'/'+dir for dir in os.listdir(opts.dir))
		#print os.listdir(opts.dir)
	elif opts.list:
		ab_list=(open(opts.list)).readlines()
		#for subdir in os.listdir(opts.dir):
		#	list=os.listdir(subdir)
	#ab_list_sort=ab_list.sort()
	
	if opts.tipo == 'abi':
		merge_ab1_in_fastq(ab_list,fastq)
	elif opts.tipo == 'phd':
		merge_phd_in_fastq(ab_list,fastq)
	fastq.close()

main()



