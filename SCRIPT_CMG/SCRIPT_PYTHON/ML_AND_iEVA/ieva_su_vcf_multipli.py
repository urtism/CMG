import argparse
import os
import subprocess

path_ieva = '/home/jarvis/git/iEVA/iEVA/iEVA.py'
ref = '/home/jarvis/NGS_TOOLS/hg19/ucsc.hg19.fasta'



if __name__ == '__main__':

	parser = argparse.ArgumentParser('ANNOTA I VCF CONTENUTI NELLA LISTA CON IEVA E CREA LE CARTELLE PER DIVIDERE LE RUN')
	parser.add_argument('-l', '--list', help="list of vcf")
	parser.add_argument('-bl', '--bamlist', help="bam list")
	parser.add_argument('-vl', '--variantlist', help="variant list")

	parser.add_argument('-o', '--out', help="CARTELLA DI OUTPUT DOVE CREARE LE CARTELLE DELLE RUN E STORARE I VCF ANNOTATI")

	opts = parser.parse_args()

	vcflist = open(opts.list,'r')

	if not os.path.exists(opts.out):
			os.makedirs(opts.out)

	# bamdict= dict()

	# variantlist = open(opts.variantlist,'r')

	# bamlist = open(opts.bamlist,'r')


	# for line in variantlist:
	# 	line = line.rstrip()
	# 	samplename = line.split('\t')[1]
	# 	bamdict[samplename] = ''
	
	# for line in bamlist:
	# 	line = line.rstrip()
	# 	bamname = line.split('/')[-1].split('.')[0]
	# 	if bamname in bamdict.keys():
	# 		print line
		
	for vcf in vcflist:
		vcf=vcf.rstrip()
		
		cartella_run = opts.out + '/' + vcf.split('/')[-2]
		nome_vcf = vcf.split('/')[-1]
		nuovo_vcf = '.'.join([nome_vcf.split('.')[0],'iEVA',nome_vcf.split('.')[-1]])
		
		if not os.path.exists(cartella_run):
			os.makedirs(cartella_run)

		args = ['python',path_ieva,'-Ref',ref,'-I',vcf,'-O',cartella_run+ '/' + nuovo_vcf,'-v']

		args += ['-AS']
		args += ['-AB','-L',opts.bamlist]
		args += ['-AG']

		success = subprocess.call(args)