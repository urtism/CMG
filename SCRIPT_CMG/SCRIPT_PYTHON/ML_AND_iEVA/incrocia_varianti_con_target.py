import argparse
import re
import string
import sys


#-----------------------------------------MAIN------------------------------------


def main():
	
	parser = argparse.ArgumentParser('\n\nQuesto tool estrae le varianti che sono in target nel manifest dato in ingresso. Quelle escluse verranno riportate a parte in un altro file.\n')
	parser.add_argument('-I','--input',help="File tab delimited di input contenente le varianti (implementeremo vcf piu in la). REQUEST: Deve esserci CHROM e POS")
	parser.add_argument('-O','--outfile',help="File di output sempre in tsv format. Indicare come path/to/nome_file. In uscita avremo due file: path/to/nome_file e path/to/nome_file_excluded. Per l'estensione, leggi comando -F o --format")
	parser.add_argument('-T','--target',help="File di target. Sia in formato list (1-based) che bed (0-based)")


	global opts

	opts = parser.parse_args()

	out = open(opts.outfile+'.tsv','w')
	out2 = open(opts.outfile+'_out_target.tsv','w')

	target_list = []
	input_Header = []
	variant = 0

	#Salvo il target prima di tutto
	with open(opts.target,'r') as tsv:
		for line in tsv:

			if '.list' in opts.target:
				if line.startswith('chr'):
					line = line.rstrip()
					line = line.split('\t')
					target_list += [line]

			elif '.bed' in opts.target:
				if line.startswith('chr'):
					line = line.rstrip()
					line = line.split('\t')
					target_list += [line]


			else:
				continue


	with open(opts.input,'r') as var:
		for line in var:
			#Scrivo l'header in output
			if 'VAR_ID' in line:
				line = line.rstrip()
				line = line.split('\t')
				input_Header = line
				out2.write('\t'.join(input_Header) + '\n')
				out.write('\t'.join(input_Header) + '\n')

			elif 'ID' in line:
				line = line.rstrip()
				line = line.split('\t')
				input_Header = line
				out2.write('\t'.join(input_Header) + '\n')
				out.write('\t'.join(input_Header) + '\n')

			else:
				line = line.rstrip()
				line = line.split('\t')

				for targ in target_list:
					try:
						chrom,pos=line[input_Header.index('ID')].split('-')[1].split(':')
					except:
						chrom,pos=line[input_Header.index('VAR_ID')].split('-')[0].split(':')

					#Verifico che la variante sia compresa nel target:
					if chrom == targ[0] and pos >= targ[1] \
					and pos <= targ[2]:

						out.write('\t'.join(line) + '\n')
						variant = 1

						break


					else:
						continue

				if variant == 0:
					try:
						print 'Variant out of target: ' + line[input_Header.index('ID')]
					except:
						print 'Variant out of target: ' + line[input_Header.index('VAR_ID')]
					out2.write('\t'.join(line) + '\n')



	out.close()
	out2.close()

main()