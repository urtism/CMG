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
	parser.add_argument('-F','--format',default="tsv",action='store',type=str,choices=['tsv','csv'],help="Specificare il formato del file in uscita e in ingresso. valori possibili: csv e tsv")


	global opts

	opts = parser.parse_args()

	if opts.format.lower() == 'csv':
		out = open(opts.outfile+'.csv','w')
		out2 = open(opts.outfile+'_out_target.csv','w')
	elif opts.format.lower() == 'tsv':
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

			if opts.format.lower() == 'csv':
				
				#Scrivo l'header in output
				if 'CHROM' and 'POS' and 'REF' and 'ALT' and 'ID' in line:
					line = line.rstrip()
					line = line.split(',')
					input_Header = line
					out2.write(','.join(input_Header) + '\n')
					out.write(','.join(input_Header) + '\n')

				
				else:
					line = line.rstrip()
					line = line.split(',')

					for targ in target_list:

						#Verifico che la variante sia compresa nel target:
						if line[input_Header.index('CHROM')] == targ[0] and line[input_Header.index('POS')] >= targ[1] \
						and line[input_Header.index('POS')] <= targ[2]:


							print str(line[input_Header.index('CHROM')]) + ' ' + str(targ[1]) + ' ' + str(line[input_Header.index('POS')]) + ' ' + str(targ[2])


							out.write(','.join(line) + '\n')
							variant = 1

							break


						else:
							variant = 0
							continue

					if variant == 0:

						print 'Variant out of target: ' + str(line[input_Header.index('CHROM')]) + ' ' + str(line[input_Header.index('POS')])
						out2.write(','.join(line) + '\n')


			elif opts.format.lower() == 'tsv':

				#Scrivo l'header in output
				if 'CHROM' and 'POS' and 'REF' and 'ALT' and 'ID' in line:
					line = line.rstrip()
					line = line.split('\t')
					input_Header = line
					out2.write('\t'.join(input_Header) + '\n')
					out.write('\t'.join(input_Header) + '\n')

				else:
					line = line.rstrip()
					line = line.split(',')

					for targ in target_list:

						#Verifico che la variante sia compresa nel target:
						if line[input_Header.index('CHROM')] == targ[0] and line[input_Header.index('POS')] >= targ[1] \
						and line[input_Header.index('POS')] <= targ[2]:

							out.write('\t'.join(line) + '\n')
							variant = 1

							break


						else:
							continue

					if variant == 0:

						print 'Variant out of target: ' + str(line[input_Header.index('CHROM')]) + ' ' + str(line[input_Header.index('POS')])
						out2.write('\t'.join(line) + '\n')


			else:
				sys.exit('\n\nNon esiste il formato ' + str(opts.format) + '. Guarda la opzione -F/--format per informazioni\n')

	out.close()
	out2.close()


main()