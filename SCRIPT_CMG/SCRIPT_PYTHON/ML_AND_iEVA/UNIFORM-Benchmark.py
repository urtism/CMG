import argparse
import re
import string
import sys

#------------------------------------------FUNCTION---------------------------------

def extract_FREE(New_Head,tag):

	if 'FREE-FORMAT-DPR' == tag:
		New_Head +=['FREE-FORMAT-DPR-TOT']
		New_Head +=['FREE-FORMAT-DPR-ALT']

	elif 'FREE-FORMAT-AD' == tag:
		New_Head += ['FREE-FORMAT-AD-REF']
		New_Head += ['FREE-FORMAT-AD-ALT']

	elif 'FREE-FORMAT-GL' == tag:
		New_Head += ['FREE-FORMAT-GL-00']
		New_Head += ['FREE-FORMAT-GL-01']
		New_Head += ['FREE-FORMAT-GL-02']

	else:
		New_Head += [tag]

	return New_Head

def extract_GATK(New_Head,tag):

 	if 'GATK-FORMAT-AD' == tag:
 		New_Head += ['GATK-FORMAT-AD-REF']
		New_Head += ['GATK-FORMAT-AD-ALT']

	elif 'GATK-FORMAT-PL' == tag:
		New_Head += ['GATK-FORMAT-GL-00']
		New_Head += ['GATK-FORMAT-GL-01']
		New_Head += ['GATK-FORMAT-GL-02']

	else:
		New_Head += [tag]

 	return New_Head

def extract_iEVA(New_Head,tag):

	if 'iEVA-INFO-PNC' == tag:
		New_Head += ['iEVA-INFO-PNC-AA']
		New_Head += ['iEVA-INFO-PNC-AC']
		New_Head += ['iEVA-INFO-PNC-AG']
		New_Head += ['iEVA-INFO-PNC-AT']
		New_Head += ['iEVA-INFO-PNC-CA']
		New_Head += ['iEVA-INFO-PNC-CC']
		New_Head += ['iEVA-INFO-PNC-CG']
		New_Head += ['iEVA-INFO-PNC-CT']
		New_Head += ['iEVA-INFO-PNC-GA']
		New_Head += ['iEVA-INFO-PNC-GC']
		New_Head += ['iEVA-INFO-PNC-GG']
		New_Head += ['iEVA-INFO-PNC-GT']
		New_Head += ['iEVA-INFO-PNC-TA']
		New_Head += ['iEVA-INFO-PNC-TC']
		New_Head += ['iEVA-INFO-PNC-TG']
		New_Head += ['iEVA-INFO-PNC-TT']

	elif 'iEVA-FORMAT-iAD' == tag:
		New_Head += ['iEVA-FORMAT-iAD-REF']
		New_Head += ['iEVA-FORMAT-iAD-ALT']

	elif 'iEVA-FORMAT-iADup' == tag:
		New_Head += ['iEVA-FORMAT-iADup-REF']
		New_Head += ['iEVA-FORMAT-iADup-ALT']

	elif 'iEVA-FORMAT-iQual' == tag:
		New_Head += ['iEVA-FORMAT-iQual-REF']
		New_Head += ['iEVA-FORMAT-iQual-ALT']

	elif 'iEVA-FORMAT-iMMQ' == tag:
		New_Head += ['iEVA-FORMAT-iMMQ-REF']
		New_Head += ['iEVA-FORMAT-iMMQ-ALT']

	elif 'iEVA-FORMAT-iAS' == tag:
		New_Head += ['iEVA-FORMAT-iAS-REF']
		New_Head += ['iEVA-FORMAT-iAS-ALT']

	elif 'iEVA-FORMAT-iXS' == tag:
		New_Head += ['iEVA-FORMAT-iXS-REF']
		New_Head += ['iEVA-FORMAT-iXS-ALT']

	elif 'iEVA-FORMAT-iXS0' == tag:
		New_Head += ['iEVA-FORMAT-iXS0-REF']
		New_Head += ['iEVA-FORMAT-iXS0-ALT']

	elif 'iEVA-FORMAT-iUnMap' == tag:
		New_Head += ['iEVA-FORMAT-iUnMap-REF']
		New_Head += ['iEVA-FORMAT-iUnMap-ALT']

	elif 'iEVA-FORMAT-iMQ0' == tag:
		New_Head += ['iEVA-FORMAT-iMQ0-REF']
		New_Head += ['iEVA-FORMAT-iMQ0-ALT']

	elif 'iEVA-FORMAT-iNPA' == tag:
		New_Head += ['iEVA-FORMAT-iNPA-REF']
		New_Head += ['iEVA-FORMAT-iNPA-ALT']

	elif 'iEVA-FORMAT-iQual' == tag:
		New_Head += ['iEVA-FORMAT-iQual-REF']
		New_Head += ['iEVA-FORMAT-iQual-ALT']

	elif 'iEVA-FORMAT-iSA' == tag:
		New_Head += ['iEVA-FORMAT-iSA-REF']
		New_Head += ['iEVA-FORMAT-iSA-ALT']

	elif 'iEVA-FORMAT-iNP' == tag:
		New_Head += ['iEVA-FORMAT-iNP-REF']
		New_Head += ['iEVA-FORMAT-iNP-ALT']

	elif 'iEVA-FORMAT-iNPP' == tag:
		New_Head += ['iEVA-FORMAT-iNPP-REF']
		New_Head += ['iEVA-FORMAT-iNPP-ALT']

	elif 'iEVA-FORMAT-iCR' == tag:
		New_Head += ['iEVA-FORMAT-iCR-REF']
		New_Head += ['iEVA-FORMAT-iCR-ALT']

	else:
		New_Head += [tag]

	return New_Head





#-----------------------------------------MAIN----------------------------------


def main():
	
	parser = argparse.ArgumentParser('\n\nQuesto tool splitta e modifica il tsv ottenuto dal tool Bencharmk extractor. Quindi: splitta i campi di iEVA quali iAD, iDUP, ecc... in, ad esempio, iAD-REF e iAD-ALT etc... poi aggiusta altri campi come: il valore FREQ in VarScan ecc. Nel codice viene spiegato con commenti.\n')
	parser.add_argument('-I','--input',help="File tab delimited contenente le varianti estratte da Benchmark-Extractor. Per funzionare, il tool DEVE AVERE IL NOME DEL VARIANT CALLER (VarScan, FreeBayes, GATK, indipendentemente da lettere maiuscole o minuscole) nel nome del file. Indicare come PATH/TO/NOME_FILE")
	parser.add_argument('-O','--outfile',help="File di output in TSV format. Se opzione -W attiva, allora l'output sara' in formato CSV. Indicare come path/to/nome_file.tsv")
	parser.add_argument('-W','--weka',action="store_true",help="se attivi questa opzione, in output avrai la sostituzione dei valori con . in ? come utile nel formato per weka")

	global opts

	opts = parser.parse_args()
	out = open(opts.outfile,'w')
	New_Head = []

	with open(opts.input,'r') as tsv:
		for line in tsv:

			#Verifico che header contenga chr pos ref alt e id:
			if 'CHROM' and 'POS' and 'REF' and 'ALT' and 'ID' in line:

				#Salvo l'header
				line = line.rstrip()
				Header_input = line.split('\t')

				for tag in Header_input:

					#Aggiungo le tag che splitto in FREE nella head
					if 'FreeBayes'.lower() in opts.input.lower():

						if 'FREE' in tag:
							New_Head = extract_FREE(New_Head,tag)

						elif 'iEVA' in tag:
							New_Head = extract_iEVA(New_Head,tag)

						else:
							New_Head += [tag]


					#Aggiungo le tag che splitto in VARSCAN nella head
					if 'VarScan'.lower() in opts.input.lower():

						if 'iEVA' in tag:
							New_Head = extract_iEVA(New_Head,tag)

						else:
							New_Head += [tag]


					#Aggiungo le tag che splitto in GATK nella head
					if 'GATK'.lower() in opts.input.lower():

						if 'GATK' in tag:
							New_Head = extract_GATK(New_Head,tag)

						elif 'iEVA' in tag:
							New_Head = extract_iEVA(New_Head,tag)

						else:
							New_Head += [tag]

				if opts.weka:
						
					out.write(','.join(New_Head) + '\n')

				if opts.weka == None:

					out.write('\t'.join(New_Head) + '\n')


			#Altrimenti vado a modificare i dati
			else:

				line = line.rstrip()
				line = line.split('\t')


				if 'FreeBayes'.lower() in opts.input.lower():

					if line[Header_input.index('FREE-FORMAT-GL')] == '.':
						line[Header_input.index('FREE-FORMAT-GL')] = '.,.,.'

					if line[Header_input.index('FREE-FORMAT-AD')] == '.':
						line[Header_input.index('FREE-FORMAT-AD')] = '.,.'

					if line[Header_input.index('FREE-FORMAT-DPR')] == '.':
						line[Header_input.index('FREE-FORMAT-DPR')] = '.,.'

					new_line = [item for val in line for item in val.split(',')]



				if 'VarScan'.lower() in opts.input.lower():

					#Scrivo la freq in in decimi, tolgo la virgola e poi il simbolo %
					try:
						line[Header_input.index('VARSCAN-FORMAT-FREQ')] = str(float(line[Header_input.index('VARSCAN-FORMAT-FREQ')].rstrip('%').replace(',','.'))/float(100))
					except:
						pass

					#Sosituisco la virgola col punto nel campo PVAL
					try:
						line[Header_input.index('VARSCAN-FORMAT-PVAL')] = line[Header_input.index('VARSCAN-FORMAT-PVAL')].replace(',','.')
					except:
						pass

					new_line = [item for val in line for item in val.split(',')]



				if 'GATK'.lower() in opts.input.lower():

					new_line = [item for val in line for item in val.split(',')]


				#Se attivo opzione WEKA allora output in formato csv
				if opts.weka:
					for i in range(0,len(new_line)):
						if new_line[i] == '.':
							new_line[i] = '?'

					out.write(','.join(new_line) + '\n')

				if opts.weka == None:

					out.write('\t'.join(new_line) + '\n')

	out.close()

main()