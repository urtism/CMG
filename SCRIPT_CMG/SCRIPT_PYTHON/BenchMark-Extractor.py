import argparse
import re
import string
import sys
from collections import Counter


#------------------------------------------FUNCTION---------------------------------

def extract_variant(line,ID_PAZ,DATA_PAZ,flist,variant_caller,CALLER,Header_input,debugfile):

	header_vcf = []
	var_info = {}

	with open(flist) as vcf:
		for row in vcf:

			#Estraggo l'header del vcf e strippo il # da #CHROM
			if row.startswith('#CHROM'):
				row = row.rstrip()
				header_vcf = row.split('\t')
				header_vcf[0] = header_vcf[0].lstrip('#')
				

			row = row.rstrip()
			row = row.split('\t')

			#Se la variante ha la stessa chrom pos ref alt allora la prendo ed estraggo la var del paziente
			if row[0].startswith('chr') and line[Header_input.index('CHR')] == row[header_vcf.index('CHROM')] \
				and line[Header_input.index('POS')] == row[header_vcf.index('POS')] \
				and line[Header_input.index('REF')] == row[header_vcf.index('REF')] \
				and line[Header_input.index('ALT')] == row[header_vcf.index('ALT')]:

				#print line

				#Ora scrollo tutte le keys per assegnare i valori trovati:
				for key in CALLER.keys():

					# Assegno alla key ID il nome del sample
					if key == 'ID':
						ind = header_vcf.index(ID_PAZ)
						CALLER[key] = header_vcf[ind]

						#Faccio giusto un controllo in piu per verificare che gli ID dei sample corrispondano
						if header_vcf[ind] != ID_PAZ:
							sys.exit('\n' + 'Codice sample in tsv e codice sample in vcf non corrispondono: Run: '\
								+ DATA_PAZ + ' Sample-tsv: ' + ID_PAZ + ' Sample-vcf: ' + header_vcf[ind] + '\n')

					#con il try estraggo prima i campi che sono unici nell'header del vcf come CHROM,POS,REF,ALT,QUAL
					else:
						try:
							CALLER[key] = row[header_vcf.index(key.replace(variant_caller+'-','').replace('iEVA-',''))]
						except:
							INFO = row[header_vcf.index('INFO')].split(';')
							FORMAT = row[header_vcf.index('FORMAT')].split(':')
							SAMPLE_FORMAT = row[header_vcf.index(ID_PAZ)].split(':')


							#Passo ad estrarre i valori di INFO
							for val in INFO:
								#N.B Dato che la DP e' ripetuta sia in INFO che FORMAT, controllo che nel file delle feature da estrarre ci siano INFO
								#prima della DP --> GATK-INFO-DP
								if (key == variant_caller+'-INFO-DP') and ('DP' in val.split('=')[0] and variant_caller == 'GATK' or variant_caller == 'FREE'):
									CALLER[key] = val.split('=')[1]

								#Ora controllo tutte le altre
								elif key.replace(variant_caller+'-','').replace('iEVA-','') == val.split('=')[0]:
									CALLER[key] = val.split('=')[1]
									#print key + ': ' + str(CALLER[key])


							#Passo ad estrarre i valori di FORMAT
							for val in FORMAT:
								if key.replace(variant_caller+'-','').replace('iEVA-','') == val:
									CALLER[key] = SAMPLE_FORMAT[FORMAT.index(val)]

				return CALLER





#-----------------------------------------MAIN----------------------------------


def main():
	
	parser = argparse.ArgumentParser('Questo tool estrae le varianti dai vcf a partire da una lista di varianti in tab delimited. \
		BenchMark-Extractor prende la variante (deve avere nel tab delimited i campi CHR POS REF ALT ID e CLASS \
		con ID riferito al codice paziente e CLASS riferito alla classe che puo essere PASS o FILTER (verificata in sanger o wt) \
		e va a controllare in tutti i file forniti nel list a quale run appartiene e poi estra l inter riga del \
		vcf in un nuovo file (UNO PER OGNI VARIANT CALLER) tab delimited con tutti i campi del vcf splittati correttamente. \
		PREREQUISITO: data in ID paziente uguale a data in ID run nel nome del file. Es: ID_PAZ = 20160724_01_Conn mentre \
		nome del file = 20150716_Cardio_iEVA_GATK.vcf --> Basta che 20160724 sia presente anche nel file.')
	parser.add_argument('-I','--input',help="File tab delimited contenente le varianti confermate in sanger")
	parser.add_argument('-L','--list',default=None,help="path contenente i file in cui cercare le varianti. N.B il codice paziente \
		contiene la data della run, per cui lo script accede direttamente a quel file.")
	parser.add_argument('-P','--panel',help="Specificare il pannello utilizzato nell'analisi: Cancer | Cardio | Exome | BRCA. \
		Serve per estrarre i file giusti nel caso in cui due run di diversi pannelli abbiano corso nello stesso giorno (quindi data \
		run uguale ma ID_PAZ diversi")
	parser.add_argument('-HG','--HeaderGATK',help="path al file contenente gli attributi da estrarre dal file vcf di GATK. Un attributo per riga. \
		CONTROLLA SE ESISTONO TAG UGUALI TRA INFO E FORMAT. Se ESISTONO, allora aggiungi INFO prima della tag in questo file. Es: GATK-DP (nel campo INFO)\
		diventa GATK-INFO-DP. Poi controlla la funzione extract_variant e aggiungi un if per valutare quell'info come ho gia fatto per la DP (basati su quell'esempio)")
	parser.add_argument('-HF','--HeaderFreeBayes',help="path al file contenente gli attributi da estrarre dal file vcf di FreeBayes. Un attributo per riga. \
		CONTROLLA SE ESISTONO TAG UGUALI TRA INFO E FORMAT. Se ESISTONO, allora aggiungi INFO prima della tag in questo file. Es: GATK-DP (nel campo INFO)\
		diventa FREE-INFO-DP. Poi controlla la funzione extract_variant e aggiungi un if per valutare quell'info come ho gia fatto per la DP (basati su quell'esempio)")
	parser.add_argument('-HV','--HeaderVarScan',help="path al file contenente gli attributi da estrarre dal file vcf di VarScan. Un attributo per riga. \
		CONTROLLA SE ESISTONO TAG UGUALI TRA INFO E FORMAT. Se ESISTONO, allora aggiungi INFO prima della tag in questo file. Es: GATK-DP (nel campo INFO)\
		diventa VARSCAN-INFO-DP. Poi controlla la funzione extract_variant e aggiungi un if per valutare quell'info come ho gia fatto per la DP (basati su quell'esempio)")
	parser.add_argument('-O','--outfile',help="file di output in vcf format. Indicare come path/to/nome_file --> In uscita si avra un file \
		per ciascun variant caller (e.g. nome_file_GATK.vcf) usato (prendo il nome del variant caller dal nome del file vcf).")
	parser.add_argument('-debug','--debugfile',help="Path del file di log per vedere eventuali errori.") 

	global opts

	run = []
	opts = parser.parse_args()
	out = opts.outfile
	pan = opts.panel
	save_head = []
	HG = []
	HF = []
	HV = []
	header_GATK = 0
	header_FreeBayes = 0
	header_VarScan = 0
	GATK = {}
	FreeBayes = {}
	VarScan = {}
	variant_caller = ''
	sample = {}
	out_GATK = out+'_GATK.tsv'
	out_FREE = out+'_FREE.tsv'
	out_VARSCAN = out+'_VARSCAN.tsv'

	print out_GATK

	#Salvo le header da estrarre presenti nei file in input
	with open(opts.HeaderGATK) as head:
		for elem in head:
			HG += [elem.rstrip()]

	#Salvo le header da estrarre presenti nei file in input
	with open(opts.HeaderFreeBayes) as head:
		for elem in head:
			HF += [elem.rstrip()]

	#Salvo le header da estrarre presenti nei file in input
	with open(opts.HeaderVarScan) as head:
		for elem in head:
			HV += [elem.rstrip()]

	#Qui salvo i path dei file vcf
	with open(opts.list) as vcf_path:
		for rep in vcf_path:
			run += [rep.rstrip()]

	with open(opts.input,'r') as tsv:
		for line in tsv:

			#Verifico che header contenga chr pos ref alt e id:
			if 'CHR' and 'POS' and 'REF' and 'ALT' and 'ID' in line :

				line = line.rstrip()
				Header_input = line.split('\t')

			else:

				line = line.rstrip()
				line = line.split('\t')
				#Estraggo l'ID del paziente
				ID_PAZ = line[Header_input.index('ID')]
				DATA_PAZ = line[Header_input.index('ID')].split('_')[0]

				#Applico un for alla lista dei vcf
				for flist in run:
					#Se il paziente ha la stessa data del file vcf e il variant caller e' GATK, allora passo ad estrarre
					if DATA_PAZ in flist and 'GATK' in flist:
						variant_caller = 'GATK'
						if header_GATK == 0:
							out2 = out + '_GATK.tsv'
							outfile = open(out2, 'w')
							outfile.write('\t'.join(HG) + '\n')
							outfile.close()

							#Costruisco un hash contenente come chiavi i tag da estrarre
							for tag in HG:
								GATK[tag] = '.'

							header_GATK = 1
						
						elif header_GATK == 1:
							
							variant = extract_variant(line,ID_PAZ,DATA_PAZ,flist,variant_caller,GATK,Header_input,opts.debugfile)

							open(out_GATK)

							#LUNEDI: SCRIVERE IL WRITE DELLA FUNZIONE

							#for hd in Header_input:


							#print variant

							
							#print GATK
							
							
						

					elif DATA_PAZ in flist and 'FreeBayes' in flist:
						variant_caller = 'FREE'
						if header_FreeBayes == 0:
							out3 = out + '_FreeBayes.tsv'
							outfile = open(out3, 'w')
							outfile.write('\t'.join(HF) + '\n')
							outfile.close()

							for tag in HF:
								#tag.replace('FREE-','').replace('iEVA-','')
								FreeBayes[tag] = '.'

							header_FreeBayes = 1

						elif header_FreeBayes == 1:
							continue




					elif DATA_PAZ in flist and 'VarScan' in flist:
						variant_caller = 'VARSCAN'
						if header_VarScan == 0:
							out4 = out + '_VarScan.tsv'
							outfile = open(out4, 'w')
							outfile.write('\t'.join(HV) + '\n')
							outfile.close()

							for tag in HV:
								VarScan[tag] = '.'


							header_VarScan = 1

						elif header_VarScan == 1:
							continue



						


		#ORA AGGIUNGO IL PEZZO DI CONTROLLO FILE-VAR E POS

	# 	#Vado a lavorare sull'header per estrarre le informazione e modificare il campo INFO
	# 	if line.startswith('#'):
	# 		line = line.rstrip()

	# 		if line.startswith('##INFO=<ID=ANN') or line.startswith('##INFO=<ID=CSQ'):
	# 			new_vep = add_INFO(line)
	# 			save_head.append(new_vep)
	# 			header_ann = extract_header(new_vep)

	# 		elif line.startswith('#CHROM'):
	# 			save_head.append(line)
	# 			out.write('\n'.join(save_head) + '\n')

	# 		else:
	# 			save_head.append(line)

	# 	else:
		
	# 		line = line.rstrip()
	# 		line = line.split('\t')

	# 		if opts.clinvar:
	# 			line = Add_ClinVar(line)

	# 		if opts.Humsavar:
	# 			line = Add_Humsavar(line,header_ann)

	# 		if opts.exome_variant_server:
	# 			line = Add_ESP(line)

	# 		if opts.gerp:
	# 			line = Add_Gerp(line,pan)

	# 		if opts.phastCons:
	# 			line = Add_phastCons(line,pan)

	# 		if opts.phyloP:
	# 			line = Add_phyloP(line,pan)

	# 		print line
	# 		out.write('\t'.join(line) + '\n')

	# vcf.close()


main()