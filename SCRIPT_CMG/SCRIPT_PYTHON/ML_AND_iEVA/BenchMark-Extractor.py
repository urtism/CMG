import argparse
import re
import string
import sys

#------------------------------------------FUNCTION---------------------------------

def extract_variant(line,ID_PAZ,DATA_PAZ,flist,variant_caller,CALLER,Header_input,Header_line):

	header_vcf = []

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
							CALLER[key] = row[header_vcf.index(key)]
							#.replace(variant_caller+'-','').replace('iEVA-','')
						except:
							INFO = row[header_vcf.index('INFO')].split(';')
							FORMAT = row[header_vcf.index('FORMAT')].split(':')
							SAMPLE_FORMAT = row[header_vcf.index(ID_PAZ)].split(':')


							#Passo ad estrarre i valori di INFO
							for val in INFO:
								#N.B Dato che la DP e' ripetuta sia in INFO che FORMAT, controllo che nel file delle feature da estrarre ci siano INFO
								#prima della DP --> GATK-INFO-DP. Idem per AO in FreeBayes
								if 'INFO' in key and \
								(key.replace(variant_caller+'-INFO-','').replace('iEVA-INFO-','') == val.split('=')[0]):
									CALLER[key] = val.split('=')[1]

							#Passo ad estrarre i valori di FORMAT
							for val in FORMAT:
								if 'FORMAT' in key and key.replace(variant_caller+'-FORMAT-','').replace('iEVA-FORMAT-','') == val:
									CALLER[key] = SAMPLE_FORMAT[FORMAT.index(val)]


				return CALLER





#-----------------------------------------MAIN----------------------------------


def main():
	
	parser = argparse.ArgumentParser('\n\nQuesto tool estrae le varianti dai vcf a partire da una lista di varianti in tab delimited. BenchMark-Extractor prende la variante (deve avere nel tab delimited i campi CHR POS REF ALT ID e CLASS con ID riferito al codice paziente e CLASS riferito alla classe che puo essere PASS o FILTER (verificata in sanger o wt) e va a controllare in tutti i file forniti nel list a quale run appartiene e poi estra l inter riga del vcf in un nuovo file (UNO PER OGNI VARIANT CALLER) tab delimited con tutti i campi del vcf splittati correttamente. PREREQUISITO: data in ID paziente uguale a data in ID run nel nome del file. Es: ID_PAZ = 20160724_01_Conn mentre nome del file = 20150716_Cardio_iEVA_GATK.vcf --> Basta che 20160724 sia presente anche nel file.\n')
	parser.add_argument('-I','--input',help="File tab delimited contenente le varianti confermate in sanger")
	parser.add_argument('-L','--list',default=None,help="path contenente i file in cui cercare le varianti. N.B il codice paziente contiene la data della run, per cui lo script accede direttamente a quel file.")
	parser.add_argument('-HG','--HeaderGATK',help="path al file contenente gli attributi da estrarre dal file vcf di GATK. Un attributo per riga. CONTROLLA SE ESISTONO TAG UGUALI TRA INFO E FORMAT. Se ESISTONO, allora aggiungi INFO prima della tag in questo file. Es: GATK-DP (nel campo INFO) diventa GATK-INFO-DP. Poi controlla la funzione extract_variant e aggiungi un if per valutare quell'info come ho gia fatto per la DP (basati su quell'esempio)")
	parser.add_argument('-HF','--HeaderFreeBayes',help="path al file contenente gli attributi da estrarre dal file vcf di FreeBayes. Un attributo per riga. CONTROLLA SE ESISTONO TAG UGUALI TRA INFO E FORMAT. Se ESISTONO, allora aggiungi INFO prima della tag in questo file. Es: GATK-DP (nel campo INFO) diventa FREE-INFO-DP. Poi controlla la funzione extract_variant e aggiungi un if per valutare quell'info come ho gia fatto per la DP (basati su quell'esempio)")
	parser.add_argument('-HV','--HeaderVarScan',help="path al file contenente gli attributi da estrarre dal file vcf di VarScan. Un attributo per riga. CONTROLLA SE ESISTONO TAG UGUALI TRA INFO E FORMAT. Se ESISTONO, allora aggiungi INFO prima della tag in questo file. Es: GATK-DP (nel campo INFO) diventa VARSCAN-INFO-DP. Poi controlla la funzione extract_variant e aggiungi un if per valutare quell'info come ho gia fatto per la DP (basati su quell'esempio)")
	parser.add_argument('-O','--outfile',help="file di output in vcf format. Indicare come path/to/nome_file --> In uscita si avra un file per ciascun variant caller (e.g. nome_file_GATK.vcf) usato (prendo il nome del variant caller dal nome del file vcf).")
	parser.add_argument('-debug','--debugfile',help="Path del file di log per vedere eventuali errori.")

	global opts

	run = []
	opts = parser.parse_args()
	out = opts.outfile
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
	debug = open(opts.debugfile,'w')
	variant = {}


	#Salvo le header da estrarre presenti nei file in input
	with open(opts.HeaderGATK) as head:
		for elem in head:
			HG += [elem.rstrip()]
	HG = HG + ['GT_CLASS'] + ['CLASS']

	#Salvo le header da estrarre presenti nei file in input
	with open(opts.HeaderFreeBayes) as head:
		for elem in head:
			HF += [elem.rstrip()]
	HF = HF + ['GT_CLASS'] + ['CLASS']

	#Salvo le header da estrarre presenti nei file in input
	with open(opts.HeaderVarScan) as head:
		for elem in head:
			HV += [elem.rstrip()]
	HV = HV + ['GT_CLASS'] + ['CLASS']

	#Qui salvo i path dei file vcf
	with open(opts.list) as vcf_path:
		for rep in vcf_path:
			run += [rep.rstrip()]

	with open(opts.input,'r') as tsv:
		for line in tsv:

			print line

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
						#Scrivo l'header del file in uscita
						if header_GATK == 0:
							out2 = out + '_GATK.tsv'
							outfile = open(out2, 'w')
							outfile.write('\t'.join(HG) + '\n')

							#Costruisco un hash contenente come chiavi i tag da estrarre
							for tag in HG:
								GATK[tag] = '.'

							header_GATK = 1
						
						#Se la riga non e' l'header del file di input, allora estraggo la variante
						elif header_GATK == 1:
							
							variant = extract_variant(line,ID_PAZ,DATA_PAZ,flist,variant_caller,GATK,Header_input,HG)

							#Se la variante non e' stata trovata in GATK, allora lo segnalo scrivendolo nel file debug
							if variant is None:
								debug.write(variant_caller + ' ---> Missing variant for sample ' + line[1] + ' at position:\t' + line[2] + '-' + line[3] + '-' + line[4] + '-' + line [5] + '\n')

							#Altrimenti scrivo la variante nel file di output
							else:

								extracted_var = []

								for tag_head in HG:

									if tag_head == 'GT_CLASS':
										extracted_var += [line[Header_input.index(tag_head)]]
									elif tag_head == 'CLASS':
										extracted_var += [line[Header_input.index(tag_head)]]
									else:
										extracted_var += [variant.get(tag_head)]
								
								outfile.write('\t'.join(extracted_var) + '\n')
								extracted_var = []
								for key in GATK.keys():
									GATK[key] = '.'


					elif DATA_PAZ in flist and 'FreeBayes' in flist:
						variant_caller = 'FREE'
						if header_FreeBayes == 0:
							out3 = out + '_FreeBayes.tsv'
							outfile2 = open(out3, 'w')
							outfile2.write('\t'.join(HF) + '\n')

							for tag in HF:
								FreeBayes[tag] = '.'

							header_FreeBayes = 1

						elif header_FreeBayes == 1:
							
							variant = extract_variant(line,ID_PAZ,DATA_PAZ,flist,variant_caller,FreeBayes,Header_input,HG)

							#Se la variante non e' stata trovata in GATK, allora lo segnalo scrivendolo nel file debug
							if variant is None:
								debug.write(variant_caller + ' ---> Missing variant for sample ' + line[1] + ' at position:\t' + line[2] + '-' + line[3] + '-' + line[4] + '-' + line [5] + '\n')

							#Altrimenti scrivo la variante nel file di output
							else:

								extracted_var = []

								for tag_head in HF:

									if tag_head == 'GT_CLASS':
										extracted_var += [line[Header_input.index(tag_head)]]
									elif tag_head == 'CLASS':
										extracted_var += [line[Header_input.index(tag_head)]]
									else:
										extracted_var += [variant.get(tag_head)]
										FreeBayes[tag] = '.'
								
								outfile2.write('\t'.join(extracted_var) + '\n')
								extracted_var = []
								for key in FreeBayes.keys():
									FreeBayes[key] = '.'


					elif DATA_PAZ in flist and 'VarScan' in flist:
						variant_caller = 'VARSCAN'
						if header_VarScan == 0:
							out4 = out + '_VarScan.tsv'
							outfile3 = open(out4, 'w')
							outfile3.write('\t'.join(HV) + '\n')

							for tag in HV:
								VarScan[tag] = '.'

							header_VarScan = 1

						elif header_VarScan == 1:
							
							variant = extract_variant(line,ID_PAZ,DATA_PAZ,flist,variant_caller,VarScan,Header_input,HG)

							#Se la variante non e' stata trovata in GATK, allora lo segnalo scrivendolo nel file debug
							if variant is None:
								debug.write(variant_caller + ' ---> Missing variant for sample ' + line[1] + ' at position:\t' + line[2] + '-' + line[3] + '-' + line[4] + '-' + line [5] + '\n')

							#Altrimenti scrivo la variante nel file di output
							else:

								extracted_var = []

								for tag_head in HV:

									if tag_head == 'GT_CLASS':
										extracted_var += [line[Header_input.index(tag_head)]]
									elif tag_head == 'CLASS':
										extracted_var += [line[Header_input.index(tag_head)]]
									else:
										extracted_var += [variant.get(tag_head)]
								
								outfile3.write('\t'.join(extracted_var) + '\n')
								extracted_var = []
								for key in VarScan.keys():
									VarScan[key] = '.'

				#if check == 0:
				#	debug.write('Variant missing: --> Sample: ' + str(line[Header_input.index('ID')]) + ' at position ' + str(line[Header_input.index('CHR')]) + '-' + str(line[Header_input.index('POS')]) + '-' + str(line[Header_input.index('REF')]) + '-' + str(line[Header_input.index('ALT')]) + '\n')
				#	continue
				extracted_var = []
				variant={}
			debug.write('\n')


		outfile.close()
		outfile2.close()
		outfile3.close()
		debug.close()


main()