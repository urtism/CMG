import argparse
import re

def add_ann(vcf,file_list,file_coor):
	out = open(opts.outfile,'w')
	file_list = open(opts.list,'r')
	file_coor = open(opts.file,'r')
	#Con questo leggo il file delle annotazioni da prendere e assegno ad ogni elemento del vettore uan riga letta
	ann_list = file_list.readlines()
	vettore = []
	
	for line in vcf:
		line = line.rstrip()
		if line.startswith('##INFO=<ID=ANN') or line.startswith('##INFO=<ID=CSQ'):
			
			#Con questo estraggo la header contenente INFO dell'annotazione, del tipo Allele|Consequence|IMPACT etc...
			#Con start e end definisco la parola da pescare in mezzo all'header
			start = line.find('Allele')
			end = line.find('">')
			#Con questo assegno ad header la stringa trovata da start a stop
			
			header_ann = line[start:end]
			header_ann = header_ann.split('|')
			#Se non metto continue mi salva in vettore[] anche la header che inizia con ##INFO=<ID=ANN
			continue

		elif line.startswith('##'):
			#Sto dicendo di andare avanti se incontra le altre header
			continue

		elif line.startswith('#CHROM'):
			header_chrom = line.split('\t')
			#Salvo solo l'header fino al campo INFO + le annotazioni
			header = header_chrom[0:8] + header_ann
			vettore = vettore + [header]
			continue

		else:
			line = re.split('\||ANN=|;\t|\t|,', line)
			for i in line:
				if i == '':
					#qui sostituisco al valore vuoto il simbolo -
					line[line.index(i)] = '-'
		# Salvo in ogni elemento del vettore il vettore contenente ciascuna riga del vcf annotato
		vettore += [line[0:len(header)]]

	# Con i cicli for di seguito si procede cosi:
	# 1) il primo for processa tutte le varianti nel file tab delimited
	# 2) il secondo processa tutte le varianti nel vettore contenente in ogni elemento le righe del vcf precedentemente modificate
	# 3) verifico che CHROM POS REF ALT siano uguali tra i due file elem[0].lstrip('#') == var[0] and elem[1] ecc..
	# 4) Se nella lista annotazioni leggi # prima della parola, allora salta quell'annotazione, altrimenti annota

	for var in file_coor:

		var = var.rstrip()
		var = var.split('\t')

		for elem in vettore:

			if elem[0].lstrip('#') == var[0] and elem[1] == var[1] and elem[3] == var[3] and elem[4] == var[4]:

				for tag in ann_list:
					
					tag = tag.rstrip()	
					
					if tag.startswith('#'):
						continue

					if tag not in header:
				#		print tag + 'non e nella lista'
						continue

				#	else:
					#tag = tag.rstrip()

					ann = header.index(tag)
					var.append(elem[ann])
					#print var
				out.write('\t'.join(var) + '\n')
			#else:
			#	out.write('\t'.join(var) + '\n' + 'Annotazione Non Trovata')




def main():
	parser = argparse.ArgumentParser('aggiunge le annotazioni fornite nel file -f (una per riga) al file di input -i in un formato tab delimited in output -o')
	
	parser.add_argument('-i','--input',help="file di input in formato vcf")
	parser.add_argument('-l','--list',help="lista di annotazioni: una per riga")
	parser.add_argument('-f','--file',help="file tab delimited da cui pescare le coordinate del cromosoma")
	parser.add_argument('-o','--outfile',help="file di output tab delimited")
	
	global opts
	
	opts = parser.parse_args()
	
	vcf = open(opts.input,'r')
	file_list = open(opts.list,'r')
	file_coor = open(opts.file,'r')
	add_ann(vcf,file_list,file_coor)

main()