import argparse
import re

def add_hum(tab,hum):

	out = open(opts.outfile,'w')

	#Con questo leggo il file delle annotazioni da prendere e assegno ad ogni elemento del vettore uan riga letta

	vettore = []
	vect = []


	for line in tab:

		line = line.rstrip()

		if line.startswith('CHROM'):

			header = line + '\tVariant_type'
			header = header.split('\t')
			vettore += [header]

		else:

			line = line.split('\t')
			vettore += [line]

	
	for var in hum:

		var = var.rstrip()

		if var.startswith('Main_gene_name'):
			header_var = var.split('\t')

		else:

			var = var.split('\t')
			vect += [var]

	for raw in vettore:
		if raw[0] == 'CHROM':
			out.write('\t'.join(raw) + '\n')
			continue
		else:
			a='-'
			for hum in vect:
				if hum[header_var.index('Main_gene_name')] in raw[header.index('SYMBOL')] and hum[header_var.index('AA_Change')] in raw[header.index('HGVSp')]:
					a=hum[header_var.index('Variant_type')]
					break
				else:
					continue
		out.write('\t'.join(raw+[a]) + '\n')


			
def main():
	parser = argparse.ArgumentParser('aggiunge le annotazioni fornite nel file -f (una per riga) al file di input -i in un formato tab delimited in output -o')
	
	parser.add_argument('-i','--input',help="file di input annotato contente le varianti")
	parser.add_argument('-f','--file',help="file tab delimited humsavar da cui matchare polimorfismi e varianti disease causing")
	parser.add_argument('-o','--outfile',help="file di output tab delimited: concatenazione dei tag humsavar sulla variante contenuta nel file --input")
	
	global opts
	
	#Con questo script vado ad aggiungere al file tab delimited annotato l'evebntuale presenza della variante nella lista di annotazioni disease causing
	# o polimorphism 

	opts = parser.parse_args()
	
	tab = open(opts.input,'r')
	hum = open(opts.file,'r')

	zzadd_hum(tab,hum)

main()