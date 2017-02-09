import argparse
import os


def main():

	parser = argparse.ArgumentParser('aggiunge le annotazioni fornite nel file -f (una per riga) al file di input -i in un formato tab delimited in output -o')

	parser.add_argument('-p','--paz_list',help="lista di pazienti da cercare nei database")
	parser.add_argument('-d','--dbpath',help="path del DB diviso per geni")
	parser.add_argument('-f','--minfreq',help="soglia di frequenza allelica nel database per riportare una variante",default=1)
	parser.add_argument('-o','--out',help="path di output")

	global opts

	opts = parser.parse_args()
	paz_list = open(opts.paz_list,'r')
	for paz in paz_list:
		to_print=[]
		paz=paz.rstrip()
		paz_name='_'.join(paz.split(' '))
		header_to_print='\t'.join(["CHROM","POS","REF","ALT","GENE","HGVSc","HGVSp","CONSEQUENCE","FRAZ_ALLELICA","FRAZ_PERC","MAF","ETERO","OMO","NUM_PAZ_MUTATI","PAZIENTI"])
		for genedir in os.listdir(opts.dbpath):
			if 'Conferme' in genedir:
				continue
			else:
				genepath = '/'.join([opts.dbpath,genedir])
				for file in os.listdir(genepath):
					if '.xls' in file or 'Other' in file or 'GATK' in file:
						continue
					else:
						gene=open('/'.join([genepath,file]),'r')
						#print '/'.join([genepath,file])
						for line in gene:
							line = line.rstrip()
							line_split = line.split('\t')
							if line.startswith('CHROM'):
								header=line_split
							elif paz in line_split[-1].split(';') and float(line_split[header.index('MAF')]) < float(opts.minfreq): 
								line_split[header.index('ALT')+1:header.index('ALT')+1]=[genedir]
								to_print.append('\t'.join(line_split))
						gene.close()


		if to_print != []:

			out=open(opts.out+'/' + paz_name+ '.tsv','w')
			out.write(header_to_print + '\n')
			for elem in to_print:
				out.write(elem + '\n')

	paz_list.close()


main()