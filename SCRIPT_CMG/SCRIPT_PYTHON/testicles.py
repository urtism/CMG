import argparse
import os
import glob2


def main():


	
	parser = argparse.ArgumentParser('Cerca quali geni sono mutati nei pazienti presenti nella lista pazienti')
	parser.add_argument('-p', '--lista_pazienti', help="lista di pazienti contente CognomeNome	id_illumina")
	parser.add_argument('-d', '--tsv_dir', help="directory che contiene i tsv")
	parser.add_argument('-o', '--out', help="dir output")
	parser.add_argument('-l', '--lista_geni', help="Lista di geni")
	parser.add_argument('-m','--mut_escluse', help="mutazioni da escludere separati da -")
	parser.add_argument('--omim', help="Lista di geni con info da omim")


	
	global opts
	opts = parser.parse_args()
	
	
	geni = dict()
	mut_escluse= opts.mut_escluse.split('-')
	lista_geni=open(opts.lista_geni,'r')
	
	for line in lista_geni:
		gene=line.rstrip()
		file_gene=open(opts.out+'/'+gene+'.tsv','a')
		file_gene.write('CHROM\tPOS\tID\tREF\tALT\tPAZIENTE\tCONSEQUENCE\tHGVSc\tHGVSp\tGMAF\tGT\tAD\n')
		file_gene.close()
		geni[gene]=[0,[],'.','.']

	lista_geni.close()

	omim=open(opts.omim,'r')
	for line in omim:
		line=line.rstrip()
		gene=line.split('\t')[0]
		cod=line.split('\t')[2]
		try:
			phen=line.split('\t')[3]
		except:
			print line
			phen='.'
		if gene in geni.keys():
			dati=geni.get(gene)
			dati[2]=cod
			dati[3]=phen
			geni[gene]=dati
	omim.close()

	lista_paz=open(opts.lista_pazienti,'r')
	num_paz=0
	for paz in lista_paz:
		num_paz+=1
		paz=paz.rstrip()
		nome=paz.split('\t')[0]
		id=paz.split('\t')[1]

		tsv_paz=open(opts.tsv_dir+'/'+id+'_GATK_Conferme.tsv','r')
		for line in tsv_paz:
			line=line.rstrip()
			if line.startswith('CHROM'):
				header=line.split('\t')
			else:
				chrom=line.split('\t')[header.index('CHROM')]
				pos=line.split('\t')[header.index('POS')]
				ref=line.split('\t')[header.index('REF')]
				alt=line.split('\t')[header.index('ALT')]
				gene=line.split('\t')[header.index('SYMBOL')]
				cons=line.split('\t')[header.index('Consequence')]
				try:
					HGVSc=line.split('\t')[header.index('HGVSc')]
					HGVSp=line.split('\t')[header.index('HGVSp')]
					GMAF=line.split('\t')[header.index('GMAF')]
					gt=line.split('\t')[header.index('GT')]
					ad=line.split('\t')[header.index('AD')]
				except:
					print line

				if gene in geni.keys() and cons not in mut_escluse:
					file_gene=open(opts.out+'/'+gene+'.tsv','a')					
					file_gene.write(chrom+'\t'+pos+'\t'+id+'\t'+ref+'\t'+alt+'\t'+nome+'\t'+cons+'\t'+HGVSc+'\t'+HGVSp+'\t'+GMAF+'\t'+gt+'\t'+ad+'\n')

					dati=geni.get(gene)
					if nome in dati[1]:
						continue
					else:
						dati[0]+=1
						dati[1]+=[nome]
					geni[gene]=dati
		tsv_paz.close()
	lista_paz.close()


	out=open(opts.out+'/GENI_IN_COMUNE.tsv','w')
	out.write('GENE\tNUM_PAZ\tPAZIENTI\tCODIFICA\tFENOTIPO\n')
	
	for gene in geni.keys():
		dati=geni.get(gene)
		if dati[0]==0:
			continue
		else:
			out.write(gene+'\t'+str(dati[0])+'/'+str(num_paz)+'\t'+';'.join(dati[1])+'\t'+dati[2]+'\t'+dati[3]+'\n')


		
main()
