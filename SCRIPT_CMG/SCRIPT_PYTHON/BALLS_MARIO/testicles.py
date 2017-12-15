import argparse
import os
import glob2


def main():


	
	parser = argparse.ArgumentParser('Cerca quali geni sono mutati nei pazienti presenti nella lista pazienti')
	parser.add_argument('-p', '--lista_pazienti', help="lista di pazienti contente CognomeNome	id_illumina",default=None)
	parser.add_argument('-d', '--tsv_dir', help="directory che contiene i tsv",default=None)
	parser.add_argument('-o', '--out', help="dir output",default=None)
	parser.add_argument('-l', '--lista_geni', help="Lista di geni",default=None)
	parser.add_argument('-m','--mut_escluse', help="mutazioni da escludere separati da -",default=None)
	parser.add_argument('--omim', help="Lista di geni con info da omim",default=None)
	parser.add_argument('--rs', help="Lista DI RS da cercare nei pazienti della lista_pazienti",default=None)



	
	global opts
	opts = parser.parse_args()
	
	if opts.rs:
		rs=open(opts.rs,'r')
		lista_paz=open(opts.lista_pazienti,'r')
		rs_list=dict()
		num_paz_testicoli=0
		num_paz_other=0
		for paz in lista_paz:
			paz=paz.rstrip()
			nome=paz.split('\t')[0]
			id=paz.split('\t')[1]
			tum=paz.split('\t')[-1]
			tsv_paz=open(opts.tsv_dir+'/'+id+'_GATK_Conferme.tsv','r')
			if tum=='SI':
				num_paz_testicoli+=1
			else:
				num_paz_other+=1
			for line in rs:
				line=line.rstrip()
				if line.startswith('RS'):
					continue
				else:
					rs_id=line.split('\t')[0]
					geni_vicini=line.split('\t')[-1]
					rs_list[rs_id]=[[],[],geni_vicini,0,0]
			
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
						existing_variation=line.split('\t')[header.index('Existing_variation')]
					except:
						continue
					for rs_id in rs_list.keys():
						print rs_id
						if rs_id in existing_variation:
							#print rs_id
							s=rs_list.get(rs_id)
							
							if tum=='SI':
								s[0]+=[nome]
								s[3]+=1
							else:
								s[1]+=[nome]
								s[4]+=1
							rs_list[rs_id]=s
		
		out=open(opts.out,'w')
		out.write('RS\tMAF_TUM_TESTICOLI\tMAF_ALTRI_TUM\tPAZIENTI_TUM_TESTICOLI\tPAZIENTI_ALTRI_TUM\tGENI_VICINI\n')
		for rs_id in rs_list.keys():

			if rs_list[rs_id][0]==[]:
				rs_list[rs_id][0]=['.']
			out.write(rs_id+'\t'+str(rs_list[rs_id][3])+'/'+str(num_paz_testicoli)+'\t'+str(rs_list[rs_id][4])+'/'+str(num_paz_other)+'\t'+'; '.join(rs_list[rs_id][0])+'\t' +'; '.join(rs_list[rs_id][1])+ '\t'+rs_list[rs_id][2]+'\n')
	else:
								
	
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
