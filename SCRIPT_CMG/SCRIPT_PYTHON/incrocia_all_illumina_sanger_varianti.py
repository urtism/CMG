import argparse
import os




def main():
	parser = argparse.ArgumentParser('Adds number of exon for each interval')
	
	parser.add_argument('-i','--ill_path',help="path del file di varianti da illumina")
	parser.add_argument('-s','--sang_path',help="path del file di varianti da sanger")
	parser.add_argument('-l','--paz_list',help="file lista pazienti illumina con sc id")
	parser.add_argument('-o','--out',help="path di output")

	global opts
	opts = parser.parse_args()
	lista_snp=dict()

	lista_pazienti=open(opts.paz_list,'r')
	
	paz=dict()
	for line in lista_pazienti:
		line=line.rstrip()
		if line.startswith('SC_ID'):
			header=line.split('\t')
		else:
			sc_id=line.split('\t')[header.index('SC_ID')]
			nome=line.split('\t')[header.index('NOME_PAZ')]
			id=line.split('\t')[header.index('ID')]
			if sc_id != '.':
				paz[sc_id]=id

	for dir_sample in os.listdir(opts.sang_path):
		lista_snp=dict()
		sc_id_sang=dir_sample.split('.')[0]

		if sc_id_sang in paz.keys():
			id=paz[sc_id_sang]
			print sc_id_sang,id
			try:
				file1=open(opts.ill_path+'/'+id+'_GATK_Conferme.tsv','r')
			except:
				try:
					file1=open(opts.ill_path+'/'+id+'_GATK_Conferme.ods','r')
				except:
					print opts.ill_path+'/'+id+'_GATK_Conferme.ods','file non trovato'
					continue
			
			file2=open(opts.sang_path+'/'+sc_id_sang+'.tsv','r')
			inters=open(opts.out+'/'+id+'.intersect.tsv','w')
			inters_all=open(opts.out+'/ALL.intersect.tsv','a+')
			sanger=open(opts.out+'/'+id+'.sanger.tsv','w')
			#inters.write('\t'.join(['CHROM','POS','ID','REF','ALT'])+'\n')
			sanger.write('\t'.join(['CHROM','POS','REF','ALT'])+'\n')

			for snp in file1:
				snp=snp.rstrip()
				if snp.startswith('CHROM'):
					header=snp.split('\t')
					inters.write(snp +'\n')
				else:
					chrom=snp.split('\t')[header.index('CHROM')]
					pos=snp.split('\t')[header.index('POS')]
					ref=snp.split('\t')[header.index('REF')]
					id=snp.split('\t')[header.index('ID')]
					alt=snp.split('\t')[header.index('ALT')]
					id_var='\t'.join([chrom,pos,ref,alt])
					lista_snp[id_var]=snp


			if not any(x.startswith('CHROM') for x in inters_all):
				inters_all.write('\t'.join(header) + '\n')

			for snp in file2:
				snp=snp.rstrip()
				if snp.startswith('CHROM'):
					header=snp.split('\t')
				else:
					chrom=snp.split('\t')[header.index('CHROM')]
					pos=snp.split('\t')[header.index('POS')]
					ref=snp.split('\t')[header.index('REF')]
					alt=snp.split('\t')[header.index('ALT')]
					id_var='\t'.join([chrom,pos,ref,alt])
					if id_var in lista_snp.keys():
						inters.write(lista_snp[id_var]+'\n')
						inters_all.write(lista_snp[id_var]+'\n')
					else:
						sanger.write('\t'.join([chrom,pos,ref,alt])+'\n')

			file2.close()
			file1.close()
			inters.close()
			inters_all.close()
			sanger.close()
main()
