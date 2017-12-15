import argparse


def main():
	parser = argparse.ArgumentParser('Parse VCF HEADER from FreeBayes,VarScan,GATK to fix it.  Output is to stdout.')
	parser.add_argument('-f', '--infile', help="file di varianti da controllare")
	parser.add_argument('-d','--brca_exc', help="database brca exchange",default=None)
	parser.add_argument('-t','--tipo', help="tipo di pathogenicity da ricercare",default='unknown')
	parser.add_argument('-o', '--out', help="out_file")
	parser.add_argument('--clear', help="flag to clear file, delete double variants",action='store_true')
	parser.add_argument('-a', '--anno', help="anno di analisi del paziente per filtrare i dati")



	global opts
	opts = parser.parse_args()
	#print opts
	if opts.clear:
		infile=open(opts.infile,'r')
		out=open(opts.out+'.txt','w')

		for line in infile:
			var1=[]
			line=line.rstrip()
			if line.startswith('Paziente'):
				head=line.split('\t')
				out.write('\t'.join(head)+'\n')
			else:
				line_split=line.split('\t')
				paz=line.split('\t')[0]
				varianti=line_split[1:]
				for var in varianti:
					if var in var1:
						continue
					else:
						var1+=[var]
				out.write('\t'.join([paz]+var1)+'\n')
	else:
		probandi=dict()
		#probandi_list=open('/home/jarvis/Scrivania/BRCA_STATS/familiari.txt','r')
		probandi_list=open('/home/jarvis/Scrivania/BRCA_STATS/probandi.txt','r')
		for line in probandi_list:
			line=line.rstrip()
			if line.startswith('reg_ricet_ricevuta'):
				head=line.split('\t')
				
			else:
				line_split=line.split('\t')
				nome=line_split[head.index('reg_anag_nome')]
				cognome=line_split[head.index('reg_anag_cognome')]
				data=((line_split[head.index('reg_ricet_dataprenotazione')]).split(' ')[0]).split('/')[2]
				probandi[cognome+' '+nome]=data


		var=open(opts.infile,'r')
		out=open(opts.out+'.txt','w')
		out2=open(opts.out+'_non_trovati.txt','w')
		out_tipo=open(opts.out+'_'+opts.tipo+'.txt','w')
		out_pathogenic=open(opts.out+'_pathogenic.txt','w')
		out_likely_pathogenic=open(opts.out+'_likely_pathogenic.txt','w')
		out_benign=open(opts.out+'_Benign.txt','w')
		out_likely_benign=open(opts.out+'_likely_benign.txt','w')
		
		to_print='\t'.join(['CHR','POS','ID','REF','ALT','GENE','HGVSc','HGVSp','IVS','MAF','PATHOGENICITY'])
		out.write(to_print+'\n')
		out_pathogenic.write(to_print+'\n')
		out_likely_pathogenic.write(to_print+'\n')
		out_tipo.write(to_print+'\n')
		out_benign.write(to_print+'\n')
		out_likely_benign.write(to_print+'\n')
		out2.write('\t'.join(['PAZ','HGVSc','HGVSp'])+'\n')
		num_prob=0
		num_paz=0
		num_patho=0
		num_tipo=0
		num_likely_patho=0
		num_var=0
		num_likely_benign=0
		num_benign=0
		paz_patho=dict()
		paz_likely_patho=dict()
		paz_tipo=dict()
		paz_benign=dict()
		paz_likely_benign=dict()
		
		for line in var:
			line=line.rstrip()
			if line.startswith('Paziente'):
				head=line.split('\t')
			else:
				num_paz+=1
				line_split=line.split('\t')
				anno=line.split('\t')[1]
				paz=line.split('\t')[0]
				probando=line.split('\t')[2]
				varianti=line_split[3:]
				#print paz,paz in probandi.keys()
				if paz in probandi.keys() and opts.anno == probandi[paz]:
					num_prob+=1
					for elem in varianti:
						if elem == '':
							continue
						num_var+=1
						find=0
						HGVSc=elem.split(' ')[0]
						try:
							HGVSp=elem.split(' ')[1]
						except:
							HGVSp='.'
						brca=open(opts.brca_exc,'r')
						for line in brca:
							line=line.rstrip()
							if line.startswith('id'):
								header=line.split('\t')
							else:
								brca_HGVSc=(line.split('\t')[header.index('HGVS_cDNA')]).split(':')[1]
								try:
									brca_HGVSp='p.'+(((line.split('\t')[header.index('HGVS_Protein')]).split(':')[1]).split('(')[1]).split(')')[0]
								except:
									brca_HGVSp=''
								bic_nom=(line.split('\t')[header.index('BIC_Nomenclature')])
								if brca_HGVSc==HGVSc or brca_HGVSp==HGVSp or HGVSc==bic_nom or HGVSp==bic_nom:
									freq=line.split('\t')[header.index('Allele_Frequency')]
									chr=(line.split('\t')[header.index('Genomic_Coordinate_hg37')]).split(':')[0]
									pos=((line.split('\t')[header.index('Genomic_Coordinate_hg37')]).split(':')[1]).split('.')[1]
									ref,alt=((line.split('\t')[header.index('Genomic_Coordinate_hg37')]).split(':')[2]).split('>')
									pathogen=line.split('\t')[header.index('Pathogenicity_all')]
									gene=line.split('\t')[header.index('Gene_Symbol')]
									out.write('\t'.join([chr,pos,paz,ref,alt,gene,HGVSc,HGVSp,bic_nom,freq,pathogen])+'\n')
									if 'Pathogenic' in pathogen:
										out_pathogenic.write('\t'.join([chr,pos,paz,ref,alt,gene,HGVSc,HGVSp,bic_nom,freq,pathogen])+'\n')
										paz_patho[paz]=''
										num_patho+=1
									elif 'Likely_pathogenic'in pathogen:
										out_likely_pathogenic.write('\t'.join([chr,pos,paz,ref,alt,gene,HGVSc,HGVSp,bic_nom,freq,pathogen])+'\n')
										paz_likely_patho[paz]=''
										num_likely_patho+=1
									elif 'Benign' in pathogen:
										out_benign.write('\t'.join([chr,pos,paz,ref,alt,gene,HGVSc,HGVSp,bic_nom,freq,pathogen])+'\n')
										paz_benign[paz]=''
										num_benign+=1
									elif 'Likely_benign'in pathogen:
										out_likely_benign.write('\t'.join([chr,pos,paz,ref,alt,gene,HGVSc,HGVSp,bic_nom,freq,pathogen])+'\n')
										paz_likely_benign[paz]=''
										num_likely_benign+=1
									elif opts.tipo in pathogen:
										out_tipo.write('\t'.join([chr,pos,paz,ref,alt,gene,HGVSc,HGVSp,bic_nom,freq,pathogen])+'\n')
										paz_tipo[paz]=''
										num_tipo+=1


									find=1
									brca.close()
									break
						if find== 0:
							out2.write('\t'.join([paz,HGVSc,HGVSp])+'\tNON TROVATO\n')
							out.write('\t'.join(['.','.',paz,'.','.','.',HGVSc,HGVSp,'.','.','.'])+'\n')
				else:
					out2.write(line+'\tNON IN PROBANDI\n')

		print 'Anno 20'+opts.anno
		print 'Numero Varianti: ',num_var
		print 'Numero Varianti Patologiche: ',num_patho
		print 'Numero Varianti Likely Patologiche: ',num_likely_patho
		print 'Numero Varianti Benign: ',num_benign
		print 'Numero Varianti Likely Benign: ',num_likely_benign
		print 'Numero Varianti '+ opts.tipo + ': ',num_tipo
		print 'Numero Pazienti: ',num_paz
		print 'Numero Probandi: ',num_prob
		print 'Numero Pazienti con varianti patologiche: ',len(paz_patho.keys())
		print 'Numero Pazienti con varianti Likely patologiche: ',len(paz_likely_patho.keys())
		print 'Numero Pazienti con varianti Benign: ',len(paz_benign.keys())
		print 'Numero Pazienti con varianti Likely Benign: ',len(paz_likely_benign.keys())
		print 'Numero Pazienti con varianti '+ opts.tipo + ': ',len(paz_tipo.keys())
					



main()


