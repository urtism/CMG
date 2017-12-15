import argparse



def main():
	parser = argparse.ArgumentParser('Intersect files from conn and cardio DB')
	
	parser.add_argument('--conn',help="path del file di varianti dei connettivi")
	parser.add_argument('--cardio',help="path del file di varianti dei cardio")
	parser.add_argument('-o','--out',help="path di output")

	global opts
	opts = parser.parse_args()
	
	cardio=dict()
	conn=dict()
	condivisi=dict()

	file_conn=open(opts.conn,'r')
	file_cardio=open(opts.cardio,'r')
	
	out_cardio=open(opts.out+'.Cardio.tsv','w')
	out_conn=open(opts.out+'.Conn.tsv','w')
	out_condivisi=open(opts.out+'.Cardio+Conn.tsv','w')

	
	for var in file_cardio:
		var=var.rstrip()
		if var.startswith('CHROM'):
			header=var.split('\t')
		else:
			chrom=var.split('\t')[header.index('CHROM')]
			pos=var.split('\t')[header.index('POS')]
			ref=var.split('\t')[header.index('REF')]
			alt=var.split('\t')[header.index('ALT')]
			id_var='\t'.join([chrom,pos,ref,alt])
			cardio[id_var]=var

	for var in file_conn:
		var=var.rstrip()
		if var.startswith('CHROM'):
			header=var.split('\t')
		else:
			chrom=var.split('\t')[header.index('CHROM')]
			pos=var.split('\t')[header.index('POS')]
			ref=var.split('\t')[header.index('REF')]
			alt=var.split('\t')[header.index('ALT')]
			id_var='\t'.join([chrom,pos,ref,alt])
			conn[id_var]=var		

	out_cardio.write('\t'.join(header)+'\n')
	out_conn.write('\t'.join(header)+'\n')

	header_cond=["CHROM","POS","REF","ALT","HGVSc","HGVSp","CONSEQUENCE","FRAZ_ALLELICA_CONN","FRAZ_PERC_CONN","MAF_CONN","ETERO_CONN","OMO_CONN","NUM_PAZ_MUTATI_CONN","PAZIENTI_HET_CONN","PAZIENTI_HOM_CONN","FRAZ_ALLELICA_CARDIO","FRAZ_PERC_CARDIO","MAF_CARDIO","ETERO_CARDIO","OMO_CARDIO","NUM_PAZ_MUTATI_CARDIO","PAZIENTI_HET_CARDIO","PAZIENTI_HOM_CARDIO"] 
	out_condivisi.write('\t'.join(header_cond)+'\n')



	for id_var in cardio.keys():
		if id_var in conn.keys():
			condivisi[id_var]='\t'.join(conn[id_var].split('\t')+cardio[id_var].split('\t')[6:])
			del cardio[id_var]
			del conn[id_var]
		else:
			out_cardio.write(cardio[id_var]+'\n')

	for id_var in conn.keys():
		if id_var in cardio.keys():
			condivisi[id_var]='\t'.join(conn[id_var].split('\t')+cardio[id_var].split('\t')[6:])
			del cardio[id_var]
			del conn[id_var]
		else:
			out_conn.write(conn[id_var]+'\n')

	for id_var in condivisi.keys():
		out_condivisi.write(condivisi[id_var]+'\n')


main()
