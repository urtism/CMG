import argparse

def extract_from_targets(db,targets):
	#new_targ = open(opts.out,'w')
	bd_nostri_geni = open(opts.out,'w')
	esone = 0
	gene_prev = ''
	lista_geni = []
	lista_geni_esoni = []
	for line in targets:
		line = line.rstrip()
		if line.startswith('@'):
			continue
		targ = line.split('\t')[-1]
		lista_geni = lista_geni + [targ.split('.')[0]]
		lista_geni_esoni = lista_geni_esoni+ targ.split('\t')

	for line in db:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		if line.split('\t')[-1] in lista_geni:
			#buffer_line = buffer_line + [line]
			bd_nostri_geni.write(line + '\n')



def extract(geni):
	buff = ''
	trascr_can = open(opts.out,'w')
	for line in geni:
		line = line.rstrip()
		gene = line.split('\t')[-1]
		chrom = line.split('\t')[1]
		start = line.split('\t')[3]
		stop = line.split('\t')[4]
		n_exons = line.split('\t')[5]
		start_exons = line.split('\t')[6]
		stop_exons = line.split('\t')[7]
		if gene != buff.split('\t')[-1]:
			if buff != '':
				trascr_can.write(buff + '\n')
			buff = line
		else:
			buff_l = int(buff.split('\t')[4]) - int(buff.split('\t')[3])
			var_l = int(stop) - int(start)
			if var_l > buff_l:
				buff = line


def print_exons(exons_list):
	new_exons_list = open(opts.out,'w')

	for line in exons_list:
		line = line.rstrip()
		gene = line.split('\t')[-1]
		chrom = line.split('\t')[1]
		strand = line.split('\t')[2]
		start = line.split('\t')[3]
		stop = line.split('\t')[4]
		n_exons = line.split('\t')[5]
		start_exons = line.split('\t')[6]
		stop_exons = line.split('\t')[7]
		i = 0 
		for pos_start in start_exons.split(',')[0:-1]:
			if pos_start !='':
				
				pos_stop = (stop_exons.split(','))[i]
				i = i+1 
				if strand == '+':
					print gene,pos_start,i,strand,"entrato nello strand +"
					new_exons_list.write('\t'.join([gene,chrom,pos_start,pos_stop,str(i)]) + '\n')
				elif strand == '-':
					#print gene,pos_start,i,strand,"entrato nello strand -"
					new_exons_list.write('\t'.join([gene,chrom,pos_start,pos_stop,str(len(start_exons.split(',')[0:-1])-i+1)]) + '\n')





	# for line in targets:
	# 	esone = esone + 1
	# 	line = line.rstrip()
	# 	if line.startswith('@'):
	# 		continue
	# 	targ = line.split('\t')[-1]
	# 	gene = targ.split('.')[0]
	# 	if gene_prev != gene:
	# 		esone = 1
	# 		gene_prev = gene
	# 	chrom = targ.split('.')[1]
	# 	start = targ.split('.')[2]
	# 	stop = targ.split('.')[3]
	# 	new_target = '\t'.join([gene,chrom,start,stop,str(esone)])
	# 	new_targ.write(('\t'.join([gene,chrom,start,stop,str(esone)])) +'\n')
	# new_targ.close()

def check_exon(exons, gene, pos):
	for exon in exons:
		exon = exon.rstrip()
		if exon.split('\t')[0] == gene:
			print gene,pos,exon.split('\t')[2],exon.split('\t')[3], int(pos) >= int(exon.split('\t')[2])
		if exon.split('\t')[0] == gene and int(pos) >= int(exon.split('\t')[2]) and int(pos) <= int(exon.split('\t')[3]):
			num_exon = exon.split('\t')[4]
			print 'ENTRATO'
			return num_exon
	return 'FUORI_TARGET'

def add_exon(varianti, gene, id_paz):
	new_varianti = open(opts.out_path + '/' + id_paz + '_' + opts.tipo + '.txt','w')
	header = 0
	for line in varianti:
		line = line.rstrip()
		if line.startswith('CHROM') and header == 0:
			new_varianti.write('NUM_EXON' + '\t' + line + '\n')
			header = 1
			continue
		elif line.startswith('CHROM') and header != 0:
			continue
		id_paz_var = line.split('\t')[2]
		pos_var = line.split('\t')[1]
		if opts.tipo == 'GATK':
			gene_var = line.split('\t')[25]
		elif opts.tipo == 'Free':
			gene_var = line.split('\t')[28]
		elif opts.tipo == 'Other':
			gene_var = line.split('\t')[25]
		elif opts.tipo == 'VarScan':
			gene_var = line.split('\t')[14]
		if id_paz_var == id_paz and gene_var == gene:
			num_exon = check_exon(open(opts.targets,'r'), gene, pos_var)
			#new_varianti.write(num_exon + '\t' + line + '\n')
			if num_exon != 'FUORI_TARGET':
				new_varianti.write(num_exon + '\t' + line + '\n')
	new_varianti.close()

def main():
	parser = argparse.ArgumentParser('Adds number of exon for each interval')
	
	parser.add_argument('-f','--file',help="variants file")
	parser.add_argument('-t','--targets',help="targets file",default= None)
	parser.add_argument('-o','--out_path',help="new variants file path",default= None)
	parser.add_argument('--out',help="new variants file path",default= None)
	parser.add_argument('-a','--analisi',help="analysis type")
	parser.add_argument('--tipo',help="GATK,Free,Other,VarScan",default= None)
	parser.add_argument('--list',help="lista di pazienti e geni",default= None)

	global opts
	opts = parser.parse_args()
	if opts.analisi == '0':
		#extract_from_targets(open(opts.file,'r'),open(opts.targets,'r'))
		#extract(open(opts.file,'r'))
		print_exons(open(opts.file,'r'))

	elif opts.analisi == '1':
		for line in open(opts.list,'r'):
			line = line.rstrip()
			gene = line.split('\t')[0]
			id_paz = line.split('\t')[1]
			add_exon(open(opts.file,'r'), gene, id_paz)
main()