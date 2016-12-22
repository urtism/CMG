import os
import glob2
import argparse

# def extract_var_from_gene(varianti,vc):
# 	var_list = open(opts.out + '/lista_varianti_' + opts.gene + '_' + vc + '.list','a+')
# 	var_array = var_list.readlines()
# 	index_of_gene = 0
# 	for line in varianti:
# 		line = line.rstrip()
# 		#print line
# 		# if line in var_array:
# 		# 	print 'header'

# 		if line.startswith('CHROM'):
# 			header = line + '\n'
# 			if header not in var_array: 
# 				var_list.write(line + '\n')

# 			index_of_gene = header.split('\t').index('SYMBOL')
# 		else:
# 			variante = line.split('\t')
# 			chrom = variante[0]
# 			pos = variante[1]
# 			id = variante[2]
# 			ref = variante[3]
# 			alt = variante[4]
# 			info = chrom = variante[4:]

# 			var_id = '_'.join([chrom,pos,id,ref,alt])
# 			hash_var_paz[var_id] = [vc,variante]
# 			hash_var[]

# 		for var in hash_var.items():
# 			if variante[index_of_gene] == opts.gene:
# 				var_list.write(line + '\n')
# 	#print var_array

def extract_info(varianti,hash_var,pazienti):
	index_of_gene = 0
	index_of_gen = 0
	buff = []
	for line in varianti:
		line = line.rstrip()
		if line.startswith('CHROM'):
			#print line
			index_of_gene = line.split('\t').index('SYMBOL')
			index_of_gen = line.split('\t').index('GT')
			
			try:
				index_of_HGVSc = line.split('\t').index('HGVSc')
			except:
				print "non trova hgsvc",line
			
			try:
				index_of_HGVSp = line.split('\t').index('HGVSp')
			except:
				print "non trova hgsvp",line
			try:
				index_of_conseq = line.split('\t').index('Consequence')
			except:
				print "non trova CONSEQUENCE",line

			#print index_of_gen
		else:
			variante = line.split('\t')
			chrom = variante[0]
			pos = variante[1]
			id = variante[2]
			ref = variante[3]
			alt = variante[4]
			HGVSc = variante[index_of_HGVSc] 
			HGVSp = variante[index_of_HGVSp]
			CONSEQUENCE = variante[index_of_conseq]
			info = variante[4:]
			var_id = '\t'.join([chrom,pos,ref,alt])
			id_buff = var_id + '\t' + id
			if id_buff not in buff:
				buff = buff + [id_buff]
				
				#print index_of_gen
				if var_id in hash_var.keys():
					if variante[index_of_gen].split('/')[0] == variante[index_of_gen].split('/')[1]:
						hash_var[var_id][0] = hash_var[var_id][0] + 2
						hash_var[var_id][2] = hash_var[var_id][2] + 1
					else: 
						hash_var[var_id][0] = hash_var[var_id][0] + 1
						hash_var[var_id][1] = hash_var[var_id][1] + 1
					try:
						hash_var[var_id][3] = hash_var[var_id][3] + [pazienti[id][0]]
					except:
						hash_var[var_id][3] = hash_var[var_id][3] + [id]
						print id
					hash_var[var_id][4] = HGVSc
					hash_var[var_id][5] = HGVSp
					hash_var[var_id][6] = CONSEQUENCE
				else:
					if variante[index_of_gen].split('/')[0] == variante[index_of_gen].split('/')[1]:
						try:
							hash_var[var_id] = [2,0,1,[pazienti[id][0]],HGVSc,HGVSp,CONSEQUENCE]
						except:
							hash_var[var_id] = [2,0,1,[id],HGVSc,HGVSp,CONSEQUENCE]
							print id	
					else:
						try:
							hash_var[var_id] = [1,1,0,[pazienti[id][0]],HGVSc,HGVSp,CONSEQUENCE]
						except:
							hash_var[var_id] = [2,0,1,[id],HGVSc,HGVSp,CONSEQUENCE]
							print id	
			else:
				continue

def extract_var_from_gene(varianti,vc):
	var_list = open(opts.out + '/lista_varianti_' + opts.gene + '_' + vc + '.list','a+')
	var_array = var_list.readlines()
	index_of_gene = 0
	for line in varianti:
		line = line.rstrip()
		if line.startswith('CHROM'):
			header = line + '\n'
			if header not in var_array: 
				var_list.write(line + '\n')
			index_of_gene = header.split('\t').index('SYMBOL')
		else:
			variante = line.split('\t')
			if variante[index_of_gene] == opts.gene:
				var_list.write(line + '\n')

def extract_paz_name(pazienti):
	id_info={}
	for paz in pazienti:
		paz = paz.rstrip()
		if paz.startswith('NOME_PAZ') or paz == '':
			continue
		else:
			#print paz
			try:
				nome = paz.split('\t')[0]
				num_paz = paz.split('\t')[1]
				num_run = paz.split('\t')[2]
				data_run = paz.split('\t')[3]
				id = paz.split('\t')[4]
				id_info[id]=[nome,num_paz,num_run,data_run]
			except:
				print "il paziente",paz.split('\t')[0],"ha un formato non valido"
	return id_info		

def main():
	parser = argparse.ArgumentParser('Adds number of exon for each interval')
	
	parser.add_argument('-p','--path',help="variants file")
	parser.add_argument('-g','--gene',help="variants file")
	parser.add_argument('-o','--out',help="variants file")
	parser.add_argument('-l','--paz_list',help="lista dei nomi dei pazienti e degli id")
	global opts
	opts = parser.parse_args()
	global pazienti
	pazienti = extract_paz_name(open(opts.paz_list,'r'))
	
	try:
		os.mkdir(opts.out)
	except:
		print opts.out,"directory esistente"

	hash_var = {}
	num_paz = 0

	for filename in glob2.glob(os.path.join(opts.path,'*.tsv')):
		num_paz = num_paz + 1
		#print filename
		if 'FREEBAYES' in filename:
			vc = 'FREEBAYES'
		elif 'GATK_Other' in filename:
			vc = 'GATK_Other'
		elif 'VarScan' in filename:
			vc = 'VarScan'
		elif 'GATK' in filename:
			vc = 'GATK'

		extract_var_from_gene(open(filename,'r'),vc)
	
	#print num_paz

	for filename in glob2.glob(os.path.join(opts.out,'*.list')):
		extract_info(open(filename,'r'),hash_var,pazienti)

	statistiche = open(opts.out + '/lista_varianti_' + opts.gene + '_STATS.txt','w')
	statistiche.write('CHROM'+'\t' + 'POS'+'\t' + 'REF'+'\t' + 'ALT' +'\t' +'HGVSc'+'\t' + 'HGVSp'+'\t'+'CONSEQUENCE'+'\t'+ 'FRAZ_ALLELICA'+'\t' + 'FRAZ_PERC'+'\t' + 'MAF'+'\t'+'ETERO'+'\t'+ 'OMO'+'\t'+'PAZIENTI'+'\n')
	
	for var in hash_var.keys():
		maf = float("{0:.4f}".format(float(hash_var[var][0])/float(num_paz/2)))
		fraz = str(hash_var[var][0]) + '/' + str(num_paz/2)
		perc = str(maf*100) + '%'
		info_var = str(hash_var[var][4]) + '\t' + str(hash_var[var][5])+ '\t' + str(hash_var[var][6])
		stats = fraz + '\t' + perc + '\t' +  str(maf) + '\t' + str(hash_var[var][1]) + '\t' + str(hash_var[var][2]) + '\t' + ';'.join(hash_var[var][3])
		statistiche.write(var + '\t' +info_var + '\t' + stats + '\n')

main()
