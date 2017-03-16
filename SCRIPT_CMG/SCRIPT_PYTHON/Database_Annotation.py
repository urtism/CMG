import argparse
import re
import string




#------------------------------------------------------ Funzione per aggiungere ClinVar -----------------------------------------------------




def Add_ClinVar(line):

	#Estraggo le tage dopo il path del db: es. path/to/ClinVar,CLNSIG,CLNDB --> prendo solo CLNSIG e CLNDB e li metto in tag
	#iIn clnv invece leggo le linee di clinvar
	clnv = open((opts.clinvar).split(',')[0],'r')
	tag = (opts.clinvar).split(',')[1:]

	if line[0] == 'CHROM':
		#Aggiungo all'header le tag che mi servono
		head = line + tag
		return head

	else:
		val_add=[]

		for elem in tag:
			val_add = val_add + ['.']

		for raw in clnv:
			
			raw = raw.rstrip()
			raw = raw.split('\t')

			if raw[0].startswith('##'):
				continue

			if raw[0].startswith('#CHROM'):
				header_clnv = raw
				continue

			#Controllo che abbiano stesso chrom pos ref e alt
			if raw[0] == line[0] and line[1] == raw[1] and line[3] == raw[3] and line[4] == raw[4]:
				INFO=raw[header_clnv.index('INFO')].split(';')
				#print INFO
				val_add=[]

				for elem in tag:
					val = '.'
					for stat in INFO:
						if elem in stat:
							val = stat.split('=')[-1]
							#print elem
							#print val

							#Sostituisco ai numeri i valori del dizionario

							if elem == 'CLNSIG':
								val = string.replace(val, '255', 'Other')
								val = string.replace(val, '0', 'Uncertain significance')
								val = string.replace(val, '1', 'Not provided')
								val = string.replace(val, '3', 'Likely benign')
								val = string.replace(val, '4', 'Likely pathogenic')
								val = string.replace(val, '6', 'Drug response')
								val = string.replace(val, '7', 'histocompatibility')
								val = string.replace(val, '2', 'Benign')
								val = string.replace(val, '5', 'Pathogenic')
								continue

							if elem == 'CLNORIGIN':
								val = string.replace(val, '1073741824', 'Other')
								val = string.replace(val, '512', 'Tested-inconclusive')
								val = string.replace(val, '256', 'Not-tested')
								val = string.replace(val, '128', 'Uniparental')
								val = string.replace(val, '64', 'Biparental')
								val = string.replace(val, '32', 'De Novo')
								val = string.replace(val, '16', 'Maternal')
								val = string.replace(val, '8', 'Paternal')
								val = string.replace(val, '4', 'Inherited')
								val = string.replace(val, '2', 'Somatic')
								val = string.replace(val, '1', 'Germline')
								val = string.replace(val, '0', 'Unknown')
								continue

							if elem == 'SAO':
								val = string.replace(val, '0', 'Unspecified')
								val = string.replace(val, '1', 'Germline')
								val = string.replace(val, '2', 'Somatic')
								val = string.replace(val, '3', 'Both')
								continue

							if elem == 'SSR':
								val = string.replace(val, '1024', 'Other')
								val = string.replace(val, '16', '1kg_failed')
								val = string.replace(val, '8', 'Para_EST')
								val = string.replace(val, '4', 'oldAlign')
								val = string.replace(val, '2', 'byEST')
								val = string.replace(val, '1', 'Paralog')
								val = string.replace(val, '0', 'Unspecified')
								continue

							else:
								continue

					val_add.append(val)
				#Una volta trovata la variante e processata, allora esco dal for e passo alla linea dopo
				break
					
		new_line = line + val_add
		return new_line
	clnv.close()	
		



#------------------------------------------------------ Funzione per aggiungere HumSavar -----------------------------------------------------




def Add_Humsavar(line,header_tsv):

	humsavar = open((opts.Humsavar).split(',')[0],'r')
	tag = (opts.Humsavar).split(',')[1:]

	if line[0] == 'CHROM':
		#Aggiungo all'header le tag che mi servono
		head = line + tag
		return head

	else:

		val_add=[]

		for elem in tag:
			val_add = val_add + ['.']

		for raw in humsavar:
			
			raw = raw.rstrip()
			raw = raw.split('\t')

			if raw[0].startswith('Main_gene_name'):
				header_humsavar = raw
				continue

			#Controllo che abbiano stesso gene e lo stesso p. :
			if raw[header_humsavar.index('Main_gene_name')] == line[header_tsv.index('SYMBOL')] and raw[header_humsavar.index('AA_Change')] == line[header_tsv.index('HGVSp')].split(':')[-1]:
				
				val_add=[]

				for elem in tag:
					val_add = val_add + [raw[header_humsavar.index(str(elem))]]

					break
		
		line = line + val_add
		return line
	humsavar.close()



#------------------------------------------------------ Funzione per aggiungere ESP -----------------------------------------------------




def Add_ESP(line,header_tsv):

	#Estraggo le tage dopo il path del db: es. path/to/ClinVar,CLNSIG,CLNDB --> prendo solo CLNSIG e CLNDB e li metto in tag
	#iIn clnv invece leggo le linee di clinvar
	ESP = open((opts.exome_variant_server).split(',')[0],'r')
	tag = (opts.exome_variant_server).split(',')[1:]

	if line[0] == 'CHROM':
		#Aggiungo all'header le tag che mi servono
		head = line + tag
		return head


	else:

		val_add=[]

		for elem in tag:
			val_add = val_add + ['.']

		for raw in ESP:
			
			raw = raw.rstrip()
			raw = raw.split('\t')

			if raw[0].startswith('##'):
				continue

			if raw[0].startswith('#CHROM'):
				header_ESP = raw
				continue

			#Controllo che abbiano stesso chrom pos ref e alt
			if raw[0] == line[0] and line[1] == raw[1] and line[3] == raw[3] and line[4] == raw[4]:
				INFO=raw[header_ESP.index('INFO')].split(';')
				
				for x in INFO:
					if 'MAF=' in x:
						indice = INFO.index(x)
						new_tag = x.split(',')
						new_tag[0] = string.replace(new_tag[0],'MAF=','')
						new_tag[0] = '{0:f}'.format(float(new_tag[0])/100, 5)
						new_tag[0] = 'ESP_EA=' + str(new_tag[0])
						new_tag[1] = '{0:f}'.format(float(new_tag[1])/100, 5)
						new_tag[1] = 'ESP_AA=' + str(new_tag[1])
						new_tag[2] = '{0:f}'.format(float(new_tag[2])/100, 5)
						new_tag[2] = 'ESP_All=' + str(new_tag[2])
				
				INFO[indice:indice+1] = new_tag
				
				val_add=[]

				for elem in tag:
					val = '.'
					for stat in INFO:
						if elem in stat:
							val = stat.split('=')[-1]

						else:
							continue

					val_add.append(val)

				#Una volta trovata la variante e processata, allora esco dal for e passo alla linea successiva
				break
					
		new_line = line + val_add
		return new_line	
	ESP.close()




#------------------------------------------------------ Funzione per aggiungere GERP++ -----------------------------------------------------




def Add_Gerp(line,header_tsv,pan):

	tag = (opts.gerp).split(',')[1:]

	if line[0] == 'CHROM':
		#Aggiungo all'header le tag che mi servono
		head = line + tag
		return head

	else:
		chr = line[0]
		gerp_code = open((opts.gerp).split(',')[0]+'/gerp.'+chr+'.'+pan+'.rates','r')
		val_add=[]

		for elem in tag:
			val_add = val_add + ['.']

		for raw in gerp_code:
			
			raw = raw.rstrip()
			raw = raw.split('\t')

			if raw[0].startswith('CHROM'):
				header_gerp = raw
				continue

			#Controllo che abbiano stesso gene e lo stesso p. :
			if raw[header_gerp.index('CHROM')] == line[header_tsv.index('CHROM')] and raw[header_gerp.index('POS')] == line[header_tsv.index('POS')]:
				
				val_add=[]

				for elem in tag:
					val_add = val_add + [raw[header_gerp.index(str(elem))]]

					break
		
		line = line + val_add
		return line
		gerp_code.close()




#------------------------------------------------------ Funzione per aggiungere phastCons -------------------------------------------------




def Add_phastCons(line,header_tsv,pan):

	species = (opts.phastCons).split(',')[1:]
	new_spec = ['phastCons_' + spec for spec in species]
	

	if line[0] == 'CHROM':
		head = line + new_spec
		return head

	else:
		val_add=[]

		for val in new_spec:
			val_add = val_add + ['.']

		for spec in new_spec:
			chr = line[0]
			phastCons = open((opts.phastCons).split(',')[0]+spec+'_'+chr+'_'+pan+'.chr.tsv','r')

			for raw in phastCons:	
				raw = raw.rstrip()
				raw = raw.split('\t')

				if raw[0].startswith('CHROM'):
					header_phast = raw
					continue

				#Controllo che abbiano stesso gene e lo stesso p. :
				if raw[header_phast.index('CHROM')] == line[header_tsv.index('CHROM')] and raw[header_phast.index('POS')] == line[header_tsv.index('POS')]:
					val_add[new_spec.index(spec)] = raw[header_phast.index(str(spec))]
					break
			
			phastCons.close()

		line = line + val_add
		return line




#------------------------------------------------------ Funzione per aggiungere phyloP ------------------------------------------------------




def Add_phyloP(line,header_tsv,pan):

	species = (opts.phyloP).split(',')[1:]
	new_spec = ['phyloP_' + spec for spec in species]
	

	if line[0] == 'CHROM':
		head = line + new_spec
		return head

	else:
		val_add=[]

		for val in new_spec:
			val_add = val_add + ['.']

		for spec in new_spec:
			chr = line[0]
			phyloP = open((opts.phyloP).split(',')[0]+spec+'_'+chr+'_'+pan+'.chr.tsv','r')

			for raw in phyloP:	
				raw = raw.rstrip()
				raw = raw.split('\t')

				if raw[0].startswith('CHROM'):
					header_phyloP = raw
					continue

				#Controllo che abbiano stesso gene e lo stesso p. :
				if raw[header_phyloP.index('CHROM')] == line[header_tsv.index('CHROM')] and raw[header_phyloP.index('POS')] == line[header_tsv.index('POS')]:
					val_add[new_spec.index(spec)] = raw[header_phyloP.index(str(spec))]
					break
			
			phyloP.close()

		line = line + val_add
		return line















def main():
	parser = argparse.ArgumentParser('aggiunge le annotazioni fornite nel file -f (una per riga) al file di input -i in un formato tab delimited in output -o')
	
	parser.add_argument('-i','--input',help="file tab delimited annotato")
	parser.add_argument('-clnv','--clinvar',default=None,help="path del database clinvar in formato vcf. Utilizzo: path/Clinvar.vcf,CLNSIG,CLNDB con CLNDB e CLNSIG aggiungo i valori rispettivi nel vcf ClinVar")
	parser.add_argument('-ESP','--exome_variant_server',default=None,help="path database ESP. Utilizzo: path/ESP.vcf,ESP_EU")
	parser.add_argument('-hvar','--Humsavar',default=None,help="path database Humsavar. Utilizzo: path/Humsavar.txt,Variant_type")
	parser.add_argument('-gerp','--gerp',default=None,help="path database Gerp++. Utilizzo: path_to_GERP_folder,RS_Score")
	parser.add_argument('-phastCons','--phastCons',default=None,help="path database phastCons. Utilizzo: path_to_folder,species,species ---> es: ~/phastCons/,primate,vertebrate,")
	parser.add_argument('-phyloP','--phyloP',default=None,help="path database phyloP46way_placental. Utilizzo: path_to_folder,species ---> es: ~/phyloP/placental,vertebrate,")
	parser.add_argument('-COSMIC','--cosmic',default=None,help="path database COSMIC es: ~/COSMIC/valore_da_aggiungere,")	
	parser.add_argument('-P','--pannello',help="Specificare il pannello utilizzato nell'analisi: Cancer | Cardio | Exome | BRCA")
	parser.add_argument('-o','--outfile',help="file di output tab delimited")

	global opts
	
	opts = parser.parse_args()
	tab = open(opts.input,'r')
	out = open(opts.outfile,'w')
	pan = opts.pannello


	
	for line in tab:

		line = line.rstrip()
		line = line.split('\t')

		if line[0] == 'CHROM':
			header_tsv = line
		
		if opts.clinvar:
			line = Add_ClinVar(line)

		if opts.Humsavar:
			line = Add_Humsavar(line,header_tsv)

		if opts.exome_variant_server:
			line = Add_ESP(line,header_tsv)

		if opts.gerp:
			line = Add_Gerp(line,header_tsv,pan)

		if opts.phastCons:
			line = Add_phastCons(line,header_tsv,pan)

		if opts.phyloP:
			line = Add_phyloP(line,header_tsv,pan)

		if opts.cosmic:
			line = Add_cosmic(line,header_tsv,pan)

		print line

		out.write('\t'.join(line) + '\n')

main()