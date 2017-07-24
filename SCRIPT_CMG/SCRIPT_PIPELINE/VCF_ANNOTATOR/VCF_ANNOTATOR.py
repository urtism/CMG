import argparse
import re
import string




#------------------------------------------------------ Funzione per aggiungere INFO -----------------------------------------------------




def extract_header(new_vep):
	start_vep = new_vep.find('Allele')
	end_vep = new_vep.find('">')
	header_ann = new_vep[start_vep:end_vep]
	header_ann = header_ann.split('|')
	return header_ann



def add_INFO(vep):

	#Modifico l'header che inizia con ##INFO contenente l'annotazione di VEP aggiungendo in coda le tag inserite nella command line
	vep = vep.replace('">', '|">')

	if opts.clinvar:
		tag_cl = (opts.clinvar).split(',')[1:]
		new_tag_cl = '|'.join(tag_cl)
		vep = vep.replace('">', str(new_tag_cl)+'|">')

	if opts.Humsavar:
		tag_hum = (opts.Humsavar).split(',')[1:]
		new_tag_hum = '|'.join(tag_hum)
		vep = vep.replace('">', str(new_tag_hum)+'|">')

	if opts.exome_variant_server:
		tag_esp = (opts.exome_variant_server).split(',')[1:]
		new_tag_esp = '|'.join(tag_esp)
		vep = vep.replace('">', str(new_tag_esp)+'|">')

	if opts.gerp:
		tag_gerp = (opts.gerp).split(',')[1:]
		new_tag_gerp = '|'.join(tag_gerp)
		vep = vep.replace('">', str(new_tag_gerp)+'|">')

	if opts.phastCons:
		tag_ph = (opts.phastCons).split(',')[1:]
		new_ph = ['phastCons_' + spec for spec in tag_ph]
		new_tag_ph = '|'.join(new_ph)
		vep = vep.replace('">', str(new_tag_ph)+'|">')

	if opts.phyloP:
		tag_phy = (opts.phyloP).split(',')[1:]
		new_phy = ['phyloP_' + spe for spe in tag_phy]
		new_tag_phy = '|'.join(new_phy)
		vep = vep.replace('">', str(new_tag_phy)+'|">')

	vep = vep.replace('|">', '">')

	return vep




#------------------------------------------------------ Funzione per aggiungere ClinVar -----------------------------------------------------




def Add_ClinVar(line):

	#Estraggo le tage dopo il path del db: es. path/to/ClinVar,CLNSIG,CLNDB --> prendo solo CLNSIG e CLNDB e li metto in tag
	#iIn clnv invece leggo le linee di clinvar
	clnv = open((opts.clinvar).split(',')[0],'r')
	tag = (opts.clinvar).split(',')[1:]
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
			val_add=[]

			for elem in tag:
				val = '.'
				for stat in INFO:
					stat = stat.replace('|', '/')
					if elem in stat:
						stat = stat.split('=')[-1]
							#print val

							#Sostituisco ai numeri i valori del dict

						if 'CLNSIG' in elem:
							stat = stat.replace('255', 'Other')
							stat = stat.replace('0', 'Uncertain significance')
							stat = stat.replace('1', 'Not provided')
							stat = stat.replace('3', 'Likely benign')
							stat = stat.replace('4', 'Likely pathogenic')
							stat = stat.replace('6', 'Drug response')
							stat = stat.replace('7', 'histocompatibility')
							stat = stat.replace('2', 'Benign')
							stat = stat.replace('5', 'Pathogenic')
							val = stat
							continue

						if 'CLNORIGIN' in elem:
							stat = stat.replace('1073741824', 'Other')
							stat = stat.replace('512', 'Tested-inconclusive')
							stat = stat.replace('256', 'Not-tested')
							stat = stat.replace('128', 'Uniparental')
							stat = stat.replace('64', 'Biparental')
							stat = stat.replace('32', 'De Novo')
							stat = stat.replace('16', 'Maternal')
							stat = stat.replace('8', 'Paternal')
							stat = stat.replace('4', 'Inherited')
							stat = stat.replace('3', 'Germline')
							stat = stat.replace('2', 'Somatic')
							stat = stat.replace('1', 'Germline')
							stat = stat.replace('0', 'Unknown')
							val = stat
							continue

						if 'SAO' in elem:
							stat = stat.replace('0', 'Unspecified')
							stat = stat.replace('1', 'Germline')
							stat = stat.replace('2', 'Somatic')
							stat = stat.replace('3', 'Both')
							val = stat
							continue

						if 'SSR' in elem:
							stat = stat.replace('1024', 'Other')
							stat = stat.replace('16', '1kg_failed')
							stat = stat.replace('8', 'Para_EST')
							stat = stat.replace('4', 'oldAlign')
							stat = stat.replace('2', 'byEST')
							stat = stat.replace('1', 'Paralog')
							stat = stat.replace('0', 'Unspecified')
							val = stat
							continue

						else:
							val = stat

				val_add.append(val)
			#Una volta trovata la variante e processata, allora esco dal for e passo alla linea dopo
			break
	line[7] = str(line[7] + '|' +'|'.join(val_add))
	return line
	clnv.close()	
		



#------------------------------------------------------ Funzione per aggiungere HumSavar -----------------------------------------------------




def Add_Humsavar(line,header_ann):

	humsavar = open((opts.Humsavar).split(',')[0],'r')
	tag = (opts.Humsavar).split(',')[1:]

	val_add=[]

	for elem in tag:
		val_add = val_add + ['.']

	tab_INFO = line[7].split('|')

	for raw in humsavar:
		
		raw = raw.rstrip()
		raw = raw.split('\t')

		if raw[0].startswith('Main_gene_name'):
			header_humsavar = raw
			continue

		#Controllo che abbiano stesso gene e lo stesso p. :
		if raw[header_humsavar.index('Main_gene_name')] == tab_INFO[header_ann.index('SYMBOL')] and raw[header_humsavar.index('AA_Change')] == tab_INFO[header_ann.index('HGVSp')].split(':')[-1]:
			val_add=[]

			for elem in tag:
				val_add = val_add + [raw[header_humsavar.index(str(elem))]]

			break
		
	line[7] = str(line[7] + '|' +'|'.join(val_add))
	return line
	humsavar.close()



#------------------------------------------------------ Funzione per aggiungere ESP -----------------------------------------------------




def Add_ESP(line):

	#Estraggo le tag dopo il path del db: es. path/to/ESP,ESP_All,ESP_EA --> prendo si valori e li metto in tag
	#In ESP invece leggo le linee di ESP
	ESP = open((opts.exome_variant_server).split(',')[0],'r')
	tag = (opts.exome_variant_server).split(',')[1:]

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

			#print new_tag
				
			val_add=[]

			for elem in tag:
				val = '.'
				for stat in INFO:
					if elem in stat:
						val = stat

					else:
						continue

				val_add.append(val)

			break
	line[7] = str(line[7] + '|' + '|'.join(val_add))
	return line	
	ESP.close()




#------------------------------------------------------ Funzione per aggiungere GERP++ -----------------------------------------------------




def Add_Gerp(line,pan):

	tag = (opts.gerp).split(',')[1:]
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
		if raw[header_gerp.index('CHROM')] == line[0] and raw[header_gerp.index('POS')] == line[1]:
				
			val_add=[]

			for elem in tag:
				val_add = val_add + [raw[header_gerp.index(str(elem))]]

			break

	line[7] = str(line[7] + '|' + '|'.join(val_add))
	return line
	gerp_code.close()




#------------------------------------------------------ Funzione per aggiungere phastCons -------------------------------------------------




def Add_phastCons(line,pan):

	species = (opts.phastCons).split(',')[1:]
	new_spec = ['phastCons_' + spec for spec in species]
	
	val_add=[]

	for val in new_spec:
		val_add = val_add + ['.']

	chr = line[0]
	phastCons = open((opts.phastCons).split(',')[0] + '/' + 'phastCons_' + chr + '_' + pan + '.tsv','r')

	for raw in phastCons:	
		raw = raw.rstrip()
		raw = raw.split('\t')

		if raw[0].startswith('CHROM'):
			header_phast = raw
			continue

		if raw[header_phast.index('CHROM')] == line[0] and raw[header_phast.index('POS')] == line[1]:
			val_add=[]
			for spe in new_spec:
				val_add = val_add + [raw[header_phast.index(str(spe))]]
			break

		else:
			continue
		
	phastCons.close()

	line[7] = str(line[7] + '|' + '|'.join(val_add))
	return line




#------------------------------------------------------ Funzione per aggiungere phyloP ------------------------------------------------------




def Add_phyloP(line,pan):

	species = (opts.phyloP).split(',')[1:]
	new_spec = ['phyloP_' + spec for spec in species]

	val_add=[]

	for val in new_spec:
		val_add = val_add + ['.']

	chr = line[0]
	phyloP = open((opts.phyloP).split(',')[0] + '/' + 'phyloP_' + chr + '_' + pan + '.tsv','r')

	for raw in phyloP:	
		raw = raw.rstrip()
		raw = raw.split('\t')

		if raw[0].startswith('CHROM'):
			header_phyloP = raw
			continue

		if raw[header_phyloP.index('CHROM')] == line[0] and raw[header_phyloP.index('POS')] == line[1]:
			val_add=[]
			for spe in new_spec:
				val_add = val_add + [raw[header_phyloP.index(str(spe))]]
			break
		
		else:
			continue

	phyloP.close()

	line[7] = str(line[7] + '|' + '|'.join(val_add))
	return line





def main():
	
	parser = argparse.ArgumentParser('aggiunge le annotazioni fornite nel file -f (una per riga) al file di input -i in un formato tab delimited in output -o')
	parser.add_argument('-i','--input',help="file vcf annotato con VEP")
	parser.add_argument('-clnv','--clinvar',default=None,help="path del database clinvar in formato vcf. Utilizzo: path/Clinvar.vcf,CLNSIG,CLNDB con CLNDB e CLNSIG aggiungo i valori rispettivi nel vcf ClinVar")
	parser.add_argument('-ESP','--exome_variant_server',default=None,help="path database ESP. Utilizzo: path/ESP.vcf,ESP_EU")
	parser.add_argument('-hvar','--Humsavar',default=None,help="path database Humsavar. Utilizzo: path/Humsavar.txt,Variant_type")
	parser.add_argument('-gerp','--gerp',default=None,help="path database Gerp++. Utilizzo: path_to_GERP_folder,RS_Score")
	parser.add_argument('-phastCons','--phastCons',default=None,help="path database phastCons. Utilizzo: path_to_folder,species,species ---> es: ~/phastCons/primate,vertebrate,")
	parser.add_argument('-phyloP','--phyloP',default=None,help="path database phyloP46way_placental. Utilizzo: path_to_folder,species ---> es: ~/phyloP/placental,vertebrate,")
	parser.add_argument('-P','--pannello',help="Specificare il pannello utilizzato nell'analisi: Cancer | Cardio | Exome | BRCA")
	parser.add_argument('-o','--outfile',help="file di output in vcf format")

	global opts
	
	opts = parser.parse_args()
	vcf = open(opts.input,'r')
	out = open(opts.outfile,'w')
	pan = opts.pannello
	save_head = []

	
	for line in vcf:

		#Vado a lavorare sull'header per estrarre le informazione e modificare il campo INFO
		if line.startswith('#'):
			line = line.rstrip()

			if line.startswith('##INFO=<ID=ANN') or line.startswith('##INFO=<ID=CSQ'):
				new_vep = add_INFO(line)
				save_head.append(new_vep)
				header_ann = extract_header(new_vep)

			elif line.startswith('#CHROM'):
				save_head.append(line)
				out.write('\n'.join(save_head) + '\n')

			else:
				save_head.append(line)

		else:
		
			line = line.rstrip()
			line = line.split('\t')

			if opts.clinvar:
				line = Add_ClinVar(line)

			if opts.Humsavar:
				line = Add_Humsavar(line,header_ann)

			if opts.exome_variant_server:
				line = Add_ESP(line)

			if opts.gerp:
				line = Add_Gerp(line,pan)

			if opts.phastCons:
				line = Add_phastCons(line,pan)

			if opts.phyloP:
				line = Add_phyloP(line,pan)

			print line
			out.write('\t'.join(line) + '\n')

	vcf.close()
	out.close()

main()