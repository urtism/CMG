import argparse

def extract_var_data(value,line, file, HEAD, i):

	new_line = []

	new_line += [line[HEAD.index('#CHROM')]]
	new_line += [line[HEAD.index('POS')]]
	new_line += [HEAD[i]]
	new_line += [line[HEAD.index('REF')]]
	new_line += [line[HEAD.index('ALT')]]
	new_line += [str(value)]

	if 'FreeBayes' in file:
		new_line += ['FreeBayes']
	elif 'VarScan' in file:
		new_line += ['VarScan']
	elif 'GATK' in file:
		new_line += ['GATK']
	else:
		new_line += ['.']

	new_line += [line[i].split(':')[line[HEAD.index('FORMAT')].split(':').index('GT')]]
	new_line += [line[i].split(':')[line[HEAD.index('FORMAT')].split(':').index('DP')]]

	#try:
	if 'FreeBayes' in file:
		AD = float(line[i].split(':')[line[HEAD.index('FORMAT')].split(':').index('AO')])
		RD = float(line[i].split(':')[line[HEAD.index('FORMAT')].split(':').index('RO')])
		try:
			FREQ = round(AD/(AD+RD),4)
		except:
			FREQ = '.'

		new_line += [str(FREQ)]

	elif 'VarScan' in file:
		RD = float(line[i].split(':')[line[HEAD.index('FORMAT')].split(':').index('RD')])
		AD = float(line[i].split(':')[line[HEAD.index('FORMAT')].split(':').index('AD')])
		try:
			FREQ = round(AD/(AD+RD),4)
		except:
			FREQ = '.'

		new_line += [str(FREQ)]

	elif 'GATK' in file:
		RD = float(line[i].split(':')[line[HEAD.index('FORMAT')].split(':').index('AD')].split(',')[0])
		AD = float(line[i].split(':')[line[HEAD.index('FORMAT')].split(':').index('AD')].split(',')[1])
		try:
			FREQ = round(AD/(AD+RD),4)
		except:
			FREQ = '.'

		new_line += [str(FREQ)]

	return new_line















def main():

	parser = argparse.ArgumentParser('\n\nQuesto tool estrae le varianti dai vcf che presentano caratteristiche di iEVA che ci interessano\n')
	parser.add_argument('-I','--input',help="Input file in formato .list contenente il path dei vcf da analizzare e da cui estrarre le varianti. NOTA BENE: non usare il tool su diversi variant caller!!")
	parser.add_argument('-O','--outfile',help="file di output in vcf format. Indicare come path/to/nome_file --> In uscita si avra un file per ciascun variant caller (e.g. nome_file_GATK.vcf) usato (prendo il nome del variant caller dal nome del file vcf).")

	global opts
	opts = parser.parse_args()
	out = opts.outfile
	Variant_list = {}
	Variant_table = {}

	with open(opts.input,'r') as tsv:
		
		for file in tsv:
			file = file.rstrip()
			
			with open(file,'r') as vcf:
				
				for line in vcf:
					line = line.rstrip()

					if line.startswith('#'):
						if line.startswith('#CHROM'):
							Header_vcf = line.split('\t')
							Header_out1 = ['CHROM','POS','ID','REF','ALT','MMQ','CALLER','GT','DP','AD']
							Header_out2 = ['CHROM','POS','ID','REF','ALT','MQ0','CALLER','GT','DP','AD']
							Header_out3 = ['CHROM','POS','ID','REF','ALT','NPA','CALLER','GT','DP','AD']
							Header_out4 = ['CHROM','POS','ID','REF','ALT','SA','CALLER','GT','DP','AD']
							Header_out5 = ['CHROM','POS','ID','REF','ALT','NP','CALLER','GT','DP','AD']
							Header_out6 = ['CHROM','POS','ID','REF','ALT','NPP','CALLER','GT','DP','AD']

							out1 = out + 'MMQ<35.tsv'
							outfile1 = open(out1, 'w')
							outfile1.write('\t'.join(Header_out1) + '\n')
							outfile1.close()

						else:
							continue

					else:
						line = line.split('\t')

						variant_ID = str(line[Header_vcf.index('#CHROM')])+'-'+str(line[Header_vcf.index('POS')])+'-'+str(line[Header_vcf.index('REF')])+'-'+str(line[Header_vcf.index('ALT')])
						
						for i in range(9,(len(line)-1)):

							try:
								Extract_MMQ = float(line[i].split(':')[line[8].split(':').index('MMQ')])
							except:
								pass

							try:
								Extract_MQ0 = float(line[i].split(':')[line[8].split(':').index('MQ0')])
							except:
								pass

							try:
								Extract_NPA = float(line[i].split(':')[line[8].split(':').index('NPA')])
							except:
								pass

							try:
								Extract_SA = float(line[i].split(':')[line[8].split(':').index('SA')])
							except:
								pass

							try:
								Extract_NP = float(line[i].split(':')[line[8].split(':').index('NP')])
							except:
								pass

							try:
								Extract_NPP = float(line[i].split(':')[line[8].split(':').index('NPP')])
							except:
								pass

							GT = line[i].split(':')[line[Header_vcf.index('FORMAT')].split(':').index('GT')]

							if Extract_MMQ < 35 and (GT == '0/1' or GT == '1/1'):

								new_line = extract_var_data(Extract_MMQ, line, file, Header_vcf, i)

								outfile1 = open(out1, 'a')
								outfile1.write('\t'.join(new_line) + '\n')
								outfile1.close()

							new_line = []

							#if line[i].split(':') is '0':
							#	print 'ok'

						#if variant_ID in Variant_list:
						#	print 'OLD'
						#	continue
						#else:
						#	Variant_list[variant_ID] = 0
						#	print 'new_var'						
						





main()