import argparse
import re
import string
import sys


def extract_feat_from_vcf(file,var_list,features,sample):
	v=[]
	vars=[]
	
	for var_id,CLASS in var_list:
		var_id.split('\t')
		var = [sample,var_id.split('\t')[0]+':'+var_id.split('\t')[1]+'-'+var_id.split('\t')[3]+'-'+var_id.split('\t')[4],'NON TROVATA']		
		with open(file,'r') as vcf:
			for line in vcf:
				line = line.rstrip()
				if line.startswith('#CHROM'):
					header = line.rstrip().split('\t')
				if line.startswith(var_id) and sample in header:
					chr,pos,id,ref,alt,filter,qual,info,format = line.split('\t')[:9]
					var = [sample,chr+':'+pos+'-'+ref+'-'+alt]
					try:
						format_sample = line.split('\t')[header.index(sample)].split(':')
					except:
						new_sample = string.replace(sample, 'Cardio', 'Conn')
						format_sample = line.split('\t')[header.index(new_sample)].split(':')
					format_split = format.split(':')
					info_split = info.split(';')
					
					for f in features[:-1]:
						newv = ['?']
						htag,tag = f.split('-')
						if htag == 'INFO' or htag == 'info':
							if VariantCaller == 'SNVer':	
								if tag == 'AF':
									dp = ''
									ac = ''
									for itag in info_split:
										if itag.startswith('DP='):
											dp = float(itag.split('=')[-1])
										if itag.startswith('AC='):
											ac = float(itag.split('=')[-1])
									newv =  [str(round(ac/dp, 3))]
								elif tag == 'PV':
									for itag in info_split:
										if itag.startswith(tag + '='):
											newv = [string.replace(itag.split('=')[-1], 'e', 'E')]
								else:
									for itag in info_split:
										if itag.startswith(tag + '='):
											newv = [itag.split('=')[-1]]
							else:
								for itag in info_split:
									if itag.startswith(tag + '='):
										newv = [itag.split('=')[-1]]
										break
						elif htag == 'FORMAT' or htag == 'format':
							if VariantCaller == 'SNVer':
								if tag == 'PL':
									newv = format_sample[format_split.index(tag)].split(',')
								else:
									if format_sample[format_split.index(tag)] == '.':
										newv = ['?']
									elif format_sample[format_split.index(tag)] == '1/0':
										newv = ['0/1']
									else:
										newv = [format_sample[format_split.index(tag)]]
							elif VariantCaller == 'VarScan' or VariantCaller == 'Varscan':
								if tag == 'PVAL':
									if ',' in format_sample[format_split.index(tag)]:
										try:
											newv = [str((format_sample[format_split.index(tag)]).replace(',', '.'))]
										except:
											newv = ['?']

								elif tag == 'FREQ':
									try:
										newv = [str(round(float(string.replace(format_sample[format_split.index(tag)].split('%')[0], ',', '.'))/100.0,3))]
									except:
										newv = ['?']
								else:
									newv = [format_sample[format_split.index(tag)]]
									if format_sample[format_split.index(tag)] == '.':
										newv = ['?']
									if format_sample[format_split.index(tag)] == '1/0':
										newv = ['0/1']
									if format_sample[format_split.index(tag)] == './.':
										newv = ['Undef']

							elif VariantCaller == 'GATK':
								if tag == 'AF':
									ad = float(format_sample[format_split.index('AD')].split(',')[1])
									dp = float(format_sample[format_split.index('AD')].split(',')[0]) + float(format_sample[format_split.index('AD')].split(',')[1])
									try:
										newv = [str(round(ad/dp,3))]
									except:
										newv = ['?']
								elif tag == 'PL':
									if ',' in format_sample[format_split.index(tag)]:
										oldv = format_sample[format_split.index(tag)].split(',')
										newv = [x if x != '.' else '?' for x in oldv]
									else:
										format_sample[format_split.index(tag)] = '.,.,.'
										oldv = format_sample[format_split.index(tag)].split(',')
										newv = [x if x != '.' else '?' for x in oldv]
								else:
									if format_sample[format_split.index(tag)] == '.':
										newv = ['?']
									if format_sample[format_split.index(tag)] == '1/0':
										newv = ['0/1']
									if format_sample[format_split.index(tag)] == './.':
										newv = ['Undef']
									else:
										oldv = [format_sample[format_split.index(tag)]]
										newv = [x if x != '.' else '?' for x in oldv]
							elif VariantCaller == 'FreeBayes':
								if tag == 'AF':
									ao = float(format_sample[format_split.index('AO')])
									ro = float(format_sample[format_split.index('RO')])
									try:
										newv = [str(round(ao/(ao+ro),3))]
									except:
										newv = ['?']
								elif tag == 'GL':
<<<<<<< HEAD
									try:
										newv = format_sample[format_split.index(tag)].split(',')
									except:
										newv = ['.','.','.']
=======
									oldv = format_sample[format_split.index(tag)].split(',')
									newv = [x if x != '.' else '?' for x in oldv]
>>>>>>> branch 'master' of https://github.com/urtism/CMG.git
								else:
									if format_sample[format_split.index(tag)] == '.':
										newv = ['?']
									if format_sample[format_split.index(tag)] == '1/0':
										newv = ['0/1']
									if format_sample[format_split.index(tag)] == './.':
										newv = ['Undef']
									else:
										oldv = [format_sample[format_split.index(tag)]]
										newv = [x if x != '.' else '?' for x in oldv]
							elif VariantCaller == 'Platypus':
								if tag == 'GL':
									newv = format_sample[format_split.index(tag)].split(',')
								elif tag == 'AF':
									ad = float(format_sample[format_split.index('NV')])
									dp = float(format_sample[format_split.index('NR')])
									try:
										newv = [str(round(ad/dp,3))]
									except:
										newv = ['?']
								else:
									if format_sample[format_split.index(tag)] == '.':
										newv = ['?']
									if format_sample[format_split.index(tag)] == '1/0':
										newv = ['0/1']
									if format_sample[format_split.index(tag)] == './.':
										newv = ['Undef']
									else:
										oldv = [format_sample[format_split.index(tag)]]
										newv = [x if x != '.' else '?' for x in oldv]
							else:
<<<<<<< HEAD
								try:
									if format_sample[format_split.index(tag)] == '.':
										newv = ['?']
									if format_sample[format_split.index(tag)] == '1/0':
										newv = ['0/1']
									else:
										newv = [format_sample[format_split.index(tag)]]
								except:
									print newv,VariantCaller
=======
								if format_sample[format_split.index(tag)] == '.':
									newv = ['?']
								if format_sample[format_split.index(tag)] == '1/0':
									newv = ['0/1']
								if format_sample[format_split.index(tag)] == './.':
									newv = ['Undef']
								else:
									oldv = [format_sample[format_split.index(tag)]]
									newv = [x if x != '.' else '?' for x in oldv]
>>>>>>> branch 'master' of https://github.com/urtism/CMG.git

						elif htag == 'IEVA' or htag == 'iEVA':
<<<<<<< HEAD
							
=======
>>>>>>> branch 'master' of https://github.com/urtism/CMG.git
							if tag == 'iSR' or tag == 'iSRL' or tag == 'iPNC' or tag == 'iSRU' or tag == 'iRM' or tag == 'iGC' or tag == 'iVC':
								for itag in info_split:
									#print itag
									if itag.startswith(tag + '='):
<<<<<<< HEAD
										
=======
>>>>>>> branch 'master' of https://github.com/urtism/CMG.git
										if tag == 'iPNC':
											newv = itag.split('=')[-1].split(',')

										else:
<<<<<<< HEAD
											newv = [itag.split('=')[-1]]

=======
											oldv = [itag.split('=')[-1]]
											newv = [x if x != '.' else '?' for x in oldv]
>>>>>>> branch 'master' of https://github.com/urtism/CMG.git
										break
							else:
								if tag == 'iAD' or tag == 'iSBD' or tag == 'iQual' or tag == 'iAMMQ' or tag == 'iAAS' or tag == 'iAXS' or tag == 'iAXS0' or tag == 'iAMQ0' or tag == 'iACR':
									oldv = format_sample[format_split.index(tag)].split(',')
									newv = [x if x != '.' else '?' for x in oldv]
								else:
<<<<<<< HEAD
									#print file,format_sample
									newv = [format_sample[format_split.index(tag)]]
									
						newvv=[]
						for w in newv:
							if w == '.':
								w='?'
							newvv += [w]
						var += newvv
=======
									oldv = [format_sample[format_split.index(tag)]]
									newv = [x if x != '.' else '?' for x in oldv]

						var += newv
>>>>>>> branch 'master' of https://github.com/urtism/CMG.git
		var+=[CLASS]
		if var[2] != 'NON TROVATA':
<<<<<<< HEAD
			vars += ['\t'.join(var)]

=======
			vars += [','.join(var)]	
>>>>>>> branch 'master' of https://github.com/urtism/CMG.git
					
	return  vars


if __name__ == '__main__':

	parser = argparse.ArgumentParser('\n\nQuesto tool estrae le varianti dai vcf a partire da una lista di varianti in tab delimited. BenchMark-Extractor prende la variante (deve avere nel tab delimited i campi CHR POS REF ALT ID e CLASS con ID riferito al codice paziente e CLASS riferito alla classe che puo essere PASS o FILTER (verificata in sanger o wt) e va a controllare in tutti i file forniti nel list a quale run appartiene e poi estra l inter riga del vcf in un nuovo file (UNO PER OGNI VARIANT CALLER) tab delimited con tutti i campi del vcf splittati correttamente. PREREQUISITO: data in ID paziente uguale a data in ID run nel nome del file. Es: ID_PAZ = 20160724_01_Conn mentre nome del file = 20150716_Cardio_iEVA_GATK.vcf --> Basta che 20160724 sia presente anche nel file.\n')
	parser.add_argument('-I','--input',help="File tab delimited contenente le varianti confermate in sanger")
<<<<<<< HEAD
	parser.add_argument('-L','--list',default=None,help="file contenente i path ai vcf file in cui cercare le varianti.I file devono essere nominati DATA_NUMPAZ_PANNELLO_VARIANTCALLER. se singlesample e DATA_PANNELLO_VARIANTCALLER. se multisample")
	parser.add_argument('-f','--features_list',help="file contenente gli attributi da estrarre dai file vcf. Un attributo per riga preceduto da INFO- se confenuto nelle info , da FORMAT- se contenuto nel formato o IEVA- se una feature di iEVA.")
=======
	parser.add_argument('-L','--list',default=None,help="file contenente i path ai vcf file in cui cercare le varianti.I file devono essere nominati DATA_NUMPAZ_PANNELLO_VARIANTCALLER.* se e' singlesample e DATA_PANNELLO_VARIANTCALLER.* se multisample")
	parser.add_argument('-f','--features_list',help="file contenente gli attributi da estrarre dai file vcf. Un attributo per riga preceduto da INFO- se e' confenuto nelle info , da FORMAT- se contenuto nel formato o IEVA- se e' una feature di iEVA.")
>>>>>>> branch 'master' of https://github.com/urtism/CMG.git
	parser.add_argument('-O','--outfile',help="file di output in tab delimited format.")

	global opts
	opts = parser.parse_args()
	global VariantCaller

	features = []
	with open(opts.features_list) as head:
		for elem in head:
			features += [elem.rstrip()]
	features = features + ['CLASS']

	ss_files = dict()
	ms_files = dict()

	with open(opts.list) as fl:
		for file in fl:
			try:
				data,num,pannello,VariantCaller= file.rstrip().split('/')[-1].split('.')[0].split('_')[:4]
				sample = '_'.join([data,num,pannello])
				#print sample,VariantCaller
				ss_files[sample] = file.rstrip()
			except:
				sample = '_'.join(file.rstrip().split('/')[-1].split('.')[0].split('_')[:1])
				VariantCaller = file.rstrip().split('/')[-1].split('.')[0].split('_')[-1]
				#print sample,VariantCaller
				ms_files[sample] = file.rstrip()
	
	samples = dict()
	with open(opts.input) as bm:
		for line in bm:
			line=line.rstrip()
			if line.startswith('Nome'):
				header = line.split('\t')
			else:
				Nome,id,chr,pos,ref,alt,HGVSc,HGVSp,GT_CLASS,CLASS=line.split('\t')
				var_id = '\t'.join([chr,pos,'.',ref,alt])
				if id in samples.keys():
					vars = samples[id]
					vars += [[var_id,CLASS]]
					samples[id] = vars
				else:
					samples[id] = [[var_id,CLASS]]

	out_bm = []
	for sample in samples.keys():
		#print ss_files.keys()
		if sample in ss_files.keys():

			vcf_sample = ss_files[sample]
			vars_sample = samples[sample]
			var = extract_feat_from_vcf(vcf_sample,vars_sample,features,sample)
			if var != []:
				out_bm += ['\n'.join(var)]
		else:
			data = sample.split('_')[0]

			for file in ms_files.keys():
				if file.startswith(data) and sample in samples.keys():
					vcf_sample = ms_files[file]
					vars_sample = samples[sample]
					var = extract_feat_from_vcf(vcf_sample,vars_sample,features,sample)
					if var != []:
						out_bm += ['\n'.join(var)]
	

	out = open(opts.outfile,'w')

	if VariantCaller == "SNVer" or VariantCaller == 'GATK' or VariantCaller == 'Samtools':
		features[features.index('FORMAT-PL'):features.index('FORMAT-PL')] = ['FORMAT-PL-0/0', 'FORMAT-PL-0/1','FORMAT-PL-1/1']
		del features[features.index('FORMAT-PL')]
	elif VariantCaller == 'Platypus' or VariantCaller == 'FreeBayes':
		features[features.index('FORMAT-GL'):features.index('FORMAT-GL')] = ['FORMAT-GL-0/0', 'FORMAT-GL-0/1','FORMAT-GL-1/1']
		del features[features.index('FORMAT-GL')]

	if "IEVA-iPNC" in features:
<<<<<<< HEAD
		features[features.index('IEVA-iPNC'):features.index('IEVA-iPNC')] = ['IEVA-iPNC AA','IEVA-iPNC AC', 'IEVA-iPNC AG', 'IEVA-iPNC AT', 'IEVA-iPNC CA', 'IEVA-iPNC CC', 'IEVA-iPNC CG', 'IEVA-iPNC CT', 'IEVA-iPNC GA', 'IEVA-iPNC GC', 'IEVA-iPNC GG', 'IEVA-iPNC GT', 'IEVA-iPNC TA', 'IEVA-iPNC TC', 'IEVA-iPNC TG', 'IEVA-iPNC TT']
=======
		features[features.index('IEVA-iPNC'):features.index('IEVA-iPNC')] = ['IEVA-iPNC-AA','IEVA-iPNC-AC', 'IEVA-iPNC-AG', 'IEVA-iPNC-AT', 'IEVA-iPNC-CA', 'IEVA-iPNC-CC', 'IEVA-iPNC-CG', 'IEVA-iPNC-CT', 'IEVA-iPNC-GA', 'IEVA-iPNC-GC', 'IEVA-iPNC-GG', 'IEVA-iPNC-GT', 'IEVA-iPNC-TA', 'IEVA-iPNC-TC', 'IEVA-iPNC-TG', 'IEVA-iPNC-TT']
>>>>>>> branch 'master' of https://github.com/urtism/CMG.git
		del features[features.index('IEVA-iPNC')]
	if "IEVA-iAD" in features:
		features[features.index('IEVA-iAD'):features.index('IEVA-iAD')] = ['IEVA-iAD-REF','IEVA-iAD-ALT']
		del features[features.index('IEVA-iAD')]
	if "IEVA-iSBD" in features:
		features[features.index('IEVA-iSBD'):features.index('IEVA-iSBD')] = ['IEVA-iSBD-RF','IEVA-iSBD-RR', 'IEVA-iSBD-AF' , 'IEVA-iSBD-AR']
		del features[features.index('IEVA-iSBD')]
	if "IEVA-iQual" in features:
		features[features.index('IEVA-iQual'):features.index('IEVA-iQual')] = ['IEVA-iQual-REF','IEVA-iQual-ALT']
		del features[features.index('IEVA-iQual')]
	if "IEVA-iAAS" in features:
		features[features.index('IEVA-iAMMQ'):features.index('IEVA-iAMMQ')] = ['IEVA-iAMMQ-REF','IEVA-iAMMQ-ALT']
		del features[features.index('IEVA-iAMMQ')]
	if "IEVA-iAAS" in features:
		features[features.index('IEVA-iAAS'):features.index('IEVA-iAAS')] = ['IEVA-iAAS-REF','IEVA-iAAS-ALT']
		del features[features.index('IEVA-iAAS')]
	if "IEVA-iAXS" in features:
		features[features.index('IEVA-iAXS'):features.index('IEVA-iAXS')] = ['IEVA-iAXS-REF','IEVA-iAXS-ALT']
		del features[features.index('IEVA-iAXS')]
	if "IEVA-iAXS0" in features:
		features[features.index('IEVA-iAXS0'):features.index('IEVA-iAXS0')] = ['IEVA-iAXS0-REF','IEVA-iAXS0-ALT']
		del features[features.index('IEVA-iAXS0')]
	if "IEVA-iAMQ0" in features:
		features[features.index('IEVA-iAMQ0'):features.index('IEVA-iAMQ0')] = ['IEVA-iAMQ0-REF','IEVA-iAMQ0-ALT']
		del features[features.index('IEVA-iAMQ0')]
	if "IEVA-iACR" in features:
		features[features.index('IEVA-iACR'):features.index('IEVA-iACR')] = ['IEVA-iACR-REF','IEVA-iACR-ALT']
		del features[features.index('IEVA-iACR')]


	out.write(','.join(['SAMPLE_ID']+['VAR_ID']+features)+'\n')

	for line in out_bm:
			out.write(line +'\n')

