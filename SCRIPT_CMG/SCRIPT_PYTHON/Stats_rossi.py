import argparse



class Var_confidance():
	truncation_mutations_lf=0
	inframe_indel_lf=0
	splice_lf=0
	miss_lf=0
	truncation_mutations_hf=0
	inframe_indel_hf=0
	splice_hf=0
	miss_hf=0
	var_lf=0
	var_hf=0
	intr_lf=0
	intr_hf=0
	syn_lf=0
	syn_hf=0
	others_lf=0
	others_hf=0
	vars_lf=[]
	vars_hf=[]


def make_stats(chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media,confidenza,HC,HC_vars_lf,HC_vars_hf,MC,MC_vars_lf,MC_vars_hf,LC,LC_vars_lf,LC_vars_hf):
	if confidenza == 'high':
						
		if float(AF_tum_media) < 0.05 and float(AF_tum_media) > 0.01:
			HC.var_lf += 1
			
			if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
				HC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				HC.truncation_mutations_lf += 1
			elif 'inframe' in Consequence:
				HC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				HC.inframe_indel_lf += 1
			elif 'splice' in Consequence:
				HC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				HC.splice_lf += 1
			elif 'missense' in Consequence:
				HC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				HC.miss_lf +=1
			elif 'intron' in Consequence:
				HC.intr_lf +=1
			elif 'synonymous' in Consequence:
				HC.syn_lf +=1
			else:
				HC.others_lf +=1


		elif float(AF_media) >= 0.05: 
			HC.var_hf += 1
			
			if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
				HC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				HC.truncation_mutations_hf += 1
			elif 'inframe' in Consequence:
				HC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				HC.inframe_indel_hf += 1
			elif 'splice' in Consequence:
				HC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				HC.splice_hf += 1
			elif 'missense' in Consequence:
				HC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				HC.miss_hf +=1
			elif 'intron' in Consequence:
				HC.intr_hf +=1
			elif 'synonymous' in Consequence:
				HC.syn_hf +=1
			else:
				HC.others_hf +=1

	elif confidenza == 'medium':
		if float(AF_media) < 0.05 and float(AF_media) > 0.01:
			MC.var_lf += 1
			if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
				MC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				MC.truncation_mutations_lf += 1
			elif 'inframe' in Consequence:
				MC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				MC.inframe_indel_lf += 1
			elif 'splice' in Consequence:
				MC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				MC.splice_lf += 1
			elif 'missense' in Consequence:
				MC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				MC.miss_lf +=1
			elif 'intron' in Consequence:
				MC.intr_lf +=1
			elif 'synonymous' in Consequence:
				MC.syn_lf +=1
			else:
				MC.others_lf +=1

		elif float(AF_media) >= 0.05: 
			MC.var_hf += 1
			if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
				MC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				MC.truncation_mutations_hf += 1
			elif 'inframe' in Consequence:
				MC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				MC.inframe_indel_hf += 1
			elif 'splice' in Consequence:
				MC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				MC.splice_hf += 1
			elif 'missense' in Consequence:
				MC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				MC.miss_hf +=1
			elif 'intron' in Consequence:
				MC.intr_hf +=1
			elif 'synonymous' in Consequence:
				MC.syn_hf +=1
			else:
				MC.others_hf +=1

	elif confidenza == 'low':
		if float(AF_media) < 0.05 and float(AF_media) > 0.01:
			LC.var_lf += 1
			if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
				LC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				LC.truncation_mutations_lf += 1
			elif 'inframe' in Consequence:
				LC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				LC.inframe_indel_lf += 1
			elif 'splice' in Consequence:
				LC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				LC.splice_lf += 1
			elif 'missense' in Consequence:
				LC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				LC.miss_lf += 1
			elif 'intron' in Consequence:
				LC.intr_lf += 1
			elif 'synonymous' in Consequence:
				LC.syn_lf += 1
			else:
				LC.others_lf += 1

		elif float(AF_media) >= 0.05: 
			LC.var_hf += 1
			if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
				LC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				LC.truncation_mutations_hf += 1
			elif 'inframe' in Consequence:
				LC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				LC.inframe_indel_hf += 1
			elif 'splice' in Consequence:
				LC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				LC.splice_hf += 1
			elif 'missense' in Consequence:
				LC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_media])]
				LC.miss_hf +=1
			elif 'intron' in Consequence:
				LC.intr_hf +=1
			elif 'synonymous' in Consequence:
				LC.syn_hf +=1
			else:
				LC.others_hf +=1

	return HC,HC_vars_lf,HC_vars_hf,MC,MC_vars_lf,MC_vars_hf,LC,LC_vars_lf,LC_vars_hf


def assegna_conf_som(somVars,somVard,filterMu):
	
	varcal=0
	som=0
	conf='?'
	if somVars != '.':
		varcal += 1
		if somVars == '1':
			som += 1
	if somVard != '.':
		varcal += 1
		if somVard == '1':
			som += 1
	if filterMu != '.':
		varcal += 1
		if filterMu == 'PASS':
			som += 1

	print varcal,som

	if varcal > 1 and som >=2 :
		conf = 'high'
	elif varcal > 1 and som == 1:
		conf = 'medium'
	elif varcal == 1 and som == 1:
		conf = 'low'

	return conf


def assegna_conf_cirr(gt_tum,gt_cirr):
	
	OMO_tum=gt_tum.count('0/0')
	ET_tum=gt_tum.count('0/1')+gt_tum.count('1/0')
	OMO_cirr=gt_cirr.count('0/1')+ gt_cirr.count('1/0')
	ET_cirr=gt_cirr.count('0/0')

	if ET_cirr >= 2 and OMO_tum >= 2 and OMO_cirr == 0 and ET_tum == 0 :
		conf = 'high'
	elif ET_tum > 0 and  ET_cirr >= 2 :
		conf = 'medium'
	else:
		conf = 'low'

	return conf



if __name__ == '__main__':

	parser = argparse.ArgumentParser('Genera un report delle varianti trovate in un campione o in un gruppo di campioni')
	parser.add_argument('--tsv',help="file  tsv di varianti somatiche con il tag 'TESSUTO'",default=None)
	parser.add_argument('-o','--out',help="path di out in cui stampare le statistiche")
	parser.add_argument('--list',help="lista di file  tsv di varianti somatiche con il tag 'TESSUTO' separati da ,",default=None)
	parser.add_argument('-t','--tess',help="tessuto su cui fare le statistiche",default=None)
	
	global opts
	opts = parser.parse_args()


	if opts.tsv is not None:
		tsv=open(opts.tsv,'r')
		out=open(opts.out,'w')

		HC=Var_confidance()
		MC=Var_confidance()
		LC=Var_confidance()

		HC_vars_lf=[]
		HC_vars_hf=[]
		MC_vars_lf=[]
		MC_vars_hf=[]
		LC_vars_lf=[]
		LC_vars_hf=[]

		for line in tsv:
			line=line.rstrip()
			if line.startswith('CHROM'):
				header=line.split('\t')
			else:
				chr,pos,id,ref,alt=line.split('\t')[:5]
				
				tess = line.split('\t')[header.index('TESSUTO')]
				SomaticVarscan = line.split('\t')[header.index('SomaticVarscan')]
				SomaticVardict = line.split('\t')[header.index('SomaticVardict')]
				FILTER_Mutect = line.split('\t')[header.index('FILTER_Mutect')]
				AF_tum_media = line.split('\t')[header.index('AF_tum_media')]
				DP_tum_media = line.split('\t')[header.index('DP_tum_media')]
				AO_tum_media = line.split('\t')[header.index('AO_tum_media')]

				AF_cirr_media = line.split('\t')[header.index('AF_norm_media')]
				try:
					DP_cirr_media = float(line.split('\t')[header.index('DP_norm_media')])
				except:
					DP_cirr_media = 0.0
				try:
					AO_cirr_media = float(line.split('\t')[header.index('AO_norm_media')])
				except:
					AO_cirr_media = 0.0

				gt_cirr=[line.split('\t')[header.index('GT_n_Varscan')],line.split('\t')[header.index('GT_n_Vardict')],line.split('\t')[header.index('GT_n_Mutect')]]
				gt_tum=[line.split('\t')[header.index('GT_t_Varscan')],line.split('\t')[header.index('GT_t_Vardict')],line.split('\t')[header.index('GT_t_Mutect')]]
				gt_perif=line.split('\t')[header.index('GT_PERIF')].split(',')


				Delta_media = line.split('\t')[header.index('Delta_media')]
				Delta_perc_media = line.split('\t')[header.index('Delta_perc_media')]
				Consequence = line.split('\t')[header.index('Consequence')]
				gene = line.split('\t')[header.index('SYMBOL')]
				HGVSc = line.split('\t')[header.index('HGVSc')]
				HGVSp = line.split('\t')[header.index('HGVSp')]


				
				if opts.tess =='TUMORE' and float(DP_tum_media) >= 400.0 and float(AO_tum_media) >= 10 and float(Delta_perc_media) == 100 and tess == opts.tess:

					confidenza = assegna_conf_som(SomaticVarscan,SomaticVardict,FILTER_Mutect)
					HC,HC_vars_lf,HC_vars_hf,MC,MC_vars_lf,MC_vars_hf,LC,LC_vars_lf,LC_vars_hf = make_stats(chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media,confidenza,HC,HC_vars_lf,HC_vars_hf,MC,MC_vars_lf,MC_vars_hf,LC,LC_vars_lf,LC_vars_hf)

				elif opts.tess =='CIRR+TUMORE' and float(DP_cirr_media) >= 400.0 and float(AO_cirr_media) >= 10 and '0/0' in gt_tum and '0/1' in gt_cirr and tess == opts.tess:
					#print LC_vars_lf
					confidenza = assegna_conf_cirr(gt_tum,gt_cirr)
					HC,HC_vars_lf,HC_vars_hf,MC,MC_vars_lf,MC_vars_hf,LC,LC_vars_lf,LC_vars_hf = make_stats(chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_cirr_media,confidenza,HC,HC_vars_lf,HC_vars_hf,MC,MC_vars_lf,MC_vars_hf,LC,LC_vars_lf,LC_vars_hf)
					#print LC_vars_lf
				elif opts.tess =='PERIFERICO' and float(DP_tum_media) >= 400.0 and tess == opts.tess:

					confidenza = 'high'
					HC,HC_vars_lf,HC_vars_hf,MC,MC_vars_lf,MC_vars_hf,LC,LC_vars_lf,LC_vars_hf = make_stats(chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media,confidenza,HC,HC_vars_lf,HC_vars_hf,MC,MC_vars_lf,MC_vars_hf,LC,LC_vars_lf,LC_vars_hf)


					# if confidenza == 'high':
						
					# 	if float(AF_tum_media) < 0.05 and float(AF_tum_media) > 0.01:
					# 		HC.var_lf += 1
							
					# 		if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
					# 			HC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			HC.truncation_mutations_lf += 1
					# 		elif 'inframe' in Consequence:
					# 			HC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			HC.inframe_indel_lf += 1
					# 		elif 'splice' in Consequence:
					# 			HC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			HC.splice_lf += 1
					# 		elif 'missense' in Consequence:
					# 			HC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			HC.miss_lf +=1
					# 		elif 'intron' in Consequence:
					# 			HC.intr_lf +=1
					# 		elif 'synonymous' in Consequence:
					# 			HC.syn_lf +=1
					# 		else:
					# 			HC.others_lf +=1


					# 	elif float(AF_tum_media) >= 0.05: 
					# 		HC.var_hf += 1
							
					# 		if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
					# 			HC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			HC.truncation_mutations_hf += 1
					# 		elif 'inframe' in Consequence:
					# 			HC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			HC.inframe_indel_hf += 1
					# 		elif 'splice' in Consequence:
					# 			HC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			HC.splice_hf += 1
					# 		elif 'missense' in Consequence:
					# 			HC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			HC.miss_hf +=1
					# 		elif 'intron' in Consequence:
					# 			HC.intr_hf +=1
					# 		elif 'synonymous' in Consequence:
					# 			HC.syn_hf +=1
					# 		else:
					# 			HC.others_hf +=1

					# elif confidenza == 'medium':
					# 	if float(AF_tum_media) < 0.05 and float(AF_tum_media) > 0.01:
					# 		MC.var_lf += 1
					# 		if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
					# 		MC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			MC.truncation_mutations_lf += 1
					# 		elif 'inframe' in Consequence:
					# 		MC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			MC.inframe_indel_lf += 1
					# 		elif 'splice' in Consequence:
					# 		MC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			MC.splice_lf += 1
					# 		elif 'missense' in Consequence:
					# 		MC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			MC.miss_lf +=1
					# 		elif 'intron' in Consequence:
					# 			MC.intr_lf +=1
					# 		elif 'synonymous' in Consequence:
					# 			MC.syn_lf +=1
					# 		else:
					# 			MC.others_lf +=1

					# 	elif float(AF_tum_media) >= 0.05: 
					# 		MC.var_hf += 1
					# 		if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
					# 		MC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			MC.truncation_mutations_hf += 1
					# 		elif 'inframe' in Consequence:
					# 		MC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			MC.inframe_indel_hf += 1
					# 		elif 'splice' in Consequence:
					# 		MC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			MC.splice_hf += 1
					# 		elif 'missense' in Consequence:
					# 		MC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			MC.miss_hf +=1
					# 		elif 'intron' in Consequence:
					# 			MC.intr_hf +=1
					# 		elif 'synonymous' in Consequence:
					# 			MC.syn_hf +=1
					# 		else:
					# 			MC.others_hf +=1

					# elif confidenza == 'low':
					# 	if float(AF_tum_media) < 0.05 and float(AF_tum_media) > 0.01:
					# 		LC.var_lf += 1
					# 		if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
					# 		LC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			LC.truncation_mutations_lf += 1
					# 		elif 'inframe' in Consequence:
					# 		LC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			LC.inframe_indel_lf += 1
					# 		elif 'splice' in Consequence:
					# 		LC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			LC.splice_lf += 1
					# 		elif 'missense' in Consequence:
					# 		LC_vars_lf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			LC.miss_lf += 1
					# 		elif 'intron' in Consequence:
					# 			LC.intr_lf += 1
					# 		elif 'synonymous' in Consequence:
					# 			LC.syn_lf += 1
					# 		else:
					# 			LC.others_lf += 1

					# 	elif float(AF_tum_media) >= 0.05: 
					# 		LC.var_hf += 1
					# 		if 'stop_gained' in Consequence or 'frameshift_variant' in Consequence:
					# 		LC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			LC.truncation_mutations_hf += 1
					# 		elif 'inframe' in Consequence:
					# 		LC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			LC.inframe_indel_hf += 1
					# 		elif 'splice' in Consequence:
					# 		LC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			LC.splice_hf += 1
					# 		elif 'missense' in Consequence:
					# 		LC_vars_hf+=['\t'.join([chr,pos,id,ref,alt,Consequence,gene,HGVSc,HGVSp,AF_tum_media])]
					# 			LC.miss_hf +=1
					# 		elif 'intron' in Consequence:
					# 			LC.intr_hf +=1
					# 		elif 'synonymous' in Consequence:
					# 			LC.syn_hf +=1
					# 		else:
					# 			LC.others_hf +=1
		

		#print HC.vars_lf
		#print MC.vars_lf		
		out.write('HIGH_CONFIDENCE\n')
		out.write('Tipo variante\tLow Freq\tHigh Freq\n')
		out.write('truncation_mutations\t'+str(HC.truncation_mutations_lf)+'\t'+ str(HC.truncation_mutations_hf)+'\n')		
		out.write('splicing_site_mutations\t'+str(HC.splice_lf)+'\t'+ str(HC.splice_hf)+'\n')
		out.write('inframe_indel\t'+str(HC.inframe_indel_lf)+'\t'+ str(HC.inframe_indel_hf)+'\n')
		out.write('missense\t'+str(HC.miss_lf)+'\t'+ str(HC.miss_hf)+'\n')
		out.write('intronic\t'+str(HC.intr_lf)+'\t'+ str(HC.intr_hf)+'\n')
		out.write('synonimus\t'+str(HC.syn_lf)+'\t'+ str(HC.syn_hf)+'\n')
		out.write('other\t'+str(HC.others_lf)+'\t'+ str(HC.others_hf)+'\n')
		out.write('varianti a freq > 0.01 e < 0.05:\t'+str(HC.var_lf) +'\n')
		for v in HC_vars_lf:
			out.write(v+'\n')
		out.write('varianti a freq > 0.05:\t'+str(HC.var_hf) +'\n')
		for v in HC_vars_hf:
			out.write(v+'\n')

		out.write('\nMEDIUM_CONFIDENCE:\n')
		out.write('Tipo variante\tLow Freq\tHigh Freq\n')
		out.write('truncation_mutations\t'+str(MC.truncation_mutations_lf)+'\t'+ str(MC.truncation_mutations_hf)+'\n')		
		out.write('splicing_site_mutations\t'+str(MC.splice_lf)+'\t'+ str(MC.splice_hf)+'\n')
		out.write('inframe_indel\t'+str(MC.inframe_indel_lf)+'\t'+ str(MC.inframe_indel_hf)+'\n')
		out.write('missense\t'+ str(MC.miss_lf)+'\t'+ str(MC.miss_hf)+'\n')
		out.write('intronic\t'+str(MC.intr_lf)+'\t'+ str(MC.intr_hf)+'\n')
		out.write('synonimus\t'+str(MC.syn_lf)+'\t'+ str(MC.syn_hf)+'\n')
		out.write('other\t'+str(MC.others_lf)+'\t'+ str(MC.others_hf)+'\n')
		out.write('varianti a freq > 0.01 e < 0.05:\t'+str(MC.var_lf) +'\n')
		for v in MC_vars_lf:
			out.write(v+'\n')
		out.write('varianti a freq > 0.05:\t'+str(MC.var_hf) +'\n')
		for v in MC_vars_hf:
			out.write(v+'\n')

		out.write('\nLOW_CONFIDENCE:\n')
		out.write('Tipo variante\tLow Freq\tHigh Freq\n')
		out.write('truncation_mutations\t'+ str(LC.truncation_mutations_lf)+'\t'+ str(LC.truncation_mutations_hf)+'\n')		
		out.write('splicing_site_mutations\t'+ str(LC.splice_lf)+'\t'+ str(LC.splice_hf)+'\n')
		out.write('inframe_indel\t'+ str(LC.inframe_indel_lf)+'\t'+ str(LC.inframe_indel_hf)+'\n')
		out.write('missense\t'+str(LC.miss_lf)+'\t'+ str(LC.miss_hf)+'\n')
		out.write('intronic\t'+str(LC.intr_lf)+'\t'+ str(LC.intr_hf)+'\n')
		out.write('synonimus\t'+str(LC.syn_lf)+'\t'+ str(LC.syn_hf)+'\n')
		out.write('other\t'+str(LC.others_lf)+'\t'+ str(LC.others_hf)+'\n')
		out.write('varianti a > 0.01 e < 0.05:\t'+str(HC.var_lf) +'\n')
		for v in LC_vars_lf:
			out.write(v+'\n')
		out.write('varianti a freq > 0.05:\t'+str(LC.var_hf) +'\n')
		for v in LC_vars_hf:
			out.write(v+'\n')


	elif opts.list is not None:


		out=open(opts.out,'w')

		HC=Var_confidance()
		MC=Var_confidance()
		LC=Var_confidance()

		HC_vars_lf=[]
		HC_vars_hf=[]
		MC_vars_lf=[]
		MC_vars_hf=[]
		LC_vars_lf=[]
		LC_vars_hf=[]
		
		for sample in (opts.list).split(','):			
			tsv=[]
			for a in open(sample,'r').readlines():
				tsv += [a.strip()]
				if a.startswith('HIGH_CONFIDENCE'):
					index_HC = tsv.index(a.rstrip())
				elif a.startswith('MEDIUM_CONFIDENCE'):
					index_MC = tsv.index(a.rstrip())
				elif a.startswith('LOW_CONFIDENCE'):
					index_LC = tsv.index(a.rstrip())	



			tsv_HC = tsv[:index_MC]
			tsv_MC = tsv[index_MC:index_LC]
			tsv_LC = tsv[index_LC:]

			HC.truncation_mutations_lf += int(tsv_HC[2].split('\t')[1])
			HC.splice_lf += int(tsv_HC[3].split('\t')[1])
			HC.inframe_indel_lf += int(tsv_HC[4].split('\t')[1])
			HC.miss_lf += int(tsv_HC[5].split('\t')[1])
			HC.intr_lf += int(tsv_HC[6].split('\t')[1])
			HC.syn_lf += int(tsv_HC[7].split('\t')[1])
			HC.others_lf += int(tsv_HC[8].split('\t')[1])

			HC.truncation_mutations_hf += int(tsv_HC[2].split('\t')[2])
			HC.splice_hf += int(tsv_HC[3].split('\t')[2])
			HC.inframe_indel_hf += int(tsv_HC[4].split('\t')[2])
			HC.miss_hf += int(tsv_HC[5].split('\t')[2])
			HC.intr_hf += int(tsv_HC[6].split('\t')[2])
			HC.syn_hf += int(tsv_HC[7].split('\t')[2])
			HC.others_hf += int(tsv_HC[8].split('\t')[2])

			for a in tsv_HC:
				if a.startswith('varianti a freq > 0.01 e < 0.05:'):
					index_lf = tsv_HC.index(a)
				elif a.startswith('varianti a freq > 0.05:'):
					index_hf = tsv_HC.index(a)	

			
			i = index_lf
			while i < index_hf-1:
				i+=1
				HC_vars_lf+=[tsv_HC[i].rstrip()]

			
			for v in tsv_HC[index_hf+1:-1]:
				HC_vars_hf+=[v.rstrip()]
			

			MC.truncation_mutations_lf += int(tsv_MC[2].split('\t')[1])
			MC.splice_lf += int(tsv_MC[3].split('\t')[1])
			MC.inframe_indel_lf += int(tsv_MC[4].split('\t')[1])
			MC.miss_lf += int(tsv_MC[5].split('\t')[1])
			MC.intr_lf += int(tsv_MC[6].split('\t')[1])
			MC.syn_lf += int(tsv_MC[7].split('\t')[1])
			MC.others_lf += int(tsv_MC[8].split('\t')[1])

			MC.truncation_mutations_hf += int(tsv_MC[2].split('\t')[2])
			MC.splice_hf += int(tsv_MC[3].split('\t')[2])
			MC.inframe_indel_hf += int(tsv_MC[4].split('\t')[2])
			MC.miss_hf += int(tsv_MC[5].split('\t')[2])
			MC.intr_hf += int(tsv_MC[6].split('\t')[2])
			MC.syn_hf += int(tsv_MC[7].split('\t')[2])
			MC.others_hf += int(tsv_MC[8].split('\t')[2])
			
			for a in tsv_MC:
				if a.startswith('varianti a freq > 0.01 e < 0.05:'):
					index_lf = tsv_MC.index(a)
				elif a.startswith('varianti a freq > 0.05:'):
					index_hf = tsv_MC.index(a)	
			
			i = index_lf
			while i < index_hf-1:
				i+=1
				MC_vars_lf+=[tsv_MC[i].rstrip()]

			for v in tsv_MC[index_hf+1:-1]:
				MC_vars_hf+=[v.rstrip()]


			LC.truncation_mutations_lf += int(tsv_LC[2].split('\t')[1])
			LC.splice_lf += int(tsv_LC[3].split('\t')[1])
			LC.inframe_indel_lf += int(tsv_LC[4].split('\t')[1])
			LC.miss_lf += int(tsv_LC[5].split('\t')[1])
			LC.intr_lf += int(tsv_LC[6].split('\t')[1])
			LC.syn_lf += int(tsv_LC[7].split('\t')[1])
			LC.others_lf += int(tsv_LC[8].split('\t')[1])

			LC.truncation_mutations_hf += int(tsv_LC[2].split('\t')[2])
			LC.splice_hf += int(tsv_LC[3].split('\t')[2])
			LC.inframe_indel_hf += int(tsv_LC[4].split('\t')[2])
			LC.miss_hf += int(tsv_LC[5].split('\t')[2])
			LC.intr_hf += int(tsv_LC[6].split('\t')[2])
			LC.syn_hf += int(tsv_LC[7].split('\t')[2])
			LC.others_hf += int(tsv_LC[8].split('\t')[2])

			for a in tsv_LC:
				if a.startswith('varianti a freq > 0.01 e < 0.05:'):
					index_lf = tsv_LC.index(a)
				elif a.startswith('varianti a freq > 0.05:'):
					index_hf = tsv_LC.index(a)	
			i = index_lf
			while i < index_hf-1:
				i+=1
				LC_vars_lf+=[tsv_LC[i].rstrip()]
				
			
			for v in tsv_LC[index_hf+1:]:
				LC_vars_hf+=[v.rstrip()]




		out.write('HIGH_CONFIDENCE\n')
		out.write('Tipo variante\tLow Freq\tHigh Freq\n')
		out.write('truncation_mutations\t'+str(HC.truncation_mutations_lf)+'\t'+ str(HC.truncation_mutations_hf)+'\n')		
		out.write('splicing_site_mutations\t'+str(HC.splice_lf)+'\t'+ str(HC.splice_hf)+'\n')
		out.write('inframe_indel\t'+str(HC.inframe_indel_lf)+'\t'+ str(HC.inframe_indel_hf)+'\n')
		out.write('missense\t'+str(HC.miss_lf)+'\t'+ str(HC.miss_hf)+'\n')
		out.write('intronic\t'+str(HC.intr_lf)+'\t'+ str(HC.intr_hf)+'\n')
		out.write('synonimus\t'+str(HC.syn_lf)+'\t'+ str(HC.syn_hf)+'\n')
		out.write('other\t'+str(HC.others_lf)+'\t'+ str(HC.others_hf)+'\n')
		out.write('varianti a freq > 0.01 e < 0.05:\t'+str(HC.var_lf) +'\n')
		for v in HC_vars_lf:
			out.write(v+'\n')
		out.write('varianti a freq > 0.05:\t'+str(HC.var_hf) +'\n')
		for v in HC_vars_hf:
			out.write(v+'\n')

		out.write('\nMEDIUM_CONFIDENCE:\n')
		out.write('Tipo variante\tLow Freq\tHigh Freq\n')
		out.write('truncation_mutations\t'+str(MC.truncation_mutations_lf)+'\t'+ str(MC.truncation_mutations_hf)+'\n')		
		out.write('splicing_site_mutations\t'+str(MC.splice_lf)+'\t'+ str(MC.splice_hf)+'\n')
		out.write('inframe_indel\t'+str(MC.inframe_indel_lf)+'\t'+ str(MC.inframe_indel_hf)+'\n')
		out.write('missense\t'+ str(MC.miss_lf)+'\t'+ str(MC.miss_hf)+'\n')
		out.write('intronic\t'+str(MC.intr_lf)+'\t'+ str(MC.intr_hf)+'\n')
		out.write('synonimus\t'+str(MC.syn_lf)+'\t'+ str(MC.syn_hf)+'\n')
		out.write('other\t'+str(MC.others_lf)+'\t'+ str(MC.others_hf)+'\n')
		out.write('varianti a freq > 0.01 e < 0.05:\t'+str(MC.var_lf) +'\n')
		for v in MC_vars_lf:
			out.write(v+'\n')
		out.write('varianti a freq > 0.05:\t'+str(MC.var_hf) +'\n')
		for v in MC_vars_hf:
			out.write(v+'\n')

		out.write('\nLOW_CONFIDENCE:\n')
		out.write('Tipo variante\tLow Freq\tHigh Freq\n')
		out.write('truncation_mutations\t'+ str(LC.truncation_mutations_lf)+'\t'+ str(LC.truncation_mutations_hf)+'\n')		
		out.write('splicing_site_mutations\t'+ str(LC.splice_lf)+'\t'+ str(LC.splice_hf)+'\n')
		out.write('inframe_indel\t'+ str(LC.inframe_indel_lf)+'\t'+ str(LC.inframe_indel_hf)+'\n')
		out.write('missense\t'+str(LC.miss_lf)+'\t'+ str(LC.miss_hf)+'\n')
		out.write('intronic\t'+str(LC.intr_lf)+'\t'+ str(LC.intr_hf)+'\n')
		out.write('synonimus\t'+str(LC.syn_lf)+'\t'+ str(LC.syn_hf)+'\n')
		out.write('other\t'+str(LC.others_lf)+'\t'+ str(LC.others_hf)+'\n')
		out.write('varianti a freq > 0.01 e < 0.05:\t'+str(HC.var_lf) +'\n')
		for v in LC_vars_lf:
			out.write(v+'\n')
		out.write('varianti a freq > 0.05:\t'+str(LC.var_hf) +'\n')
		for v in LC_vars_hf:
			out.write(v+'\n')