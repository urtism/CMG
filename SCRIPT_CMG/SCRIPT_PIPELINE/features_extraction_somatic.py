import re
import string
import argparse
import sys
import statistics
import scipy.stats as stats


''' //////////// CLASSI ////////////'''

class Caller():
	GT=''
	AO_t=''
	AO_n=''
	RO_t=''		
	RO_n=''
	SS=''
	AO_f=''
	AO_r=''
	DP_f=''
	DP_r=''
	DP_n=0
	DP_t=0
	DP=''
	QB=''
	QB_n=''
	Somatic=''
	Call=''
	AF_t=''
	AF_n=''
	STRBIAS=''
	STRBIAS_n=''
	StrandBias=''
	MULTIALLELE=0
	LOH=''
	FILTER=''

class Mutect(Caller):
	t_ad=''
	t_af=''
	t_ALT_F1R2=''
	t_ALT_F1R1=''
	t_FOXOG=''
	t_GQ=''
	t_PGT=''
	t_PID=''
	t_PL=''
	t_QSS=''
	t_REF_F1R2=''
	t_REF_F1R1=''

	n_ad=''
	n_af=''
	n_ALT_F1R2=''
	n_ALT_F1R1=''
	n_FOXOG=''
	n_GQ=''
	n_PGT=''
	n_PID=''
	n_PL=''
	n_QSS=''
	n_REF_F1R2=''
	n_REF_F1R1=''

	HCNT='.'
	MAX_ED='.'
	MIN_ED='.'
	NLOD='.'
	PON='.'
	RPA='.'
	RU='.'
	STR='.'
	TLOD='.'

	


class Vardict(Caller):
	StrandBias=''
	ODDRATIO=''
	SSF=''
	
	t_PMEAN=''
	t_PSTD=''
	t_QSTD=''
	t_SBF=''
	t_MQ=''
	t_SN=''
	t_HIAF=''
	t_NM=''

	n_PMEAN=''
	n_PSTD=''
	n_QSTD=''
	n_SBF=''
	n_MQ=''
	n_SN=''
	n_HIAF=''
	n_NM=''

	SHIFT3=''
	MSI=''
	MSILEN=''
	SOR=''
	LSEQ=''
	RSEQ=''
	STATUS=''

	
class Varscan(Caller):
	pass
	

class Features():
	GT_Varscan='.'
	GT_Vardict='.'
	GT_Mutect='.'
	
	DP=float(0)
	MBQT='.'
	QB_Mutect='.'
	QB_Varscan='.'
	QB_Vardict='.'
	
	AF_media="."
	Delta_perc_media='.'
	STRBIAS_media='.'
	delta_media='.'
	MBQT_median='.'
	DP_median='.'
	delta_median='.'
	Delta_perc_median='.'
	AF_median='.'
	STRBIAS_median= '.'
	DPn_t_median='.'
	DPn_t_media='.'

	
	DP_n_t_Mutect='.'
	DP_n_t_Varscan='.'
	DP_n_t_Vardict='.'
	
	AF_t_Mutect='0'
	AF_n_Mutect='.'
	AF_t_Varscan='0'
	AF_n_Varscan='.'
	AF_t_Vardict='0'
	AF_n_Vardict='.'
	
	Delta_Mutect='.'
	Delta_Varscan='.'
	Delta_Vardict='.'
	
	Delta_perc_Mutect='.'
	Delta_perc_Varscan='.'
	Delta_perc_Vardict='.'
	
	SomaticMutect='.'
	SomaticVarscan='.'
	SomaticVardict='.'
	
	CallMutect='0'
	CallVarscan='0'
	CallVardict='0'
	
	SS_Mutect='.'
	SS_Varscan='.'
	SS_Vardict='.'
	
	tumor_lod_Mutect='.'
	normal_lod_Mutect='.'
	delta_lod_Mutect='.'
	t_n_lod_Mutect='.'

	STRBIAS_Mutect='.'
	STRBIAS_Varscan='.'
	STRBIAS_Vardict='.'
	BIAS_Vardict='.'
	
	ODDRATIO_Vardict='.'
	SBF_Vardict='.'
	
	LOH_Mutect='.'
	LOH_Varscan='.'
	LOH_Vardict='.'

	INFO_Mutect='.'
	INFO_Varscan='.'
	INFO_Vardict='.'

	AF_norm_media="."
	STRBIAS_norm_media='.'
	MBQN_median='.'
	MBQN_media='.'
	DP_Norm_median='.'
	AF_norm_median='.'
	STRBIAS_norm_median= '.'

	AO_tum_media='.'
	RO_tum_media='.'
	AO_tum_median='.'
	RO_tum_median='.'

	AO_norm_media='.'
	RO_norm_media='.'
	AO_norm_median='.'
	RO_norm_median='.'

	#features solo di Mutect
	ALT_F1R2_Mutect='.'
	ALT_F2R1_Mutect='.'
	FOXOG_Mutect='.'
	GQ_Mutect='.'
	PGT_Mutect='.'
	PID_Mutect='.'
	PL_Mutect='.'
	QSS_Mutect='.'
	REF_F1R2_Mutect='.'
	REF_F2R1_Mutect='.'

	HCNT_Mutect='.'
	MAX_ED_Mutect='.'
	MIN_ED_Mutect='.'
	PON_Mutect='.'
	RPA_Mutect='.'
	RU_Mutect='.'
	STR_Mutect='.'
	FILTER_Mutect='.'

	#features solo di VarDict
	SHIFT3_Vardict='.'
	MSI_Vardict='.'
	MSILEN_Vardict='.'
	SOR_Vardict='.'
	LSEQ_Vardict='.'
	RSEQ_Vardict='.'
	STATUS_Vardict='.'

	PMEAN_Vardict='.'
	PSTD_Vardict='.'
	QSTD_Vardict='.'
	SBF_Vardict='.'
	MQ_Vardict='.'
	SN_Vardict='.'
	HIAF_Vardict='.'
	NM_Vardict='.'




''' //////////// FUNZIONI ////////////'''


def get_info_Mutect(chrom,pos,ref,alt,filter,info,format,tumor,normal,Mutect):
	'''estrae le informazioni dal vcf di Mutect'''
	Mutect.GT=tumor[format.index('GT')]
	Mutect.AO_t=float((tumor[format.index('AD')]).split(',')[1])
	Mutect.RO_t=float((tumor[format.index('AD')]).split(',')[0])
	Mutect.DP_t=Mutect.AO_t+Mutect.RO_t
	Mutect.DP=Mutect.AO_t+Mutect.RO_t	
	Mutect.DP_n=float((normal[format.index('AD')]).split(',')[0]) + float((normal[format.index('AD')]).split(',')[1])
	

	if normal is not 'null' and Mutect.DP_n is not '0':
		try:
			Mutect.AO_n=float((normal[format.index('AD')]).split(',')[1])
		except:
			Mutect.AO_n=float(0)
		try:
			Mutect.RO_n=float((normal[format.index('AD')]).split(',')[0])
		except:
			Mutect.RO_n=float(0)
		try:
			Mutect.AF_n=Mutect.AO_n/Mutect.DP_n
		except:
			Mutect.AF_n='.'
	else:
		Mutect.DP_n=0
		Mutect.AF_n='.'
		Mutect.AO_n='.'
		Mutect.RO_n='.'
	
	Mutect.STR='0'
	for ind in info:
		if ind.startswith("HCNT="):
			Mutect.HCNT=ind.split('=')[1]
		if ind.startswith("MAX_ED="):
			Mutect.MAX_ED=ind.split('=')[1]
		if ind.startswith("MIN_ED="):	
			Mutect.MIN_ED=ind.split('=')[1]
		if ind.startswith("NLOD="):	
			Mutect.normal_lod=float(ind.split('=')[1])
		if ind.startswith("PON="):	
			Mutect.PON=ind.split('=')[1]
		if ind.startswith("RPA="):	
			Mutect.RPA=ind.split('=')[1]
		if ind.startswith("RU="):	
			Mutect.RU=ind.split('=')[1]
		if ind.startswith("STR="):	
			Mutect.STR='1'
		if ind.startswith("TLOD="):	
			Mutect.tumor_lod=float(ind.split('=')[1])

	Mutect.AO_f=float(tumor[format.index('ALT_F1R2')])
	Mutect.AO_r=float(tumor[format.index('ALT_F2R1')])
	Mutect.RO_f=float(tumor[format.index('REF_F1R2')])
	Mutect.RO_r=float(tumor[format.index('REF_F2R1')])

	Mutect.DP_f=float(Mutect.AO_f)+float(Mutect.RO_f)
	Mutect.DP_r=float(Mutect.AO_r)+float(Mutect.RO_r)

	try:
		Mutect.QB=(float((tumor[format.index('QSS')]).split(',')[0]) + float((tumor[format.index('QSS')]).split(',')[1]))/Mutect.DP_t
	except:
		Mutect.QB='.'
	try:
		Mutect.QB_n=(float((normal[format.index('QSS')]).split(',')[0]) + float((normal[format.index('QSS')]).split(',')[1]))/Mutect.DP_n
	except:
		Mutect.QB_n='.'
	
	Mutect.Call=1

	
	try:
		Mutect.SS= - Mutect.tumor_lod/Mutect.normal_lod
	except:
		Mutect.SS='.'
	
	try:
		Mutect.AF_t=Mutect.AO_t/Mutect.DP_t
	except:
		Mutect.AF_t=float(0)

	Mutect.delta_lod = Mutect.tumor_lod - Mutect.normal_lod
	
	try:
		Mutect.t_n_lod = Mutect.tumor_lod / Mutect.normal_lod
	except:
		Mutect.t_n_lod = '.'
	
	
	if opts.amplicon:
		if min(Mutect.DP_r,Mutect.DP_f)/(Mutect.DP_r+Mutect.DP_f) >= 0.05:
			Mutect.StrandBias=1-stats.fisher_exact([[Mutect.RO_f, Mutect.RO_r], [Mutect.AO_f, Mutect.AO_r]])[1]
		else:
			Mutect.StrandBias='.'
	else:
		if min(Mutect.DP_r,Mutect.DP_f)/(Mutect.DP_r+Mutect.DP_f) > 0:
			Mutect.StrandBias=1-stats.fisher_exact([[Mutect.RO_f, Mutect.RO_r], [Mutect.AO_f, Mutect.AO_r]])[1]

		else:
			Mutect.StrandBias='.'

	Mutect.FOXOG=tumor[format.index('FOXOG')]
	
	try:
		Mutect.GQ=float(tumor[format.index('GQ')])
	except:
		Mutect.GQ='.'
	try:
		Mutect.PGT=float(tumor[format.index('PGT')])
	except:
		Mutect.PGT='.'
	try:
		Mutect.PID=float(tumor[format.index('PID')])
	except:
		Mutect.PID='.'
	try:
		Mutect.PL=float(tumor[format.index('PL')])
	except:
		Mutect.PL='.'
	
	Mutect.FILTER=filter

	
def get_info_varscan(chrom,pos,ref,alt,filter,info,format,tumor,normal,varscan):
	'''estrae le informazioni dal vcf di varscan'''
	varscan.GT=tumor[format.index('GT')]
	varscan.AO_t=float(tumor[format.index('AD')])
	varscan.RO_t=float(tumor[format.index('RD')])
	varscan.DP=float(tumor[format.index('DP')])
	varscan.DP_t=float(tumor[format.index('DP')])

	if normal is not 'null' and normal[format.index('DP')] is not '0':
		varscan.AO_n=float(normal[format.index('AD')])
		varscan.RO_n=float(normal[format.index('RD')])
		varscan.DP_n=float(normal[format.index('DP')])
		varscan.AF_n=float(varscan.AO_n/varscan.DP_n)
		
	else:
		varscan.AO_n='.'
		varscan.RO_n='.'
		varscan.AF_n='.'
		varscan.DP_n=0

 	for ind in info:
 		if ind.startswith("MULTIALLEL"):
			varscan.MULTIALLELE=1
		if ind.startswith("SS=3"):
			varscan.LOH=1
		else:
			varscan.LOH=0 
		if ind.startswith('SOMATIC'):
			varscan.Somatic=1
			break
		else:
			varscan.Somatic=0
				
	varscan.RO_f=float(tumor[format.index('DP4')].split(',')[0])
	varscan.RO_r=float(tumor[format.index('DP4')].split(',')[1])
	varscan.AO_f=float(tumor[format.index('DP4')].split(',')[2])
	varscan.AO_r=float(tumor[format.index('DP4')].split(',')[3])
	varscan.DP_f=float(tumor[format.index('DP4')].split(',')[2])+float(tumor[format.index('DP4')].split(',')[0])
	varscan.DP_r=float(tumor[format.index('DP4')].split(',')[3])+float(tumor[format.index('DP4')].split(',')[1])
	varscan.Call=1
	varscan.QB='.'
	varscan.QB_n='.'
	try:
		varscan.AF_t=float(varscan.AO_t/(varscan.DP_t))
	except:
		varscan.AF_t=float(0)	

	if opts.amplicon:
		if min((varscan.DP_r),(varscan.DP_f))/(varscan.DP_r+varscan.DP_f) >= 0.05:
			varscan.StrandBias = 1-stats.fisher_exact([[varscan.RO_f, varscan.RO_r], [varscan.AO_f, varscan.AO_r]])[1]
		else:
		 	varscan.StrandBias='.'
	else:
		if min((varscan.DP_r),(varscan.DP_f))/(varscan.DP_r+varscan.DP_f) > 0:
			varscan.StrandBias = 1-stats.fisher_exact([[varscan.RO_f, varscan.RO_r], [varscan.AO_f, varscan.AO_r]])[1]
		else:
			varscan.StrandBias='.'

def get_info_vardict(chrom,pos,ref,alt,filter,info,format,tumor,normal,vardict):
	'''estrae le informazioni dal vcf di vardict'''
	vardict.GT=tumor[format.index('GT')]
	vardict.AO_t=float(tumor[format.index('VD')])
	vardict.RO_t=float(tumor[format.index('RD')].split(',')[0])+float(tumor[format.index('RD')].split(',')[1])
	vardict.DP_t=float(tumor[format.index('DP')])
	vardict.DP=float(tumor[format.index('DP')])

	if normal is not 'null' and normal[format.index('DP')] is not '0' :
		
		vardict.AO_n=float(normal[format.index('VD')])
		vardict.RO_n=float(normal[format.index('RD')].split(',')[0])+float(normal[format.index('RD')].split(',')[1])
		vardict.DP_n=float(normal[format.index('DP')])
		vardict.AF_n=float(vardict.AO_n/(vardict.DP_n))
	else:
		vardict.AO_n='.'
		vardict.RO_n='.'
		vardict.DP_n=0
		vardict.AF_n='.'
	
 	for ind in info:
 		
 		if ind.startswith("SHIFT3"):
 			vardict.SHIFT3=ind.split('=')[1]
 		if ind.startswith("MSI"):	
			vardict.MSI=ind.split('=')[1]
		if ind.startswith("MSILEN"):	
			vardict.MSILEN=ind.split('=')[1]
		if ind.startswith("SOR"):	
			vardict.SOR=ind.split('=')[1]
		if ind.startswith("LSEQ"):	
			vardict.LSEQ=ind.split('=')[1]
		if ind.startswith("RSEQ"):	
			vardict.RSEQ=ind.split('=')[1]
		if ind.startswith("TYPE"):	
			vardict.TYPE=ind.split('=')[1]
		if ind.startswith("STATUS"):	
			vardict.STATUS=ind.split('=')[1]
		if "LOH" in vardict.STATUS:
			vardict.LOH=1
		else:
			vardict.LOH=0
		
		if "Somatic" in vardict.STATUS:
			vardict.Somatic=1
		else:
			vardict.Somatic=0	

	vardict.RO_f=float(tumor[format.index('RD')].split(',')[0])
	vardict.RO_r=float(tumor[format.index('RD')].split(',')[1])
	vardict.AO_f=float(tumor[format.index('ALD')].split(',')[0])
	vardict.AO_r=float(tumor[format.index('ALD')].split(',')[1])
	vardict.DP_f=float(tumor[format.index('ALD')].split(',')[0]) + float(tumor[format.index('RD')].split(',')[0])
	vardict.DP_r=float(tumor[format.index('ALD')].split(',')[1]) + float(tumor[format.index('RD')].split(',')[1])
	vardict.QB=float(tumor[format.index('QUAL')])
	vardict.QB_n=float(normal[format.index('QUAL')])

	vardict.Call=1

	try:
		vardict.AF_t=float(vardict.AO_t/(vardict.AO_t + vardict.RO_t))
	except:
		vardict.AF_t=float(0)
	
	if opts.amplicon:
		if min(vardict.DP_r,vardict.DP_f)/(vardict.DP_r+vardict.DP_f) >= 0.05:
			vardict.StrandBias = 1-stats.fisher_exact([[vardict.RO_f, vardict.RO_r], [vardict.AO_f, vardict.AO_r]])[1]
		else:
			vardict.StrandBias='.'
	else:
		if min(vardict.DP_r,vardict.DP_f)/(vardict.DP_r+vardict.DP_f) > 0:
			vardict.StrandBias = 1-stats.fisher_exact([[vardict.RO_f, vardict.RO_r], [vardict.AO_f, vardict.AO_r]])[1]
		else:
			vardict.StrandBias='.'

	vardict.t_PMEAN=float(tumor[format.index('PMEAN')])
	vardict.t_PSTD=float(tumor[format.index('PSTD')])
	vardict.t_QSTD=float(tumor[format.index('QSTD')])
	vardict.t_SBF=float(tumor[format.index('SBF')])
	vardict.t_MQ=float(tumor[format.index('MQ')])
	vardict.t_SN=float(tumor[format.index('SN')])
	vardict.t_HIAF=float(tumor[format.index('HIAF')])
	vardict.t_NM=float(tumor[format.index('NM')])


	vardict.ODDRATIO=tumor[format.index('ODDRATIO')]
	vardict.SBF=tumor[format.index('SBF')]
	vardict.STRBIAS=tumor[format.index('BIAS')]

def set_features_snp(dictionary):
	'''setta i valori delle features in base alle info estratte dai vcf'''
	for variante in dictionary.keys():
		features=Features()
		varc_array=dictionary.get(variante)
		
		vett_MBQ=[]
		vett_MBQ_n=[]
		vett_DP=[]
		vett_DP_n=[]
		vett_AO_tum=[]
		vett_AO_norm=[]
		vett_RO_tum=[]
		vett_RO_norm=[]


		index=0

		for varcall in varc_array:
			if varcall is not "":
				vett_MBQ=vett_MBQ+[varcall.QB]
				vett_MBQ_n=vett_MBQ_n+[varcall.QB_n]
				vett_DP=vett_DP+[varcall.DP]
				vett_DP_n=vett_DP_n+[varcall.DP_n]
				vett_AO_tum=vett_AO_tum+[varcall.AO_t]
				vett_RO_tum=vett_RO_tum+[varcall.RO_t]
				vett_AO_norm=vett_AO_norm+[varcall.AO_n]
				vett_RO_norm=vett_RO_norm+[varcall.RO_n]
				
				if index == 0:

					features.GT_Mutect=varc_array[0].GT
					features.QB_Mutect=varc_array[0].QB
					features.DP_n_t_Mutect=float(varc_array[0].DP_n)/float(varc_array[0].DP_t)
					features.AF_t_Mutect=varc_array[0].AF_t
					features.AF_n_Mutect=varc_array[0].AF_n
					
					if features.AF_n_Mutect is not ".":
						features.Delta_Mutect=float(features.AF_t_Mutect) - float(features.AF_n_Mutect)

					try:
						features.Delta_perc_Mutect=float(features.AF_t_Mutect) / float(features.AF_n_Mutect)
					except:
						if features.AF_n_Mutect == 0:
							features.Delta_perc_Mutect='100'
						elif features.AF_n_Mutect is '.':
							features.Delta_perc_Mutect='.'
						
					features.CallMutect=varc_array[0].Call
					features.t_n_lod_Mutect=varc_array[0].t_n_lod
					features.tumor_lod_Mutect=varc_array[0].tumor_lod
					features.normal_lod_Mutect=varc_array[0].normal_lod
					features.delta_lod_Mutect=varc_array[0].delta_lod
					#features.SomaticMutect=varc_array[0].Somatic
					features.STRBIAS_Mutect=varc_array[0].StrandBias

					features.FOXOG_Mutect=varc_array[0].FOXOG
					features.GQ_Mutect=varc_array[0].GQ
					features.PGT_Mutect=varc_array[0].PGT
					features.PID_Mutect=varc_array[0].PID
					features.PL_Mutect=varc_array[0].PL

					features.HCNT_Mutect=varc_array[0].HCNT
					features.MAX_ED_Mutect=varc_array[0].MAX_ED
					features.MIN_ED_Mutect=varc_array[0].MIN_ED
					features.PON_Mutect=varc_array[0].PON
					features.RPA_Mutect=varc_array[0].RPA
					features.RU_Mutect=varc_array[0].RU
					features.STR_Mutect=varc_array[0].STR

					features.FILTER_Mutect=varc_array[0].FILTER


				elif index == 1:
					
					features.GT_Varscan=varc_array[1].GT
					features.DP_n_t_Varscan=float(varc_array[1].DP_n)/float(varc_array[1].DP_t)
					features.AF_t_Varscan=varc_array[1].AF_t
					features.AF_n_Varscan=varc_array[1].AF_n

					if features.AF_n_Varscan is not '.':
							features.Delta_Varscan=float(features.AF_t_Varscan) - float(features.AF_n_Varscan)
									
					try:
						features.Delta_perc_Varscan=float(features.AF_t_Varscan) / float(features.AF_n_Varscan)
					except:
						if features.AF_n_Varscan == 0.0:
							features.Delta_perc_Varscan='100'
						elif features.AF_n_Varscan is '.':
							features.Delta_perc_Varscan='.'

					features.SomaticVarscan=varc_array[1].Somatic
					features.CallVarscan=varc_array[1].Call
					features.SS_Varscan=varc_array[1].SS
					features.STRBIAS_Varscan=varc_array[1].StrandBias
					features.LOH_Varscan=varc_array[1].LOH
						
				elif index == 2:

					features.GT_Vardict=varc_array[2].GT
					features.QB_Vardict=varc_array[2].QB
					features.DP_n_t_Vardict=float(varc_array[2].DP_n)/float(varc_array[2].DP_t)
					features.AF_t_Vardict=varc_array[2].AF_t
					features.AF_n_Vardict=varc_array[2].AF_n

					if features.AF_n_Vardict is not ".":
						features.Delta_Vardict=float(features.AF_t_Vardict) - float(features.AF_n_Vardict)

					try:
						features.Delta_perc_Vardict=float(features.AF_t_Vardict) / float(features.AF_n_Vardict)
					except:
						if features.AF_n_Vardict == 0.0:
							features.Delta_perc_Vardict='100'
						elif features.AF_n_Vardict is '.':
							features.Delta_perc_Vardict='.'

					features.SomaticVardict=varc_array[2].Somatic
					features.CallVardict=varc_array[2].Call
					features.SS_Vardict=varc_array[2].SS
					features.STRBIAS_Vardict=varc_array[2].StrandBias
					features.BIAS_Vardict=varc_array[2].STRBIAS
					features.ODDRATIO_Vardict=varc_array[2].ODDRATIO
					features.SBF_Vardict=varc_array[2].SBF
					features.LOH_Vardict=varc_array[2].LOH

					features.SHIFT3_Vardict=varc_array[2].SHIFT3
					features.MSI_Vardict=varc_array[2].MSI
					features.MSILEN_Vardict=varc_array[2].MSILEN
					features.SOR_Vardict=varc_array[2].SOR
					features.LSEQ_Vardict=varc_array[2].LSEQ
					features.RSEQ_Vardict=varc_array[2].RSEQ
					features.STATUS_Vardict=varc_array[2].STATUS
					features.PMEAN_Vardict=varc_array[2].t_PMEAN
					features.PSTD_Vardict=varc_array[2].t_PSTD
					features.QSTD_Vardict=varc_array[2].t_QSTD
					features.SBF_Vardict=varc_array[2].t_SBF
					features.MQ_Vardict=varc_array[2].t_MQ
					features.SN_Vardict=varc_array[2].t_SN
					features.HIAF_Vardict=varc_array[2].t_HIAF
					features.NM_Vardict=varc_array[2].t_NM

			index = index + 1	

		vett_delta=[features.Delta_Mutect,features.Delta_Varscan,features.Delta_Vardict]
		vett_delta_perc=[features.Delta_perc_Mutect,features.Delta_perc_Varscan,features.Delta_perc_Vardict]
		vett_AF_media=[features.AF_t_Mutect,features.AF_t_Varscan,features.AF_t_Vardict]
		vett_AF_norm_media=[features.AF_n_Mutect,features.AF_n_Varscan,features.AF_n_Vardict]
		vett_STRB_media=[features.STRBIAS_Mutect,features.STRBIAS_Varscan,features.STRBIAS_Vardict]
		vett_DPn_t=[features.DP_n_t_Mutect,features.DP_n_t_Varscan,features.DP_n_t_Vardict]
		delta_m=0
		delta_m_perc=0
		AF_med=0
		SB_media=0
		nDP=0
		nMBQT=0
		DP_n_t_m=0
		AO_tum_media='.'
		RO_tum_media='.'
		AO_tum_median='.'
		RO_tum_median='.'

		AO_norm_media='.'
		RO_norm_media='.'
		AO_norm_median='.'
		RO_norm_median='.'

		i=0
		v=[]
		for dp in vett_DP:
			if dp and dp is not '':
				nDP= float(nDP)+float(dp)
				v=v+[float(dp)]
				i=i+1
		try:
			features.DP=nDP/i
		except:
			features.DP='.'
		try:
			features.DP_median= int(statistics.median(v))
		except:
			features.DP_median='.'

		i=0
		nDP=0
		v=[]
		for dp in vett_DP_n:
			if dp and dp is not '':
				nDP= float(nDP)+float(dp)
				v=v+[float(dp)]
				i=i+1
		try:
			features.DP_n=nDP/i
		except:
			features.DP_n='.'
		try:
			features.DP_norm_median= int(statistics.median(v))
		except:
			features.DP_norm_median='.'
		
		i=0
		v=[]
		nAO=0
		for ao in vett_AO_tum:
			if ao is not '':
				nAO= int(nAO)+int(ao)
				v=v+[int(ao)]
				i=i+1
		try:
			features.AO_tum_media=nAO/i
		except:
			features.AO_tum_media='.'
		try:
			features.AO_tum_median= int(statistics.median(v))
		except:
			features.AO_tum_median='.'

		i=0
		v=[]
		nAO=0		
		for ao in vett_AO_norm:
			if ao is not '.':
				nAO= int(nAO)+int(ao)
				v=v+[float(ao)]
				i=i+1

		try:
			features.AO_norm_media=nAO/i
			#print 'media',features.AO_norm_media
		except:
			features.AO_norm_media='.'
		try:
			features.AO_norm_median= int(statistics.median(v))
		except:
			features.AO_norm_median='.'
		i=0
		v=[]
		nRO=0
		for ro in vett_RO_tum:
			if ro is not '.':
				nRO= int(nRO)+int(ro)
				v=v+[int(ro)]
				i=i+1
		try:
			features.RO_tum_media=nRO/i
		except:
			features.RO_tum_media='.'
		try:
			features.RO_tum_median= statistics.median(v)
		except:
			features.RO_tum_median='.'

		i=0
		v=[]
		nRO=0
		for ro in vett_RO_norm:
			if ro is not '.':
				nRO= int(nRO)+int(ro)
				v=v+[int(ro)]
				i=i+1
		try:
			features.RO_norm_media=nRO/i
		except:
			features.RO_norm_media='.'
		try:
			features.RO_norm_median= statistics.median(v)
		except:
			features.RO_norm_median='.'

		i=0
		nMBQT=0
		v=[]
		for bq in vett_MBQ:
			#print variante,vett_MBQ,bq
			if bq and bq is not '.':
				nMBQT= float(nMBQT)+float(bq)
				if bq is not '0':
					v=v+[float(bq)]
				i=i+1
		try:
			features.MBQT=nMBQT/i
		except:
			features.MBQT='.'
		try:
			features.MBQT_median=round(statistics.median(v))
		except:
			features.MBQT_median='.'	

		i=0
		nMBQN=0
		v=[]
		for bq in vett_MBQ_n:
			#print variante,vett_MBQ,bq
			if bq is not '.':
				try:
					nMBQN= float(nMBQN)+float(bq)
				except:
					print variante,vett_MBQ_n,bq
				if bq is not '0':
					v=v+[float(bq)]
				i=i+1
		try:
			features.MBQN_media=nMBQN/i
		except:
			features.MBQN_media='.'
		try:
			features.MBQN_median=round(statistics.median(v))
		except:
			features.MBQN_median='.'

		i=0
		v=[]
		for delta in vett_delta:
			if delta is not '.':
				delta_m = float(delta_m)+ float(delta)
				v=v+[float(delta)]
				i=i+1
		try:
			features.delta_media = delta_m/i
		except:
			features.delta_media='.'
		try:
			features.delta_median=statistics.median(v)
		except:
			features.delta_median='.'
		
		i=0
		v=[]
		for delta_perc in vett_delta_perc:
			if delta_perc and delta_perc is not '.':
				delta_m_perc = float(delta_m_perc) + float(delta_perc)
				if delta_perc is not '0':
					v=v+[float(delta_perc)]
				i=i+1
		try:
			features.Delta_perc_media= delta_m_perc/i
		except:
			features.Delta_perc_media='.'
		try:
			features.Delta_perc_median= statistics.median(v)
		except:
			features.Delta_perc_median='.'

		i=0
		v=[]
		for af in vett_AF_media:
			if af and af is not '0':
				AF_med=float(AF_med) + float(af)
				v=v+[float(af)]
				i=i+1
		try:
			features.AF_media= AF_med/i
		except:
			features.AF_media='.'
		try:
			features.AF_median= statistics.median(v)
		except:
			features.AF_median='.'

		i=0
		AF_med=0
		v=[]
		for af in vett_AF_norm_media:
			if af is not '.':
				AF_med=float(AF_med) + float(af)
				v=v+[float(af)]
				i=i+1
		try:
			features.AF_norm_media= AF_med/i
		except:
			features.AF_norm_media='.'
		try:
			features.AF_norm_median= statistics.median(v)
		except:
			features.AF_norm_median='.'

		i=0
		v=[]
		#print vett_STRB_media
		for strb in vett_STRB_media:
			if strb is not '.':
				SB_media=float(SB_media) + float(strb)
				v=v+[float(strb)]
				i=i+1
		#print v,statistics.median(v)
		try:
			features.STRBIAS_media= SB_media/i
		except:
			features.STRBIAS_media='.'
		try:
			features.STRBIAS_median = statistics.median(v)
		except:
			features.STRBIAS_median='.'

		i=0
		v=[]
		for dp_n_t in vett_DPn_t:
			if dp_n_t is not '.':
				DP_n_t_m= float(DP_n_t_m)+float(dp_n_t)
				v=v+[float(dp_n_t)]
				i=i+1
		try:
			features.DPn_t_media=DP_n_t_m/i
		except:
			features.DPn_t_media='.'
		try:
			features.DPn_t_median= statistics.median(v)
		except:
			features.DPn_t_median='.'

		#print vett_delta,"\tdeltamedia",features.delta_media,"\n",vett_delta_perc,"\tdelta_perc_media",features.Delta_perc_media,"\n",vett_AF_media,"\taf_media",features.AF_media,"\n",vett_STRB_media,"\tstrb",features.STRBIAS_media
		dictionary[variante]= varc_array + [features]

def switch_snp(dictionary,ID,index,chrom,pos,ref,alt,filter,info,format,tumor,normal):
	'''tramite index richiama la funzione di estrazione delle informazioni del variant caller associato all'indice'''
	if dictionary.has_key(ID):
		vettore=dictionary[ID]
	else:
		vettore=['','','']

	if index==0:
		# print 'Mutect'
		mutect=Mutect()
		get_info_Mutect(chrom,pos,ref,alt,filter,info,format,tumor,normal,mutect)
		if Mutect.AF_t != 0.0 or Mutect.AF_t != '.': 
			vettore[0]=mutect
	elif index==3:
		# print 'vardict'
		vardict=Vardict()
		get_info_vardict(chrom,pos,ref,alt,filter,info,format,tumor,normal,vardict)
		if vardict.AF_t != 0.0 or vardict.AF_t != '.': 
			vettore[2]=vardict
	elif index==1 or index==2:
		# print 'varscan'
		varscan=Varscan()
		get_info_varscan(chrom,pos,ref,alt,filter,info,format,tumor,normal,varscan)
		if varscan.AF_t != 0.0 or varscan.AF_t != '.': 
			vettore[1]=varscan
	dictionary[ID]=vettore

def read(iterable,index,dictionary):
	'''legge il vcf e splitta le varie sezioni'''
	normal_tumor=0
	chrom=''
	alt=''
	pos=''
	ref=''
	filter=''
	format=''
	info=''
	for line in iterable:
		line.rstrip()
		#print line
		if line.startswith('#CHROM'):
			parts = line.split("\t")
			if 'NORMAL' in parts[9] or parts[9] == opts.normal:
				#print "prima normal",parts[9]
				normal_tumor=1
			else: 
				#print "prima tumor"
				normal_tumor=0
		elif line.startswith('#'):
			continue
		else:
			#print line
			ind=0
			parts = line.split("\t")
			#print parts[1]
			chrom=parts[0]
			pos=parts[1]
			ref=parts[3]
			alt=parts[4]
			filter=parts[6]
			ID='\t'.join([chrom,pos,ref,alt])
			INFO=parts[7]
			info=INFO.split(";") 
			FORMAT=parts[8]
			format=FORMAT.split(":")
			tumor=''
			normal=''
			if normal_tumor==1:
				NORMAL=parts[9]
				TUMOR=parts[10].rstrip()
				tumor=TUMOR.split(":")
				normal=NORMAL.split(":")
			else:
				NORMAL=parts[10].rstrip()
				TUMOR=parts[9]
				tumor=TUMOR.split(":")
				normal=NORMAL.split(":")
			
			if NORMAL.startswith('.:.:.:') or len(normal)==1:
				normal='null'
			if TUMOR.startswith('.:.:.:') or len(tumor)==1 :
				continue
			else:
				if "VD" in format:
					if tumor[format.index('VD')] is '0':
						continue
					else:
						switch_snp(dictionary,ID,index,chrom,pos,ref,alt,filter,info,format,tumor,normal)
				else:	
					switch_snp(dictionary,ID,index,chrom,pos,ref,alt,filter,info,format,tumor,normal)

def max_delta_perc(dictionary):
	''' calcolo il max delta percentuale'''
	buffer_delta_perc='0'
	max_delta_perc='0'
	
	for variante in dictionary.keys():
		features = dictionary.get(variante)[5]
		if features.Delta_perc_Mutect == '.':
			mutect=0
		else:
			mutect=features.Delta_perc_Mutect
		if features.Delta_perc_Vardict == '.':
			vard=0
		else:
			vard=features.Delta_perc_Vardict
		if features.Delta_perc_Varscan == '.':
			varsc=0
		else:
			varsc=features.Delta_perc_Varscan
		if features.Delta_perc_Radia == '.':
			rad=0
		else:
			rad=features.Delta_perc_Radia
		if features.Delta_perc_Somaticsniper == '.':
			som=0
		else:
			som=features.Delta_perc_Somaticsniper
		buffer_delta_perc= max([float(mutect),float(varsc),float(vard),float(rad),float(som)])
		if float(max_delta_perc) < float(buffer_delta_perc):
			max_delta_perc=buffer_delta_perc
	#print max_delta_perc		

def control(dictionary):
	''' esegue un controllo sulle varianti, se non hanno variant caller che le chiama vengono eliminate'''
	for variante in dictionary.keys():
		if dictionary[variante][:3] == ['','','']:
			#print "sto cancellando:",variante
			del dictionary[variante]
				
def print_var_snp_complete(dictionary):
	varianti_tsv=open(opts.out+ '.tsv','w')
	varianti_tsv.write('\t'.join(["SAMPLE_NORMAL_ID","SAMPLE_TUMOR_ID","CHROM","POS","REF","ALT","CallMutect","CallVarscan","CallVardict",
			"SomaticMutect","SomaticVarscan","SomaticVardict",
			"GT_Mutect","GT_Varscan","GT_Vardict","DP_median","BQ_Mutect","BQ_Vardict","MBQT_media","MBQT_mediana",
			"AF_Mutect","AF_Varscan","AF_Vardict","AF_media","AF_median",
			"Delta_Mutect","Delta_Varscan","Delta_Vardict","Delta_media","Delta_median",
			"Delta_perc_Mutect","Delta_perc_Varscan","Delta_perc_Vardict","Delta_perc_media","Delta_perc_median",
			"DP_n/t_Mutect","DP_n/t_varscan","DP_n/t_vardict","DP_n/t_media","DP_n/t_median",
			"T_lod_Mutect","N_lod_Mutect","Delta_lod_Mutect","T/N_lod_Mutect",
			"FOXOG_Mutect","GQ_Mutect","PGT_Mutect","PID_Mutect","PL_Mutect","HCNT_Mutect",
			"MAX_ED_Mutect","MIN_ED_Mutect","PON_Mutect","RPA_Mutect","RU_Mutect","STR_Mutect",
			"STRBIAS_Mutect","STRBIAS_Varscan","STRBIAS_Vardict","STRBIAS_medio","STRBIAS_median",
			"ODDRATIO_Vardict","SBF_Vardict","SHIFT3_Vardict","MSI_Vardict","MSILEN_Vardict","SOR_Vardict",
			"LSEQ_Vardict","RSEQ_Vardict","STATUS_Vardict","PMEAN_Vardict","PSTD_Vardict","QSTD_Vardict",
			"MQ_Vardict","SN_Vardict","HIAF_Vardict","NM_Vardict",
			"LOH_Varscan","LOH_Vardict"]) +'\n')


	for variante in dictionary.keys():
		features = dictionary.get(variante)[-1]

		varianti_tsv.write('\t'.join([opts.normal,opts.tumor,variante,str(features.CallMutect),str(features.CallVarscan),str(features.CallVardict),
			str(features.SomaticMutect),str(features.SomaticVarscan),str(features.SomaticVardict),
			str(features.GT_Mutect),str(features.GT_Varscan),str(features.GT_Vardict),
			str(features.DP_median),str(features.QB_Mutect),str(features.QB_Vardict),str(features.MBQT),str(features.MBQT_median),
			str(features.AF_t_Mutect),str(features.AF_t_Varscan),str(features.AF_t_Vardict),str(features.AF_media),str(features.AF_median),
			str(features.Delta_Mutect),str(features.Delta_Varscan),str(features.Delta_Vardict),str(features.delta_media),str(features.delta_median),
			str(features.Delta_perc_Mutect),str(features.Delta_perc_Varscan),str(features.Delta_perc_Vardict),str(features.Delta_perc_media),str(features.Delta_perc_median),
			str(features.DP_n_t_Mutect),str(features.DP_n_t_Varscan),str(features.DP_n_t_Vardict),str(features.DPn_t_media),str(features.DPn_t_median),
			str(features.tumor_lod_Mutect),str(features.normal_lod_Mutect),str(features.delta_lod_Mutect),str(features.t_n_lod_Mutect),
			str(features.FOXOG_Mutect),str(features.GQ_Mutect),str(features.PGT_Mutect),str(features.PID_Mutect),str(features.PL_Mutect),str(features.HCNT_Mutect),
			str(features.MAX_ED_Mutect),str(features.MIN_ED_Mutect),str(features.PON_Mutect),str(features.RPA_Mutect),str(features.RU_Mutect),str(features.STR_Mutect),
			str(features.STRBIAS_Mutect),str(features.STRBIAS_Varscan),str(features.STRBIAS_Vardict),str(features.STRBIAS_media),str(features.STRBIAS_median),
			str(features.ODDRATIO_Vardict),str(features.SBF_Vardict),str(features.SHIFT3_Vardict),str(features.MSI_Vardict),str(features.MSILEN_Vardict),str(features.SOR_Vardict),
			str(features.LSEQ_Vardict),str(features.RSEQ_Vardict),str(features.STATUS_Vardict),str(features.PMEAN_Vardict),str(features.PSTD_Vardict),str(features.QSTD_Vardict),
			str(features.MQ_Vardict),str(features.SN_Vardict),str(features.HIAF_Vardict),str(features.NM_Vardict),
			str(features.LOH_Varscan),str(features.LOH_Vardict)]).rstrip() +'\n')

def print_var_snp_reduced(dictionary):
	
	varianti_tsv=open(opts.out+ '.tsv','w')
	
	varianti_tsv.write('\t'.join(["CHROM","POS","ID","REF","ALT","CallMutect","CallVarscan","CallVardict",
			"SomaticMutect","SomaticVarscan","SomaticVardict",
			"FILTER_Mutect","STATUS_Vardict",
			"GT_Mutect","GT_Varscan","GT_Vardict",
			"DP_TUM","AF_TUM","AO_TUM","RO_TUM","MBQ_TUM",
			"DP_NORM","AF_NORM","AO_NORM","RO_NORM","MBQ_NORM",
			"Delta_mediana","Delta_perc_median",
			"STRBIAS_mediana","HCNT_Mutect","MAX_ED_Mutect","MIN_ED_Mutect","STR_Mutect",
			"SHIFT3_Vardict","MSI_Vardict","MSILEN_Vardict","SOR_Vardict",
			"LSEQ_Vardict","RSEQ_Vardict","PMEAN_Vardict",
			"MQ_Vardict","SN_Vardict","NM_Vardict",
			"LOH_Varscan","LOH_Vardict"]) + '\n')

	for variante in dictionary.keys():
		features = dictionary.get(variante)[-1]

		varianti_tsv.write('\t'.join([variante.split('\t')[0],variante.split('\t')[1],opts.tumor,variante.split('\t')[2],variante.split('\t')[3],
			str(features.CallMutect),str(features.CallVarscan),str(features.CallVardict),
			str(features.SomaticMutect),str(features.SomaticVarscan),str(features.SomaticVardict),
			str(features.FILTER_Mutect),str(features.STATUS_Vardict),
			str(features.GT_Mutect),str(features.GT_Varscan),str(features.GT_Vardict),
			str(features.DP_median),str(features.AF_median),str(features.AO_tum_media),str(features.RO_tum_media),str(features.MBQT_median),
			str(features.DP_norm_median),str(features.AF_norm_median),str(features.AO_norm_media),str(features.RO_norm_media),str(features.MBQN_median),
			str(features.delta_median),str(features.Delta_perc_median),
			str(features.STRBIAS_median),str(features.HCNT_Mutect),str(features.MAX_ED_Mutect),str(features.MIN_ED_Mutect),str(features.STR_Mutect),
			str(features.SHIFT3_Vardict),str(features.MSI_Vardict),str(features.MSILEN_Vardict),str(features.SOR_Vardict),
			str(features.LSEQ_Vardict),str(features.RSEQ_Vardict),str(features.PMEAN_Vardict),
			str(features.MQ_Vardict),str(features.SN_Vardict),str(features.NM_Vardict),
			str(features.LOH_Varscan),str(features.LOH_Vardict)]).rstrip() + '\n')
	
	varianti_tsv.close()

def print_vcf(varianti):
	varianti_vcf=open(opts.out+ '.vcf','w')
	varianti_vcf.write('##fileformat=VCFv4.2\n'+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLES\n')
	for variante in varianti.keys():
		var_vcf=variante.split('\t')[0]+'\t'+variante.split('\t')[1]+'\t.\t'+variante.split('\t')[2]+'\t'+variante.split('\t')[3]+'\t.\t.\t.\t.\t.'
		varianti_vcf.write(var_vcf+ '\n')
	varianti_vcf.close()

def print_var(dictionary):

	lista_features=open(opts.listaFeatures,'r')
	varianti_tsv=open(opts.out+ '.tsv','w')
	
	header=[]
	features_variante=[]
	
	for line in lista_features:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		else:
			header=header+[line]
			features_variante=features_variante+['features.'+line]

	varianti_tsv.write('CHROM\tPOS\tID\tREF\tALT\t' + '\t'.join(header)+ '\n')
	
	for variante in dictionary.keys():
		#print variante
		features = dictionary.get(variante)[-1]
		features_variante_eval=[]
		for feat in features_variante:
			feat_eval=str(eval(feat))
			features_variante_eval=features_variante_eval + [feat_eval]
		var=variante.split('\t')[0]+'\t'+variante.split('\t')[1]+'\t' +opts.tumor +'\t'+variante.split('\t')[2]+'\t'+variante.split('\t')[3]+ '\t' + '\t'.join(features_variante_eval)
		if features.AF_media != '.':
			varianti_tsv.write(var+ '\n')
	varianti_tsv.close()
	lista_features.close()
		
def main():


	
	parser = argparse.ArgumentParser('Parse VCF output from Variant callers to output a variant_dataset.txt.  Output is to stdout.')
	parser.add_argument('-m', '--mutect', help="Mutect vcf output file name")
	parser.add_argument('-d', '--vardict', help="Vardict vcf output file name")
	parser.add_argument('-v', '--varscan_snp', help="Varscan vcf output file name")
	parser.add_argument('-i', '--varscan_indel', help="Varscan vcf output file name")
	parser.add_argument('-n','--normal',help="Name of normal sample")
	parser.add_argument('-t','--tumor',help="Name of tumor sample")
	parser.add_argument('-a','--amplicon',help="Amplicon design", action='store_true')
	parser.add_argument('-c','--complete',help="Print complete info", action='store_true')
	parser.add_argument('-o', '--out', help="file name in output. It returns file_name.features.tsv and file_name.vcf ")
	parser.add_argument('-l', '--listaFeatures', help="Lista di features da stampare",default=None)
	
	global opts 
	opts = parser.parse_args()
	
	callers = [opts.mutect,opts.varscan_snp,opts.varscan_indel,opts.vardict]
	varianti = dict() 
	index=0;
	for vcf in callers:
		#print vcf
		in_file = open(vcf,'r')
		vcfreader = read(in_file,index,varianti)
		index = index + 1
	set_features_snp(varianti)
	control(varianti)
	
	if opts.listaFeatures:
		print_var(varianti)
	else:
		if opts.complete:
			print_var_snp_complete(varianti)
		else:
			print_var_snp_reduced(varianti)

	print_vcf(varianti)
main()
