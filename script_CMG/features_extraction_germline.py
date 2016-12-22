import re
import string
import argparse
import sys
import statistics
import os
import scipy.stats as stats


''' //////////// CLASSI ////////////'''

class Caller():
	GT=''
	AO=''
	RO=''
	AO_f=''
	AO_r=''
	DP_f=''
	DP_r=''
	DP=''
	QB=''
	GQ=''
	Call='0'
	AF=''
	STRBIAS=''
	FILTER=''
	AC=''
	AN=''
	STRBIAS_TOT=''

class Freebayes(Caller):
	
	#PROVENIENTI DALLE INFO
	AO_f_TOT=''
	AO_r_TOT=''
	RO_f_TOT=''
	RO_r_TOT=''
	DP_f_TOT=''
	DP_r_TOT=''
	AB=''
	ABP=''
	AN=''
	AF_TOT=''
	AO_TOT=''
	CIGAR=''
	DP_TOT=''
	DPB_TOT=''
	DPRA_TOT=''
	END_TOT=''
	EPP_TOT=''
	EPPR_TOT=''
	GTI_TOT=''
	LEN=''
	MEANALT=''
	MIN=''
	MQ_A=''
	MQ_R=''
	NS=''
	NUMALT=''
	ODDS=''
	PAIRED=''
	PAIREDR=''
	PAO=''
	PQA=''
	PQR=''
	PRO=''
	QA=''
	QR=''
	RO_TOT=''
	RPL=''
	RPP=''
	RPPR=''
	RPR=''
	RUN=''
	SAP=''
	SRP=''

	#PROVENIENTI DAL FORMAT
	lod=''
	A_QB=''
	R_QB=''

class GATK(Caller):
	
	#PROVENIENTI DALLE INFO
	AF_TOT=''
	BaseQRankSum=''
	ClippingRankSum=''
	DP_TOT=''
	DS=''
	END=''
	ExcessHet=''
	FS=''
	HaplotypeScore=''
	InbreedingCoeff=''
	MLEAC=''
	MLEAF=''
	MQ=''
	MQRankSum=''
	QD=''
	RAW_MQ=''
	ReadPosRankSum=''
	SOR=''

class Varscan(Caller):

	SDP=''
	
	ADP=''
	WT=''
	HET=''
	HOM=''
	NC=''

class Features():

	GT_Varscan='.'
	GT_Freebayes='.'
	GT_GATK='.'

	AO_media='.'
	RO_media='.'
	AO_mediana='.'
	RO_mediana='.'

	DP_media='.'
	DP_mediana='.'

	MQB_GATK='.'
	MQB_Varscan='.'
	MQB_Freebayes='.'
	MQB_media='.'
	MQB_mediana='.'
	
	AF_GATK='0'
	AF_Varscan='0'
	AF_Freebayes='0'
	AF_media='.'
	AF_mediana='.'
	
	CallGATK='0'
	CallVarscan='0'
	CallFreebayes='0'
	
	STRBIAS_GATK='.'
	STRBIAS_Varscan='.'
	STRBIAS_Freebayes='.'
	STRBIAS_media='.'
	STRBIAS_mediana= '.'

	STRBIAS_TOT_GATK='.'
	STRBIAS_TOT_Varscan='.'
	STRBIAS_TOT_Freebayes='.'
	STRBIAS_TOT_media='.'
	STRBIAS_TOT_mediana= '.'
	
	FILTER_GATK='.'
	FILTER_Varscan='.'
	FILTER_Freebayes='.'
	
	AC_GATK='.'
	AC_Varscan='.'
	AC_Freebayes='.'
	AC_media=''
	AC_mediana=''

	AN_GATK='.'
	AN_Varscan='.'
	AN_Freebayes='.'


	#features solo di GATK

	AF_TOT_GATK=''
	BaseQRankSum_GATK=''
	ClippingRankSum_GATK=''
	DP_TOT_GATK=''
	DS_GATK=''
	END_GATK=''
	ExcessHet_GATK=''
	FS_GATK=''
	HaplotypeScore_GATK=''
	InbreedingCoeff_GATK=''
	MLEAC_GATK=''
	MLEAF_GATK=''
	MQ_GATK=''
	MQRankSum_GATK=''
	QD_GATK=''
	RAW_MQ_GATK=''
	ReadPosRankSum_GATK=''
	SOR_GATK=''

	#features solo di Freebayes
	lod_Freebayes=''
	AO_f_TOT_Freebayes=''
	AO_r_TOT_Freebayes=''
	RO_f_TOT_Freebayes=''
	RO_r_TOT_Freebayes=''
	DP_f_TOT_Freebayes=''
	DP_r_TOT_Freebayes=''
	AB_Freebayes=''
	ABP_Freebayes=''
	AF_TOT_Freebayes=''
	AO_TOT_Freebayes=''
	CIGAR_Freebayes=''
	DP_TOT_Freebayes=''
	DPB_TOT_Freebayes=''
	DPRA_TOT_Freebayes=''
	END_TOT_Freebayes=''
	EPP_TOT_Freebayes=''
	EPPR_TOT_Freebayes=''
	GTI_TOT_Freebayes=''
	LEN_Freebayes=''
	MEANALT_Freebayes=''
	MIN_Freebayes=''
	MQ_A_Freebayes=''
	MQ_R_Freebayes=''
	NS_Freebayes=''
	NUMALT_Freebayes=''
	ODDS_Freebayes=''
	PAIRED_Freebayes=''
	PAIREDR_Freebayes=''
	PAO_Freebayes=''
	PQA_Freebayes=''
	PQR_Freebayes=''
	PRO_Freebayes=''
	QA_Freebayes=''
	QR_Freebayes=''
	RO_TOT_Freebayes=''
	RPL_Freebayes=''
	RPP_Freebayes=''
	RPPR_Freebayes=''
	RPR_Freebayes=''
	RUN_Freebayes=''
	SAP_Freebayes=''
	SRP_Freebayes=''

	A_QB_Freebayes=''
	R_QB_Freebayes=''

	#features solo di Varscan
	SDP_Varscan=''


''' //////////// FUNZIONI ////////////'''

def get_info_Freebayes(chrom,pos,ref,alt,filter,info,format,sample,freebayes):
	'''estrae le informazioni dal vcf di freebayes'''
	
	freebayes.GT=sample[format.index('GT')]
	if freebayes.GT=='.' :
		freebayes.GT='./.'
		#print chrom,pos,ref,alt,freebayes.GT
	else:
		freebayes.GQ=sample[format.index('GQ')]
		freebayes.AO=float(sample[format.index('AO')])
		freebayes.RO=float(sample[format.index('RO')])
		freebayes.DP=float(sample[format.index('DP')])

		for ind in info:
			if ind.startswith("SAF="):
				freebayes.AO_f_TOT=float(ind.split('=')[1])
			if ind.startswith("SAR="):
				freebayes.AO_r_TOT=float(ind.split('=')[1])
			if ind.startswith("SRF="):
				freebayes.RO_f_TOT=float(ind.split('=')[1])
			if ind.startswith("SRR="):
				freebayes.RO_r_TOT=float(ind.split('=')[1])
			if ind.startswith("AC="):
				freebayes.AC=float(ind.split('=')[1])
			if ind.startswith("AB="):
				freebayes.AB=float(ind.split('=')[1])
			if ind.startswith("ABP="):
				freebayes.ABP=float(ind.split('=')[1])
			if ind.startswith("AN="):
				freebayes.AN=float(ind.split('=')[1])
			if ind.startswith("AF="):
				freebayes.AF_TOT=float(ind.split('=')[1])
			if ind.startswith("AO="):
				freebayes.AO_TOT=float(ind.split('=')[1])
			if ind.startswith("DP="):
				freebayes.DP_TOT=float(ind.split('=')[1])
			if ind.startswith("DPB="):
				freebayes.DPB_TOT=float(ind.split('=')[1])
			if ind.startswith("DPRA="):
				freebayes.DPRA_TOT=float(ind.split('=')[1])
			if ind.startswith("END="):
				freebayes.END_TOT=float(ind.split('=')[1])
			if ind.startswith("EPP="):
				freebayes.EPP_TOT=float(ind.split('=')[1])
			if ind.startswith("EPPR="):
				freebayes.EPPR_TOT=float(ind.split('=')[1])
			if ind.startswith("GTI="):
				freebayes.GTI_TOT=float(ind.split('=')[1])
			if ind.startswith("LEN="):
				freebayes.LEN=float(ind.split('=')[1])
			if ind.startswith("MEANALT="):
				freebayes.MEANALT=float(ind.split('=')[1])
			if ind.startswith("MIN="):
				freebayes.MIN=float(ind.split('=')[1])
			if ind.startswith("MQM="):
				freebayes.MQ_A=float(ind.split('=')[1])
			if ind.startswith("MQMR="):
				freebayes.MQ_R=float(ind.split('=')[1])
			if ind.startswith("NS="):
				freebayes.NS=float(ind.split('=')[1])
			if ind.startswith("NUMALT="):
				freebayes.NUMALT=float(ind.split('=')[1])
			if ind.startswith("ODDS="):
				freebayes.ODDS=float(ind.split('=')[1])
			if ind.startswith("PAIRED="):
				freebayes.PAIRED=float(ind.split('=')[1])
			if ind.startswith("PAIREDR="):
				freebayes.PAIREDR=float(ind.split('=')[1])
			if ind.startswith("PAO="):
				freebayes.PAO=float(ind.split('=')[1])
			if ind.startswith("PQA="):
				freebayes.PQA=float(ind.split('=')[1])
			if ind.startswith("PQR="):
				freebayes.PQR=float(ind.split('=')[1])
			if ind.startswith("PRO="):
				freebayes.PRO=float(ind.split('=')[1])
			if ind.startswith("QA="):
				freebayes.QA=float(ind.split('=')[1])
			if ind.startswith("QR="):
				freebayes.QR=float(ind.split('=')[1])
			if ind.startswith("RO="):
				freebayes.RO_TOT=float(ind.split('=')[1])
			if ind.startswith("RPL="):
				freebayes.RPL=float(ind.split('=')[1])
			if ind.startswith("RPP="):
				freebayes.RPP=float(ind.split('=')[1])
			if ind.startswith("RPPR="):
				freebayes.RPPR=float(ind.split('=')[1])
			if ind.startswith("RPR="):
				freebayes.RPR=float(ind.split('=')[1])
			if ind.startswith("RUN="):
				freebayes.RUN=float(ind.split('=')[1])
			if ind.startswith("SAP="):
				freebayes.SAP=float(ind.split('=')[1])
			if ind.startswith("SRP="):
				freebayes.SRP=float(ind.split('=')[1])


		freebayes.DP_f_TOT=float(freebayes.AO_f_TOT)+float(freebayes.RO_f_TOT)
		freebayes.DP_r_TOT=float(freebayes.AO_r_TOT)+float(freebayes.RO_r_TOT)

		try:
			freebayes.A_QB=float(sample[format.index('QA')])/freebayes.AO
		except:
			freebayes.A_QB='.'
		try:	
			freebayes.R_QB=float(sample[format.index('QR')])/freebayes.RO
		except:
			freebayes.R_QB='.'
		try:	
			freebayes.QB=(float(sample[format.index('QA')]) + float(sample[format.index('QR')]))/freebayes.DP
		except:
			if sample[format.index('QA')] == '0':
				freebayes.QB=freebayes.R_QB
			elif sample[format.index('QR')] == '0':
				freebayes.QB=freebayes.R_QA
		
	


		try:
			gls = [float(x) for x in sample[format.index("GL")].split(",")]
			freebayes.lod = max(gls[i] - gls[0] for i in range(1, len(gls)))
		except:
			freebayes.lod = -1.0
		
		freebayes.AF=freebayes.AO/freebayes.DP

		if opts.amplicon:
			if min(freebayes.DP_r_TOT,freebayes.DP_f_TOT)/(freebayes.DP_r_TOT+freebayes.DP_f_TOT) >= 0.05:
				freebayes.STRBIAS_TOT=1-stats.fisher_exact([[freebayes.RO_f_TOT, freebayes.RO_r_TOT], [freebayes.AO_f_TOT, freebayes.AO_r_TOT]])[1]
			else:
				freebayes.STRBIAS_TOT='.'
		else:
			if min(freebayes.DP_r_TOT,freebayes.DP_f_TOT)/(freebayes.DP_r_TOT+freebayes.DP_f_TOT) > 0:
				freebayes.STRBIAS_TOT=1-stats.fisher_exact([[freebayes.RO_f_TOT, freebayes.RO_r_TOT], [freebayes.AO_f_TOT, freebayes.AO_r_TOT]])[1]

			else:
				freebayes.STRBIAS_TOT='.'

		freebayes.FILTER=filter
		#print chrom,pos,ref,alt,freebayes.GT
	
	freebayes.Call=1

def get_info_GATK(chrom,pos,ref,alt,filter,info,format,sample,GATK):
	'''estrae le informazioni dal vcf di GATK'''
	
	GATK.GT=sample[format.index('GT')]
	if GATK.GT=='./.':
		#print chrom,pos,ref,alt,GATK.GT
		pass
	else:
		GATK.AO=float((sample[format.index('AD')]).split(',')[1])
		GATK.RO=float((sample[format.index('AD')]).split(',')[0])
		GATK.DP=GATK.AO+GATK.RO	
		
		GATK.STR='0'
		for ind in info:
			if ind.startswith("AC="):
				GATK.AC=ind.split('=')[1]
			if ind.startswith("AF="):
				GATK.AF_TOT=ind.split('=')[1]
			if ind.startswith("AN="):
				GATK.AN=ind.split('=')[1]
			if ind.startswith("BaseQRankSum="):
				GATK.BaseQRankSum=ind.split('=')[1]
			if ind.startswith("ClippingRankSum="):
				GATK.ClippingRankSum=ind.split('=')[1]
			if ind.startswith("DP="):
				GATK.DP_TOT=ind.split('=')[1]
			if ind.startswith("DS="):
				GATK.DS=ind.split('=')[1]
			if ind.startswith("END="):
				GATK.END=ind.split('=')[1]
			if ind.startswith("ExcessHet="):
				GATK.ExcessHet=ind.split('=')[1]
			if ind.startswith("FS="):
				GATK.FS=ind.split('=')[1]
			if ind.startswith("HaplotypeScore="):
				GATK.HaplotypeScore=ind.split('=')[1]
			if ind.startswith("InbreedingCoeff="):
				GATK.InbreedingCoeff=ind.split('=')[1]
			if ind.startswith("MLEAC="):
				GATK.MLEAC=ind.split('=')[1]
			if ind.startswith("MLEAF="):
				GATK.MLEAF=ind.split('=')[1]
			if ind.startswith("MQ="):
				GATK.MQ=ind.split('=')[1]
			if ind.startswith("MQRankSum="):
				GATK.MQRankSum=ind.split('=')[1]
			if ind.startswith("QD="):
				GATK.QD=ind.split('=')[1]
			if ind.startswith("RAW_MQ="):
				GATK.RAW_MQ=ind.split('=')[1]
			if ind.startswith("ReadPosRankSum="):
				GATK.ReadPosRankSum=ind.split('=')[1]
			if ind.startswith("SOR="):
				GATK.SOR=ind.split('=')[1]

	
		
		try:
			GATK.GQ=float(sample[format.index('GQ')])
		except:
			GATK.GQ='.'
		
		GATK.FILTER=filter
		#print chrom,pos,ref,alt,GATK.GT
	GATK.Call=1
	
def get_info_varscan(chrom,pos,ref,alt,filter,info,format,sample,varscan):
	'''estrae le informazioni dal vcf di varscan'''
	
	varscan.GT=sample[format.index('GT')]
	if varscan.GT== './.':
		#print chrom,pos,ref,alt,varscan.GT
		pass
	else:
		varscan.AO=float(sample[format.index('AD')])
		varscan.RO=float(sample[format.index('RD')])
		varscan.DP=varscan.RO + varscan.AO
	
	 	for ind in info:
	 		if ind.startswith("ADP"):
				varscan.ADP=ind.split('=')[1]
			if ind.startswith("WT"):
				varscan.WT=ind.split('=')[1]
			if ind.startswith("HET"):
				varscan.HET=ind.split('=')[1]
			if ind.startswith("HOM"):
				varscan.HOM=ind.split('=')[1]
			if ind.startswith("NC"):
				varscan.NC=ind.split('=')[1]
		
		varscan.AC=	(varscan.HET+2*varscan.HOM)
		varscan.AN=	2*(varscan.HET+varscan.HOM+varscan.WT)
					
		varscan.RO_f=float(sample[format.index('RDF')])
		varscan.RO_r=float(sample[format.index('RDR')])
		varscan.AO_f=float(sample[format.index('ADF')])
		varscan.AO_r=float(sample[format.index('ADR')])
		varscan.DP_f=varscan.RO_f + varscan.AO_f
		varscan.DP_r=varscan.RO_r + varscan.AO_r
		varscan.Call=1
		varscan.QB_R=float(sample[format.index('RBQ')])
		varscan.QB_A=float(sample[format.index('ABQ')])
		Varscan.GQ=float(sample[format.index('GQ')])
		
		varscan.AF=float(varscan.AO/(varscan.DP))
		

		if opts.amplicon:
			if min((varscan.DP_r),(varscan.DP_f))/(varscan.DP_r+varscan.DP_f) >= 0.05:
				varscan.STRBIAS = 1-stats.fisher_exact([[varscan.RO_f, varscan.RO_r], [varscan.AO_f, varscan.AO_r]])[1]
			else:
			 	varscan.STRBIAS='.'
		else:
			if min((varscan.DP_r),(varscan.DP_f))/(varscan.DP_r+varscan.DP_f) > 0:
				varscan.STRBIAS = 1-stats.fisher_exact([[varscan.RO_f, varscan.RO_r], [varscan.AO_f, varscan.AO_r]])[1]
			else:
				varscan.STRBIAS='.'

		Varscan.FILTER=filter
		#print chrom,pos,ref,alt,varscan.GT
	Varscan.Call=1

def set_features(dictionary):
	'''setta i valori delle features in base alle info estratte dai vcf'''
	for variante in dictionary.keys():
		features=Features()
		varc_array=dictionary.get(variante)
		
		vett_MBQ=[]
		vett_DP=[]
		vett_AO=[]
		vett_RO=[]
		vett_AC=[]
		index=0
		#print variante
		#print varc_array
		for varcall in varc_array:
			if varcall is not "":
				vett_MBQ=vett_MBQ+[varcall.QB]
				vett_DP=vett_DP+[varcall.DP]
				vett_AO=vett_AO+[varcall.AO]
				vett_RO=vett_RO+[varcall.RO]
				vett_AC=vett_AC+[varcall.AC]
				
				if index == 0:

					features.lod_Freebayes=varc_array[0].lod			
					features.GT_Freebayes=varc_array[0].GT
					features.AO_Freebayes=varc_array[0].AO
					features.RO_Freebayes=varc_array[0].RO
					features.AO_f_Freebayes=varc_array[0].AO_f
					features.AO_r_Freebayes=varc_array[0].AO_r
					features.DP_f_Freebayes=varc_array[0].DP_f
					features.DP_r_Freebayes=varc_array[0].DP_r
					features.DP_Freebayes=varc_array[0].DP
					features.QB_Freebayes=varc_array[0].QB
					features.GQ_Freebayes=varc_array[0].GQ
					features.CallFreebayes=varc_array[0].Call
					features.AF_Freebayes=varc_array[0].AF
					features.STRBIAS_Freebayes=varc_array[0].STRBIAS
					features.FILTER_Freebayes=varc_array[0].FILTER
					features.AC_Freebayes=varc_array[0].AC
					features.AN_Freebayes=varc_array[0].AN	
					features.STRBIAS_TOT_Freebayes=varc_array[0].STRBIAS_TOT

					features.AO_f_TOT_Freebayes=varc_array[0].AO_f_TOT
					features.AO_r_TOT_Freebayes=varc_array[0].AO_r_TOT
					features.RO_f_TOT_Freebayes=varc_array[0].RO_f_TOT
					features.RO_r_TOT_Freebayes=varc_array[0].RO_r_TOT
					features.DP_f_TOT_Freebayes=varc_array[0].DP_f_TOT
					features.DP_r_TOT_Freebayes=varc_array[0].DP_r_TOT
					features.AB_Freebayes=varc_array[0].AB
					features.ABP_Freebayes=varc_array[0].ABP
					features.AF_TOT_Freebayes=varc_array[0].AF_TOT
					features.AO_TOT_Freebayes=varc_array[0].AO_TOT
					features.CIGAR_Freebayes=varc_array[0].CIGAR
					features.DP_TOT_Freebayes=varc_array[0].DP_TOT
					features.DPB_TOT_Freebayes=varc_array[0].DPB_TOT
					features.DPRA_TOT_Freebayes=varc_array[0].DPRA_TOT
					features.END_TOT_Freebayes=varc_array[0].END_TOT
					features.EPP_TOT_Freebayes=varc_array[0].EPP_TOT
					features.EPPR_TOT_Freebayes=varc_array[0].EPPR_TOT
					features.GTI_TOT_Freebayes=varc_array[0].GTI_TOT
					features.LEN_Freebayes=varc_array[0].LEN
					features.MEANALT_Freebayes=varc_array[0].MEANALT
					features.MIN_Freebayes=varc_array[0].MIN
					features.MQ_A_Freebayes=varc_array[0].MQ_A
					features.MQ_R_Freebayes=varc_array[0].MQ_R
					features.NS_Freebayes=varc_array[0].NS
					features.NUMALT_Freebayes=varc_array[0].NUMALT
					features.ODDS_Freebayes=varc_array[0].ODDS
					features.PAIRED_Freebayes=varc_array[0].PAIRED
					features.PAIREDR_Freebayes=varc_array[0].PAIREDR
					features.PAO_Freebayes=varc_array[0].PAO
					features.PQA_Freebayes=varc_array[0].PQA
					features.PQR_Freebayes=varc_array[0].PQR
					features.PRO_Freebayes=varc_array[0].PRO
					features.QA_Freebayes=varc_array[0].QA
					features.QR_Freebayes=varc_array[0].QR
					features.RO_TOT_Freebayes=varc_array[0].RO_TOT
					features.RPL_Freebayes=varc_array[0].RPL
					features.RPP_Freebayes=varc_array[0].RPP
					features.RPPR_Freebayes=varc_array[0].RPPR
					features.RPR_Freebayes=varc_array[0].RPR
					features.RUN_Freebayes=varc_array[0].RUN
					features.SAP_Freebayes=varc_array[0].SAP
					features.SRP_Freebayes=varc_array[0].SRP
					features.A_QB_Freebayes=varc_array[0].A_QB
					features.R_QB_Freebayes=varc_array[0].R_QB

				if index == 2:

					features.GT_GATK=varc_array[2].GT
					features.AO_GATK=varc_array[2].AO
					features.RO_GATK=varc_array[2].RO
					features.AO_f_GATK=varc_array[2].AO_f
					features.AO_r_GATK=varc_array[2].AO_r
					features.DP_f_GATK=varc_array[2].DP_f
					features.DP_r_GATK=varc_array[2].DP_r
					features.DP_GATK=varc_array[2].DP
					features.QB_GATK=varc_array[2].QB
					features.GQ_GATK=varc_array[2].GQ
					features.CallGATK=varc_array[2].Call
					features.AF_GATK=varc_array[2].AF
					features.STRBIAS_GATK=varc_array[2].STRBIAS
					features.FILTER_GATK=varc_array[2].FILTER
					features.AC_GATK=varc_array[2].AC
					features.AN_GATK=varc_array[2].AN
					features.STRBIAS_TOT_GATK=varc_array[2].STRBIAS_TOT

					features.AF_TOT_GATK=varc_array[2].AF_TOT
					features.BaseQRankSum_GATK=varc_array[2].BaseQRankSum
					features.ClippingRankSum_GATK=varc_array[2].ClippingRankSum
					features.DP_TOT_GATK=varc_array[2].DP_TOT
					features.DS_GATK=varc_array[2].DS
					features.END_GATK=varc_array[2].END
					features.ExcessHet_GATK=varc_array[2].ExcessHet
					features.FS_GATK=varc_array[2].FS
					features.HaplotypeScore_GATK=varc_array[2].HaplotypeScore
					features.InbreedingCoeff_GATK=varc_array[2].InbreedingCoeff
					features.MLEAC_GATK=varc_array[2].MLEAC
					features.MLEAF_GATK=varc_array[2].MLEAF
					features.MQ_GATK=varc_array[2].MQ
					features.MQRankSum_GATK=varc_array[2].MQRankSum
					features.QD_GATK=varc_array[2].QD
					features.RAW_MQ_GATK=varc_array[2].RAW_MQ
					features.ReadPosRankSum_GATK=varc_array[2].ReadPosRankSum
					features.SOR_GATK=varc_array[2].SOR


				elif index == 1:
					
					features.GT_Varscan=varc_array[1].GT
					features.AO_Varscan=varc_array[1].AO
					features.RO_Varscan=varc_array[1].RO
					features.AO_f_Varscan=varc_array[1].AO_f
					features.AO_r_Varscan=varc_array[1].AO_r
					features.DP_f_Varscan=varc_array[1].DP_f
					features.DP_r_Varscan=varc_array[1].DP_r
					features.DP_Varscan=varc_array[1].DP
					features.QB_Varscan=varc_array[1].QB
					features.GQ_Varscan=varc_array[1].GQ
					features.CallVarscan=varc_array[1].Call
					features.AF_Varscan=varc_array[1].AF
					features.STRBIAS_Varscan=varc_array[1].STRBIAS
					features.FILTER_Varscan=varc_array[1].FILTER
					features.AC_Varscan=varc_array[1].AC
					features.AN_Varscan=varc_array[1].AN
					features.STRBIAS_TOT_Varscan=varc_array[1].STRBIAS_TOT
					features.SDP_Varscan=varc_array[1].SDP


			index = index + 1	

		
		vett_AF_media=[features.AF_GATK,features.AF_Varscan,features.AF_Freebayes]
		vett_STRB_media=[features.STRBIAS_GATK,features.STRBIAS_Varscan,features.STRBIAS_Freebayes]
		
		AF_media=0
		SB_media=0
		AC_media=0
		AC_mediana=0
		DP_media='.'
		DP_mediana='.'
		AO_media='.'
		RO_media='.'
		AO_mediana='.'
		RO_mediana='.'

		v=[]
		for dp in vett_DP:
			if dp and dp is not '':
				v=v+[float(dp)]
		try:
			features.DP_media=statistics.mean(v)
		except:
			features.DP_media='.'
		try:
			features.DP_mediana=statistics.median(v)
		except:
			features.DP_mediana='.'

		v=[]
		for ao in vett_AO:
			if ao and ao is not '':
				v=v+[int(ao)]
		try:
			features.AO_media=statistics.mean(v)
		except:
			features.AO_media='.'
		try:
			features.AO_mediana=statistics.median(v)
		except:
			features.AO_mediana='.'

		v=[]
		for ro in vett_RO:
			if ro and ro is not '.' :
				v=v+[int(ro)]
		try:
			features.RO_media=statistics.mean(v)
		except:
			features.RO_media='.'
		try:
			features.RO_mediana=statistics.median(v)
		except:
			features.RO_mediana='.'

		v=[]
		for bq in vett_MBQ:
			if bq and bq is not '.':
				v=v+[float(bq)]
		try:
			features.MBQ_media=statistics.mean(v)
		except:
			features.MBQ_media='.'
		try:
			features.MBQ_mediana=statistics.median(v)
		except:
			features.MBQ_mediana='.'	

		v=[]
		for af in vett_AF_media:
			if af and af is not '0':
				v=v+[float(af)]
		try:
			features.AF_media= statistics.mean(v)
		except:
			features.AF_media='.'
		try:
			features.AF_mediana= statistics.median(v)
		except:
			features.AF_mediana='.'

		v=[]
		for strb in vett_STRB_media:
			if strb and strb is not '.':
				v=v+[float(strb)]
		try:
			features.STRBIAS_media= statistics.mean(v)
		except:
			features.STRBIAS_media='.'
		try:
			features.STRBIAS_median = statistics.median(v)
		except:
			features.STRBIAS_mediana='.'

		v=[]
		for ac in vett_AC:
			if ac and ac is not '.':
				v=v+[int(ac)]
		try:
			features.AC_media=statistics.mean(v)
		except:
			features.AC_media='.'
		try:
			features.AC_mediana= statistics.median(v)
		except:
			features.AC_mediana='.'

		dictionary[variante]= varc_array + [features]

def switch(dictionary,ID,index,chrom,pos,ref,alt,filter,info,format,sample):
	'''tramite index richiama la funzione di estrazione delle informazioni del variant caller associato all'indice'''
	if dictionary.has_key(ID):
		vettore=dictionary[ID]
	else:
		vettore=['','','']

	if index==0:
		# print 'Freebayes'
		freebayes=Freebayes()
		get_info_Freebayes(chrom,pos,ref,alt,filter,info,format,sample,freebayes)
		vettore[0]=freebayes

	elif index==2:
		# print 'gatk'
		gatk=GATK()
		get_info_GATK(chrom,pos,ref,alt,filter,info,format,sample,gatk)
		vettore[2]=gatk
	elif index==1 :
		# print 'varscan'
		varscan=Varscan()
		get_info_varscan(chrom,pos,ref,alt,filter,info,format,sample,varscan)
		vettore[1]=varscan
	dictionary[ID]=vettore

def read(iterable,index,dictionary):
	'''legge il vcf e splitta le varie sezioni'''
	for line in iterable:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		else:
			ind=0
			var = line.split("\t")
			#print var[1]
			chrom=var[0]
			pos=var[1]
			ref=var[3]
			alt=var[4]
			filter=var[6]
			ID='\t'.join([chrom,pos,ref,alt])
			INFO=var[7]
			info=INFO.split(";") 
			FORMAT=var[8]
			format=FORMAT.split(":")
			SAMPLE = var[-1]
			sample = SAMPLE.split(':')
			if alt != '*':
				#print index,chrom,pos,ref,alt
				switch(dictionary,ID,index,chrom,pos,ref,alt,filter,info,format,sample)
				
def control(dictionary):
	''' esegue un controllo sulle varianti, se non hanno variant caller che le chiama vengono eliminate'''
	for variante in dictionary.keys():
		if dictionary[variante][:3] == ['','','']:
			#print "sto cancellando:",variante
			del dictionary[variante]
				
def print_var(dictionary,out,sample_name):

	lista_features=open(opts.listaFeatures,'r')
	dataset_varianti=open(out + '/' + sample_name + '.tsv','w')
	


	header=[]
	features_variante=[]
	
	
	for line in lista_features:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		else:
			header=header+[line]
			features_variante=features_variante+['features.'+line]

	#print features_variante

	dataset_varianti.write('CHROM\tPOS\tID\tREF\tALT\t' + '\t'.join(header)+ '\n')
	for variante in dictionary.keys():
		features = dictionary.get(variante)[-1]

		#print features_variante
		features_variante_eval=[]
		for feat in features_variante:
			#print feat
		 	feat_eval=str(eval(feat))
		 	features_variante_eval=features_variante_eval + [feat_eval]
		
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	20151202_01_Cardio
		var=variante.split('\t')[0]+'\t'+variante.split('\t')[1]+'\t' +sample_name +'\t'+variante.split('\t')[2]+'\t'+variante.split('\t')[3]+ '\t' + '\t'.join(features_variante_eval)
		
		for gt in [features.GT_GATK,features.GT_Varscan,features.GT_Freebayes]:
			if gt != './.' and gt != '.' and gt != '0/0':
				dataset_varianti.write(var+ '\n')
				break
	dataset_varianti.close()
	lista_features.close()
	
def print_vcf(varianti,out):
	dataset_varianti_vcf=open(out+ '/TOTAL.vcf','w')
	dataset_varianti_vcf.write('##fileformat=VCFv4.2\n'+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLES\n')
	for variante in dictionary.keys():
		var_vcf=variante.split('\t')[0]+'\t'+variante.split('\t')[1]+'\t.\t'+variante.split('\t')[2]+'\t'+variante.split('\t')[3]+'\t.\t.\t.\t.\t.'
		dataset_varianti_vcf.write(var_vcf+ '\n')
	dataset_varianti_vcf.close()


def samples_name_extract(vcf):
	samples = []
	for line in vcf:
		line=line.rstrip()
		if line.startswith('#CHROM'):
			line_split = line.split('\t')
			samples = line_split[9:]
	return samples
		
def split_vcf(vcf_dir,samples):

	vcf_name = os.path.basename(vcf_dir)
	vcf_path = os.path.dirname(vcf_dir)
	vcf = open(vcf_dir,'r')
	header = []
	header_chrom = []
	varianti = []
	
	if 'Free' in vcf_name:
			variant_caller = 'FreeBayes'
	elif 'GATK' in vcf_name:
			variant_caller = 'GATK'
	elif 'VarScan' in vcf_name:
			variant_caller = 'VarScan'

	for line in vcf:
		line=line.rstrip()
		if line.startswith('##'):
		 	header = header + [line]
		elif line.startswith('#CHROM'):
		 	header_chrom = line.split('\t')
		else:
		 	varianti = varianti + [line]
	i=0
	
	for sample in samples:
		try: 
			os.mkdir(vcf_path +'/' + sample)
		except:
			pass

		sample_vcf = open(vcf_path +'/' + sample +'/' + sample + '_'+variant_caller +'.vcf','w')
		print 'Riscrivo le varianti in ' +vcf_path +'/'+ sample + '_'+variant_caller +'.vcf'

		sample_vcf.write('\n'.join(header) + '\n')
		sample_vcf.write('\t'.join(header_chrom[0:9] + [sample])  +'\n')
		
		for variante in varianti:
			variante_split = variante.split('\t')
			variante_common = variante_split[0:9]
			format_sample = variante_split[8 + samples.index(sample) +1  ]
			sample_vcf.write('\t'.join(variante_common + [format_sample]) + '\n')
			#if format_sample.startswith('0/0:') or format_sample.startswith('./.'):
			#	continue
			#else:
			#	sample_vcf.write('\t'.join(variante_common + [format_sample]) + '\n')
		sample_vcf.close()

def main():

	
	parser = argparse.ArgumentParser('Parse VCF output from Variant callers to output a variant_dataset.txt.  Output is to stdout.')
	parser.add_argument('-f', '--freebayes', help="Freebayes vcf output file name")
	parser.add_argument('-g', '--gatk', help="gatk vcf output file name")
	parser.add_argument('-v', '--varscan', help="Varscan vcf output file name")
	parser.add_argument('-l', '--listaFeatures', help="Lista di features da stampare")
	parser.add_argument('-s', '--split', help="Split vcf per samples", action='store_true')
	parser.add_argument('-a','--amplicon',help="Amplicon design", action='store_true')

	global opts 
	opts = parser.parse_args()
	callers = [opts.gatk,opts.varscan,opts.freebayes]
	samples = samples_name_extract(open(opts.freebayes,'r'))
	
	if opts.split:
		print 'Splittando le varianti per campione...'
		for vcf_dir in callers:
			print 'Splitto ' + vcf_dir
			split_vcf(vcf_dir,samples)

	#vcf_path =  os.path.dirname('/home/minime/Scrivania/VCF_TEST/20151202_01_Cardio/20151202_01_Cardio_FreeBayes.vcf')

	#dir_sample = '20151202_01_Cardio'
	out=os.path.dirname(opts.freebayes)+'/out'
	try:
		os.mkdir(out)
	except:
		pass

	for dir_sample in os.listdir(os.path.dirname(opts.freebayes)):
		varianti = dict()
		vcf_path = os.path.dirname(opts.freebayes) +'/' + dir_sample
		if os.path.isdir(vcf_path) and vcf_path != out:
			for vcf_name in os.listdir(vcf_path) :
				print vcf_name
				if 'Free' in vcf_name:
					index = 0
				elif 'GATK' in vcf_name:
					index = 2
				elif 'VarScan' in vcf_name:
					index = 1
				
				in_file = open(vcf_path + '/' + vcf_name)
				vcfreader = read(in_file,index,varianti)
			
			set_features(varianti)
			print_var(varianti,out,dir_sample)
	
	print_vcf(varianti,out)
main()
