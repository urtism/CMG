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
	STRBIAS='.'
	FILTER=''
	AC=''
	AN=''
	STRBIAS_TOT=''
	RF=''

class iEVA():

	SimpleRepeat=''
	SimpleRepeatLength=''
	SimpleRepeatUnit=''
	PseudoNucleotidesComposition=''
	RepeatMasker=''
	gcContent=''
	VariantClass=''

	#bam extraction

	StrandBiasReads=''
	UnMappedReads=''
	MeanMappingQuality=''
	MappingQualityZero=''
	NotPrimaryAlignment=''
	SupplementaryAlignment=''
	NotPairedReads=''
	NotProperPairedReads=''
	AlignmentScore=''
	TotalDupReads=''
	
	#Genotype extraction

	NumberReadDupRef=''
	NumberReadDupAlt=''
	DuplicateReference=''
	DuplicateAlternate=''
	DeltaDuplicate=''
	iEvaDepth=''
	iAlleleDepthiEVA=''
	ReadRef=''
	ReadAlt=''
	MeanRefQscore=''
	MeanAltQscore=''
	RefMeanMappingQuality=''
	AltMeanMappingQuality=''
	TotalDPUnfilter=''
	NumberClippedReadsRef=''
	NumberClippedReadsAlt=''
	ClippedReadsRef=''
	ClippedReadsAlt=''


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
	#PROVENIENTI DAL FORMAT
	sQD=''
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

	RF_Varscan='.'
	RF_Freebayes='.'
	RF_GATK='.'

	QB_GATK='.'
	QB_Varscan='.'
	QB_Freebayes='.'
	MBQ_media='.'
	MBQ_mediana='.'
	
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
	AC_media='.'
	AC_mediana='.'
	AC_min='.'
	AC_max='.'
	AC='.'

	AN_GATK='.'
	AN_Varscan='.'
	AN_Freebayes='.'
	AN='.'

	#features solo di GATK
	sQD_GATK='.'
	AF_TOT_GATK='.'
	BaseQRankSum_GATK='.'
	ClippingRankSum_GATK='.'
	DP_TOT_GATK='.'
	DS_GATK='.'
	END_GATK='.'
	ExcessHet_GATK='.'
	FS_GATK='.'
	HaplotypeScore_GATK='.'
	InbreedingCoeff_GATK='.'
	MLEAC_GATK='.'
	MLEAF_GATK='.'
	MQ_GATK='.'
	MQRankSum_GATK='.'
	QD_GATK='.'
	RAW_MQ_GATK='.'
	ReadPosRankSum_GATK='.'
	SOR_GATK='.'

	#features solo di Freebayes
	lod_Freebayes='.'
	AO_f_TOT_Freebayes='.'
	AO_r_TOT_Freebayes='.'
	RO_f_TOT_Freebayes='.'
	RO_r_TOT_Freebayes='.'
	DP_f_TOT_Freebayes='.'
	DP_r_TOT_Freebayes='.'
	AB_Freebayes='.'
	ABP_Freebayes='.'
	AF_TOT_Freebayes='.'
	AO_TOT_Freebayes='.'
	CIGAR_Freebayes='.'
	DP_TOT_Freebayes='.'
	DPB_TOT_Freebayes='.'
	DPRA_TOT_Freebayes='.'
	END_TOT_Freebayes='.'
	EPP_TOT_Freebayes='.'
	EPPR_TOT_Freebayes='.'
	GTI_TOT_Freebayes='.'
	LEN_Freebayes='.'
	MEANALT_Freebayes='.'
	MIN_Freebayes='.'
	MQ_A_Freebayes='.'
	MQ_R_Freebayes='.'
	NS_Freebayes='.'
	NUMALT_Freebayes='.'
	ODDS_Freebayes='.'
	PAIRED_Freebayes='.'
	PAIREDR_Freebayes='.'
	PAO_Freebayes='.'
	PQA_Freebayes='.'
	PQR_Freebayes='.'
	PRO_Freebayes='.'
	QA_Freebayes='.'
	QR_Freebayes='.'
	RO_TOT_Freebayes='.'
	RPL_Freebayes='.'
	RPP_Freebayes='.'
	RPPR_Freebayes='.'
	RPR_Freebayes='.'
	RUN_Freebayes='.'
	SAP_Freebayes='.'
	SRP_Freebayes='.'

	A_QB_Freebayes='.'
	R_QB_Freebayes='.'

	#features solo di Varscan
	SDP_Varscan=''


	#features di iEVA

	#Reference extraction arguments
	SimpleRepeat_iEVA='.'
	SimpleRepeatLength_iEVA='.'
	SimpleRepeatUnit_iEVA='.'
	PseudoNucleotidesComposition_iEVA='.'
	RepeatMasker_iEVA='.'
	gcContent_iEVA='.'
	VariantClass_iEVA='.'

	#bam extraction

	StrandBiasReads_iEVA='.'
	UnMappedReads_iEVA='.'
	MeanMappingQuality_iEVA='.'
	MappingQualityZero_iEVA='.'
	NotPrimaryAlignment_iEVA='.'
	SupplementaryAlignment_iEVA='.'
	NotPairedReads_iEVA='.'
	NotProperPairedReads_iEVA='.'
	AlignmentScore_iEVA='.'
	TotalDupReads_iEVA='.'
	
	#Genotype extraction

	NumberReadDupRef_iEVA='.'
	NumberReadDupAlt_iEVA='.'
	DuplicateReference_iEVA='.'
	DuplicateAlternate_iEVA='.'
	DeltaDuplicate_iEVA='.'
	iEvaDepth_iEVA='.'
	iAlleleDepth_iEVA='.'
	ReadRef_iEVA='.'
	ReadAlt_iEVA='.'
	MeanRefQscore_iEVA='.'
	MeanAltQscore_iEVA='.'
	RefMeanMappingQuality_iEVA='.'
	AltMeanMappingQuality_iEVA='.'
	TotalDPUnfilter_iEVA='.'
	NumberClippedReadsRef_iEVA='.'
	NumberClippedReadsAlt_iEVA='.'
	ClippedReadsRef_iEVA='.'
	ClippedReadsAlt_iEVA='.'


''' //////////// FUNZIONI ////////////'''


def get_info_iEVA(chrom,pos,ref,alt,filter,info,format,sample,ieva):

	if "SBR" in format:
		ieva.StrandBiasReads=sample[format.index("SBR")]
	if "UnMap" in format:
		ieva.UnMappedReads=sample[format.index("UnMap")]
	if "MQ0" in format:
		ieva.MappingQualityZero=sample[format.index("MQ0")]
	if "MMQ" in format:
		ieva.MeanMappingQuality=sample[format.index("MMQ")]
	if "NPA" in format:
		ieva.NotPrimaryAlignment=sample[format.index("NPA")]
	if "SA" in format:
		ieva.SupplementaryAlignment=sample[format.index("SA")]
	if "NP" in format:
		ieva.NotPairedReads=sample[format.index("NP")]
	if "NPP" in format:
		ieva.NotProperPairedReads=sample[format.index("NPP")]
	if "AS" in format:
		ieva.AlignmentScore=sample[format.index("AS")]
	if "TDR" in format:
		ieva.TotalDupReads=sample[format.index("TDR")]

	if "iNDR" in format:
		ieva.NumberReadDupRef=sample[format.index("iNDR")]
	if "iNDA" in format:
		ieva.NumberReadDupAlt=sample[format.index("iNDA")]
	if "iDR" in format:
		ieva.DuplicateReference=sample[format.index("iDR")]
	if "iDA" in format:
		ieva.DuplicateAlternate=sample[format.index("iDA")]
	if "iDDup" in format:
		ieva.DeltaDuplicate=sample[format.index("iDDup")]
	if "iDP" in format:
		ieva.iEvaDepth=sample[format.index("iDP")]
	if "iAD" in format:
		ieva.iAlleleDepth=sample[format.index("iAD")]
	if "iRR" in format:
		ieva.ReadRef=sample[format.index("iRR")]
	if "iRA" in format:
		ieva.ReadAlt=sample[format.index("iRA")]
	if "iQR" in format:
		ieva.MeanRefQscore=sample[format.index("iQR")]
	if "iQA" in format:
		ieva.MeanAltQscore=sample[format.index("iQA")]

	if "TDP" in format:
		ieva.TotalDPUnfilter=sample[format.index("TDP")]
	if "iRMQ" in format:
		ieva.RefMeanMappingQuality=sample[format.index("iRMQ")]
	if "iAMQ" in format:
		ieva.AltMeanMappingQuality=sample[format.index("iAMQ")]
	if "iNCR" in format:
		ieva.NumberClippedReadsRef=sample[format.index("NCR")]
	if "iNCA" in format:
		ieva.NumberClippedReadsAlt=sample[format.index("NCA")]
	if "iCR" in format:
		ieva.ClippedReadsRef=sample[format.index("CR")]
	if "iCA" in format:
		ieva.ClippedReadsAlt=sample[format.index("CA")]

	for el in info:
		if el.startswith('SRL='):
			ieva.SimpleRepeatLength=el.split('=')[1]
		elif el.startswith('SRU='):
			ieva.SimpleRepeatUnit=el.split('=')[1]
		elif el.startswith('SR='):
			ieva.SimpleRepeat=el.split('=')[1]
		elif el.startswith('PNC='):
			ieva.PseudoNucleotidesComposition=el.split('=')[1]
		elif el.startswith('RM='):
			ieva.RepeatMasker=el.split('=')[1]
		elif el.startswith('GC='):
			ieva.gcContent=el.split('=')[1]
		elif el.startswith('VC='):
			ieva.VariantClass=el.split('=')[1]
		
def get_info_Freebayes(chrom,pos,ref,alt,filter,info,format,sample,freebayes,merged):
	'''estrae le informazioni dal vcf di freebayes'''
	if merged==1:
		suffisso='_F'
	else:
		suffisso=''

	freebayes.GT=sample[format.index('GT'+suffisso)]

	if freebayes.GT=='.':
		freebayes.GT='./.'
		#print chrom,pos,ref,alt,freebayes.GT
	if freebayes.GT=='./.':
		pass
	else:

		freebayes.GQ=sample[format.index('GQ'+suffisso)]
		freebayes.AO=float(sample[format.index('AO'+suffisso)])
		freebayes.RO=float(sample[format.index('RO'+suffisso)])
		#freebayes.DP=freebayes.AO+freebayes.RO
		freebayes.DP=float(sample[format.index('DP'+suffisso)])

		try:
			freebayes.RF=(freebayes.AO+freebayes.RO)/freebayes.DP
		except:
			freebayes.RF='.'

		for ind in info:
			if ind.startswith("SAF"+suffisso+"="):
				freebayes.AO_f_TOT=float(ind.split('=')[1])
			if ind.startswith("SAR"+suffisso+"="):
				freebayes.AO_r_TOT=float(ind.split('=')[1])
			if ind.startswith("SRF"+suffisso+"="):
				freebayes.RO_f_TOT=float(ind.split('=')[1])
			if ind.startswith("SRR"+suffisso+"="):
				freebayes.RO_r_TOT=float(ind.split('=')[1])
			if ind.startswith("AC"+suffisso+"="):
				freebayes.AC=int(ind.split('=')[1])
			if ind.startswith("AB"+suffisso+"="):
				freebayes.AB=float(ind.split('=')[1])
			if ind.startswith("ABP"+suffisso+"="):
				freebayes.ABP=float(ind.split('=')[1])
			if ind.startswith("AN"+suffisso+"="):
				freebayes.AN=float(ind.split('=')[1])
			if ind.startswith("AF"+suffisso+"="):
				freebayes.AF_TOT=float(ind.split('=')[1])
			if ind.startswith("AO"+suffisso+"="):
				freebayes.AO_TOT=float(ind.split('=')[1])
			if ind.startswith("DP"+suffisso+"="):
				freebayes.DP_TOT=float(ind.split('=')[1])
			if ind.startswith("DPB"+suffisso+"="):
				freebayes.DPB_TOT=float(ind.split('=')[1])
			if ind.startswith("DPRA"+suffisso+"="):
				freebayes.DPRA_TOT=float(ind.split('=')[1])
			if ind.startswith("END"+suffisso+"="):
				freebayes.END_TOT=float(ind.split('=')[1])
			if ind.startswith("EPP"+suffisso+"="):
				freebayes.EPP_TOT=float(ind.split('=')[1])
			if ind.startswith("EPPR"+suffisso+"="):
				freebayes.EPPR_TOT=float(ind.split('=')[1])
			if ind.startswith("GTI"+suffisso+"="):
				freebayes.GTI_TOT=float(ind.split('=')[1])
			if ind.startswith("LEN"+suffisso+"="):
				freebayes.LEN=float(ind.split('=')[1])
			if ind.startswith("MEANALT"+suffisso+"="):
				freebayes.MEANALT=float(ind.split('=')[1])
			if ind.startswith("MIN"+suffisso+"="):
				freebayes.MIN=float(ind.split('=')[1])
			if ind.startswith("MQM"+suffisso+"="):
				freebayes.MQ_A=float(ind.split('=')[1])
			if ind.startswith("MQMR"+suffisso+"="):
				freebayes.MQ_R=float(ind.split('=')[1])
			if ind.startswith("NS"+suffisso+"="):
				freebayes.NS=float(ind.split('=')[1])
			if ind.startswith("NUMALT"+suffisso+"="):
				freebayes.NUMALT=float(ind.split('=')[1])
			if ind.startswith("ODDS"+suffisso+"="):
				freebayes.ODDS=float(ind.split('=')[1])
			if ind.startswith("PAIRED"+suffisso+"="):
				freebayes.PAIRED=float(ind.split('=')[1])
			if ind.startswith("PAIREDR"+suffisso+"="):
				freebayes.PAIREDR=float(ind.split('=')[1])
			if ind.startswith("PAO"+suffisso+"="):
				freebayes.PAO=float(ind.split('=')[1])
			if ind.startswith("PQA"+suffisso+"="):
				freebayes.PQA=float(ind.split('=')[1])
			if ind.startswith("PQR"+suffisso+"="):
				freebayes.PQR=float(ind.split('=')[1])
			if ind.startswith("PRO"+suffisso+"="):
				freebayes.PRO=float(ind.split('=')[1])
			if ind.startswith("QA"+suffisso+"="):
				freebayes.QA=float(ind.split('=')[1])
			if ind.startswith("QR"+suffisso+"="):
				freebayes.QR=float(ind.split('=')[1])
			if ind.startswith("RO"+suffisso+"="):
				freebayes.RO_TOT=float(ind.split('=')[1])
			if ind.startswith("RPL"+suffisso+"="):
				freebayes.RPL=float(ind.split('=')[1])
			if ind.startswith("RPP"+suffisso+"="):
				freebayes.RPP=float(ind.split('=')[1])
			if ind.startswith("RPPR"+suffisso+"="):
				freebayes.RPPR=float(ind.split('=')[1])
			if ind.startswith("RPR"+suffisso+"="):
				freebayes.RPR=float(ind.split('=')[1])
			if ind.startswith("RUN"+suffisso+"="):
				freebayes.RUN=float(ind.split('=')[1])
			if ind.startswith("SAP"+suffisso+"="):
				freebayes.SAP=float(ind.split('=')[1])
			if ind.startswith("SRP"+suffisso+"="):
				freebayes.SRP=float(ind.split('=')[1])


		freebayes.DP_f_TOT=float(freebayes.AO_f_TOT)+float(freebayes.RO_f_TOT)
		freebayes.DP_r_TOT=float(freebayes.AO_r_TOT)+float(freebayes.RO_r_TOT)

		try:
			freebayes.A_QB=float(sample[format.index('QA'+suffisso)])/freebayes.AO
		except:
			freebayes.A_QB='.'
		try:	
			freebayes.R_QB=float(sample[format.index('QR'+suffisso)])/freebayes.RO
		except:
			freebayes.R_QB='.'
		try:	
			freebayes.QB=float(sample[format.index('QA'+suffisso)])/freebayes.AO
		except:
			freebayes.QB='.'


		try:
			gls = [float(x) for x in sample[format.index("GL"+suffisso)].split(",")]
			freebayes.lod = max(gls[i] - gls[0] for i in range(1, len(gls)))
		except:
			freebayes.lod = -1.0
		
		try:
			freebayes.AF=freebayes.AO/freebayes.DP
		except:
			freebayes.AF

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

		try:
			freebayes.FILTER=filter.split(',')[1]
		except:
			freebayes.FILTER=filter
		#print chrom,pos,ref,alt,freebayes.GT
	
	freebayes.Call=1

def get_info_GATK(chrom,pos,ref,alt,filter,info,format,sample,GATK,merged):
	'''estrae le informazioni dal vcf di GATK'''
	if merged==1:
		suffisso='_G'
	else:
		suffisso=''

	GATK.GT=sample[format.index('GT'+suffisso)]
	if GATK.GT=='./.':
		pass
	else:
		GATK.AO=float((sample[format.index('AD'+suffisso)]).split(',')[1])
		GATK.RO=float((sample[format.index('AD'+suffisso)]).split(',')[0])
		GATK.DP=float(sample[format.index('DP'+suffisso)])
		try:
			GATK.RF=(GATK.AO+GATK.RO)/GATK.DP
		except:
			GATK.RF='.'
		try:
			GATK.AF=GATK.AO/GATK.DP
		except:
			GATK.AF='.'
		GATK.STR='0'
		
		try:
			GATK.AO_r=float((sample[format.index('SB'+suffisso)]).split(',')[3])
			GATK.AO_f=float((sample[format.index('SB'+suffisso)]).split(',')[2])
			GATK.RO_r=float((sample[format.index('SB'+suffisso)]).split(',')[1])
			GATK.RO_f=float((sample[format.index('SB'+suffisso)]).split(',')[0])
			GATK.DP_r=GATK.AO_r+GATK.RO_r
			GATK.DP_f=GATK.AO_f+GATK.RO_f
		except:
			GATK.AO_r='.'
			GATK.AO_f='.'
			GATK.RO_r='.'
			GATK.RO_f='.'
			GATK.DP_r='.'
			GATK.DP_f='.'

		GATK.QB=sample[format.index('SQD'+suffisso)]


		for ind in info:
			if ind.startswith("AC"+suffisso+"="):
				GATK.AC=ind.split('=')[1]
			if ind.startswith("AF"+suffisso+"="):
				GATK.AF_TOT=ind.split('=')[1]
			if ind.startswith("AN"+suffisso+"="):
				GATK.AN=ind.split('=')[1]
			if ind.startswith("BaseQRankSum"+suffisso+"="):
				GATK.BaseQRankSum=ind.split('=')[1]
			if ind.startswith("ClippingRankSum"+suffisso+"="):
				GATK.ClippingRankSum=ind.split('=')[1]
			if ind.startswith("DP"+suffisso+"="):
				GATK.DP_TOT=ind.split('=')[1]
			if ind.startswith("DS"+suffisso+"="):
				GATK.DS=ind.split('=')[1]
			if ind.startswith("END"+suffisso+"="):
				GATK.END=ind.split('=')[1]
			if ind.startswith("ExcessHet"+suffisso+"="):
				GATK.ExcessHet=ind.split('=')[1]
			if ind.startswith("FS"+suffisso+"="):
				GATK.FS=ind.split('=')[1]
			if ind.startswith("HaplotypeScore"+suffisso+"="):
				GATK.HaplotypeScore=ind.split('=')[1]
			if ind.startswith("InbreedingCoeff"+suffisso+"="):
				GATK.InbreedingCoeff=ind.split('=')[1]
			if ind.startswith("MLEAC"+suffisso+"="):
				GATK.MLEAC=ind.split('=')[1]
			if ind.startswith("MLEAF"+suffisso+"="):
				GATK.MLEAF=ind.split('=')[1]
			if ind.startswith("MQ"+suffisso+"="):
				GATK.MQ=ind.split('=')[1]
			if ind.startswith("MQRankSum"+suffisso+"="):
				GATK.MQRankSum=ind.split('=')[1]
			if ind.startswith("QD"+suffisso+"="):
				GATK.QD=ind.split('=')[1]
			if ind.startswith("RAW_MQ"+suffisso+"="):
				GATK.RAW_MQ=ind.split('=')[1]
			if ind.startswith("ReadPosRankSum"+suffisso+"="):
				GATK.ReadPosRankSum=ind.split('=')[1]
			if ind.startswith("SOR"+suffisso+"="):
				GATK.SOR=ind.split('=')[1]

		try:
			if opts.amplicon:
				if min(GATK.DP_r,GATK.DP_f)/(GATK.DP_r+GATK.DP_f) >= 0.05:
					GATK.STRBIAS=1-stats.fisher_exact([[GATK.RO_f, GATK.RO_r], [GATK.AO_f, GATK.AO_r]])[1]
				else:
					GATK.STRBIAS='.'
			else:
				if min(GATK.DP_r,GATK.DP_f)/(GATK.DP_r+GATK.DP_f) > 0:
					GATK.STRBIAS=1-stats.fisher_exact([[GATK.RO_f, GATK.RO_r], [GATK.AO_f, GATK.AO_r]])[1]

				else:
					GATK.STRBIAS='.'
		except:
			GATK.STRBIAS='.'


		try:
			GATK.GQ=float(sample[format.index('GQ'+suffisso)])
		except:
			GATK.GQ='.'
		
		try:
			GATK.FILTER=filter.split(',')[0]
		except:
			GATK.FILTER=filter
		#print chrom,pos,ref,alt,GATK.GT
		GATK.Call=1
	
def get_info_varscan(chrom,pos,ref,alt,filter,info,format,sample,varscan,merged):
	'''estrae le informazioni dal vcf di varscan'''
	if merged==1:
		suffisso='_V'
	else:
		suffisso=''
	varscan.GT=sample[format.index('GT'+suffisso)]
	if varscan.GT== './.':
		#print chrom,pos,ref,alt,varscan.GT
		pass
	else:
		varscan.AO=float(sample[format.index('AD'+suffisso)])
		varscan.RO=float(sample[format.index('RD'+suffisso)])
		varscan.DP=float(sample[format.index('DP'+suffisso)])
		varscan.SDP=float(sample[format.index('SDP'+suffisso)])

		varscan.RF=(varscan.RO + varscan.AO)/varscan.SDP
	
	 	for ind in info:
	 		if ind.startswith("ADP"+suffisso+"="):
				varscan.ADP=ind.split('=')[1]
			if ind.startswith("WT"+suffisso+"="):
				varscan.WT=ind.split('=')[1]
			if ind.startswith("HET"+suffisso+"="):
				varscan.HET=ind.split('=')[1]
			if ind.startswith("HOM"+suffisso+"="):
				varscan.HOM=ind.split('=')[1]
			if ind.startswith("NC"+suffisso+"="):
				varscan.NC=ind.split('=')[1]
		
		varscan.AC=	(int(varscan.HET)+int(varscan.HOM)*2)
		varscan.AN=	2*(varscan.HET+varscan.HOM+varscan.WT)
					
		varscan.RO_f=float(sample[format.index('RDF'+suffisso)])
		varscan.RO_r=float(sample[format.index('RDR'+suffisso)])
		varscan.AO_f=float(sample[format.index('ADF'+suffisso)])
		varscan.AO_r=float(sample[format.index('ADR'+suffisso)])
		varscan.DP_f=varscan.RO_f + varscan.AO_f
		varscan.DP_r=varscan.RO_r + varscan.AO_r
		varscan.QB_R=float(sample[format.index('RBQ'+suffisso)])
		varscan.QB_A=float(sample[format.index('ABQ'+suffisso)])
		try:
			varscan.QB=varscan.QB_A
		except:
			varscan.QB='.'
		Varscan.GQ=float(sample[format.index('GQ'+suffisso)])
		
		try:
			varscan.AF=float(varscan.AO/(varscan.DP))
		except:
			varscan.AF='.'

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

		try:
			Varscan.FILTER=filter.split(',')[2]
		except:
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
		for varcall in varc_array[:-1]:
			if varcall is not "":
				vett_MBQ=vett_MBQ+[varcall.QB]
				vett_DP=vett_DP+[varcall.DP]
				vett_AO=vett_AO+[varcall.AO]
				vett_RO=vett_RO+[varcall.RO]
				vett_AC=vett_AC+[varcall.AC]
				
				if index == 0:

					features.RF_Freebayes=varc_array[0].RF
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

					features.RF_GATK=varc_array[2].RF
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
					
					features.RF_Varscan=varc_array[1].RF
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

		# features.SimpleRepeat_iEVA = varc_array[3].SimpleRepeat
		# features.SimpleRepeatLength_iEVA = varc_array[3].SimpleRepeatLength
		# features.SimpleRepeatUnit_iEVA= varc_array[3].SimpleRepeatUnit
		# features.PseudoNucleotidesComposition_iEVA = '\t'.join((varc_array[3].PseudoNucleotidesComposition).split(','))
		# features.RepeatMasker_iEVA = varc_array[3].RepeatMasker
		# features.gcContent_iEVA = varc_array[3].gcContent

		# features.VariantClass_iEVA = varc_array[3].VariantClass
		# features.StrandBiasReads_iEVA = varc_array[3].StrandBiasReads
		# features.UnMappedReads_iEVA = varc_array[3].UnMappedReads
		# features.MeanMappingQuality_iEVA = varc_array[3].MeanMappingQuality
		# features.MappingQualityZero_iEVA = varc_array[3].MappingQualityZero
		# features.NotPrimaryAlignment_iEVA = varc_array[3].NotPrimaryAlignment
		# features.SupplementaryAlignment_iEVA = varc_array[3].SupplementaryAlignment
		# features.NotPairedReads_iEVA = varc_array[3].NotPairedReads
		# features.NotProperPairedReads_iEVA = varc_array[3].NotProperPairedReads
		# features.AlignmentScore_iEVA = varc_array[3].AlignmentScore
		# features.TotalDupReads_iEVA = varc_array[3].TotalDupReads

		# features.NumberReadDupRef_iEVA = varc_array[3].NumberReadDupRef
		# features.NumberReadDupAlt_iEVA = varc_array[3].NumberReadDupAlt
		# features.DuplicateReference_iEVA = varc_array[3].DuplicateReference
		# features.DuplicateAlternate_iEVA = varc_array[3].DuplicateAlternate
		# features.DeltaDuplicate_iEVA = varc_array[3].DeltaDuplicate
		# features.iEvaDepth_iEVA = varc_array[3].iEvaDepth
		# features.iAlleleDepth_iEVA = varc_array[3].iAlleleDepth
		# features.ReadRef_iEVA = varc_array[3].ReadRef
		# features.ReadAlt_iEVA = varc_array[3].ReadAlt
		# features.MeanRefQscore_iEVA = varc_array[3].MeanRefQscore
		# features.MeanAltQscore_iEVA = varc_array[3].MeanAltQscore
		# features.RefMeanMappingQuality_iEVA = varc_array[3].RefMeanMappingQuality
		# features.AltMeanMappingQuality_iEVA = varc_array[3].AltMeanMappingQuality
		# features.TotalDPUnfilter_iEVA = varc_array[3].TotalDPUnfilter
		# features.NumberClippedReadsRef_iEVA = varc_array[3].NumberClippedReadsRef
		# features.NumberClippedReadsAlt_iEVA = varc_array[3].NumberClippedReadsAlt
		# features.ClippedReadsRef_iEVA = varc_array[3].ClippedReadsRef
		# features.ClippedReadsAlt_iEVA = varc_array[3].ClippedReadsAlt


		vett_RF_media=[features.RF_GATK,features.RF_Varscan,features.RF_Freebayes]
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
		RF_media='.'
		RF_mediana='.'

		v=[]
		for dp in vett_DP:
			if dp and dp is not '.':
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
			if ao and ao is not '.':
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
			if af and af is not '0'and af is not '.':
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
			features.STRBIAS_mediana = statistics.median(v)
		except:
			features.STRBIAS_mediana='.'

		v=[]
		for ac in vett_AC:
			if ac and ac is not '.':
				v=v+[int(ac)]

		try:
			features.AC_min=statistics.amin(v)
		except:
			features.AC_min='.'
		try:
			features.AC_max=statistics.amax(v)
		except:
			features.AC_max='.'				
		try:
			features.AC_media=statistics.mean(v)
		except:
			features.AC_media='.'
		try:
			features.AC_mediana= statistics.median(v)
		except:
			features.AC_mediana='.'

		v=[]
		for rf in vett_RF_media:
			if rf and rf is not '.':
				v=v+[float(rf)]
		try:
			features.RF_media=statistics.mean(v)
		except:
			features.RF_media='.'
		try:
			features.RF_mediana=statistics.median(v)
		except:
			features.RF_mediana='.'

		dictionary[variante]= varc_array + [features]


def switch(dictionary,ID,index,chrom,pos,ref,alt,filter,info,format,sample):
	'''tramite index richiama la funzione di estrazione delle informazioni del variant caller associato all'indice'''
	if index == None:
		vettore=['','','','']

		freebayes=Freebayes()
		get_info_Freebayes(chrom,pos,ref,alt,filter,info,format,sample,freebayes,1)
		if freebayes.GT != './.':
			vettore[0]=freebayes

		gatk=GATK()
		get_info_GATK(chrom,pos,ref,alt,filter,info,format,sample,gatk,1)
		if gatk.GT != './.':
			vettore[2]=gatk

		varscan=Varscan()
		get_info_varscan(chrom,pos,ref,alt,filter,info,format,sample,varscan,1)
		if varscan.GT != './.':
			vettore[1]=varscan

		ieva=iEVA()
		get_info_iEVA(chrom,pos,ref,alt,filter,info,format,sample,ieva)
		vettore[3]=ieva

		dictionary[ID]=vettore

	else:
		if dictionary.has_key(ID):
			vettore=dictionary[ID]
		else:
			vettore=['','','']

		if index==0:
			# print 'Freebayes'
			freebayes=Freebayes()
			get_info_Freebayes(chrom,pos,ref,alt,filter,info,format,sample,freebayes,0)
			if freebayes.GT != './.':
				vettore[0]=freebayes

		elif index==2:
			# print 'gatk'
			gatk=GATK()
			get_info_GATK(chrom,pos,ref,alt,filter,info,format,sample,gatk,0)
			if gatk.GT != './.':
				vettore[2]=gatk

		elif index==1 :
			# print 'varscan'
			varscan=Varscan()
			get_info_varscan(chrom,pos,ref,alt,filter,info,format,sample,varscan,0)
			if varscan.GT != './.':
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
			if line == 'PseudoNucleotidesComposition_iEVA':
				header=header+['\t'.join(['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])]
				features_variante=features_variante+['features.'+line]
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
			if opts.cf:
				if gt != './.' and gt != '.':
					dataset_varianti.write(var+ '\n')
					break
			else:
				if gt != './.' and gt != '.' and gt != '0/0':
					dataset_varianti.write(var+ '\n')
					break
	dataset_varianti.close()
	lista_features.close()
	
def print_vcf(varianti,out):
	dataset_varianti_vcf=open(out+ '/TOTAL.vcf','w')
	dataset_varianti_vcf.write('##fileformat=VCFv4.2\n'+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLES\n')
	for variante in varianti.keys():
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
	
	if 'FreeB' in vcf_name:
		variant_caller = 'FreeBayes'
	elif 'GATK' in vcf_name:
		variant_caller = 'GATK'
	elif 'VarScan' in vcf_name:
		variant_caller = 'VarScan'
	else:
		variant_caller = 'Merged'

	print '\n'+variant_caller +'\n'

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
		print sample

		try: 
			os.mkdir(opts.out_path +'/' + sample)
		except:
			pass
		
		sample_vcf = open(opts.out_path +'/' + sample +'/' + sample + '_'+variant_caller +'.vcf','w')
		
		sample_vcf.write('\n'.join(header) + '\n')
		sample_vcf.write('\t'.join(header_chrom[0:9] + [sample])  +'\n')
		gvcf = open(opts.gvcf_path +'/' + sample +'.g.vcf','r')
		all_gvcf=gvcf.readlines()
		gvcf.close()

		sample_gvcf = []

		for riga in all_gvcf:
			riga = riga.rstrip()
			if riga.startswith('#'):
				continue
			else:
				riga_split = riga.split('\t')
				if riga_split[4] != '<NON_REF>':
					sample_gvcf += [riga]

		for variante in varianti:
			
			variante_split = variante.split('\t')
			chrom = variante_split[0]
			pos = variante_split[1]
			id = variante_split[2]
			ref = variante_split[3]
			alt = variante_split[4]
			if variant_caller == 'GATK' or variant_caller == 'Merged':
				if variant_caller == 'Merged':
					suffisso='_G'
				else:
					suffisso=''
				sSB = '.'
				sQD = '.'
				for line in sample_gvcf:
					line = line.rstrip()
					if line.startswith(chrom+'\t'+pos) and line.split('\t')[4]:
						format = line.split('\t')[-2]
						sformat = line.split('\t')[-1]
						qual = line.split('\t')[5]
						info = line.split('\t')[7].split(';')
						try:
							ad = (sformat.split(':')[(format.split(':')).index('AD')]).split(',')
						except:
							ad = ['0','0','0']

						
						ad = map(float, ad)
						ad_sum = sum(ad[1:])

						for elem in info:
							if elem.startswith('DP='):
								DP = elem.split('=')[1]
								break
						try:
							sSB = sformat.split(':')[format.split(':').index('SB')]
							if ad_sum > 0.0:
								sQD = round(float(qual)/ad_sum ,2)
							else:
								sQD = 0.0
							break
						except:
							print 'ci sono problemi',sample,chrom,pos,format,sformat
		
				variante_common = variante_split[0:8] + [variante_split[8]+':SB'+suffisso+':SQD'+suffisso]
				format_sample = variante_split[header_chrom.index(sample)]  +':'+ sSB + ':' + str(sQD)

			else:
				variante_common = variante_split[0:9]
				format_sample = variante_split[header_chrom.index(sample)   ]
				
			sample_vcf.write('\t'.join(variante_common + [format_sample]) + '\n')
			
		sample_vcf.close()

def main():

	parser = argparse.ArgumentParser('Parse VCF output from Variant callers to output a variant_dataset.txt.  Output is to stdout.')
	parser.add_argument('-f', '--freebayes', help="Freebayes vcf output file name",default=None)
	parser.add_argument('-g', '--gatk', help="gatk vcf output file name",default=None)
	parser.add_argument('-v', '--varscan', help="Varscan vcf output file name",default=None)
	parser.add_argument('-m', '--merged', help="vcf merged from GATK,Freebayes and Varscan2",default=None)
	parser.add_argument('-l', '--listaFeatures', help="Lista di features da stampare",default=None)
	parser.add_argument('-s', '--split', help="Split vcf per samples", action='store_true')
	parser.add_argument('-F', '--feat_extraction', help="Enable features extraction", action='store_true')
	parser.add_argument('-a', '--amplicon',help="Amplicon design", action='store_true')
	parser.add_argument('-o', '--out_path',help="path di output")
	parser.add_argument('-G', '--gvcf_path',help="gvcf path")
	parser.add_argument('--cf',help="cell free", action='store_true')

	global opts 
	opts = parser.parse_args()
	try:
		os.mkdir(opts.out_path)
	except:
		pass

	if opts.merged != None:
		samples = samples_name_extract(open(opts.merged,'r'))
		if opts.split:
			print 'Splitto le varianti per campione.'
			split_vcf(opts.merged,samples)
			print'\nSplitto le varianti per campione:Done'

		if opts.feat_extraction:
			print '\nFEATURES EXTRACTION.'
			varianti_total = dict()
			for dir_sample in os.listdir(opts.out_path):
				varianti = dict()
				vcf_path = opts.out_path +'/' + dir_sample
				if os.path.isdir(vcf_path):
					for vcf_name in os.listdir(vcf_path):
						in_file = open(vcf_path + '/' + vcf_name)	
						vcfreader = read(in_file,None,varianti)

					set_features(varianti)
					print_var(varianti,opts.out_path,dir_sample)
				
				for var in varianti.keys():
					#varianti_total[var] = var.split('\t')[0]+'\t'+var.split('\t')[1]+'\t.\t'+var.split('\t')[2]+'\t'+var.split('\t')[3]+'\t.\t.\t.\t.\t.'
					varianti_total[var] = ''
			print_vcf(varianti_total,opts.out_path)

	else:
		callers = [opts.gatk,opts.varscan,opts.freebayes]
		samples = samples_name_extract(open(opts.freebayes,'r'))
		print callers
		
		
		if opts.split:
			print 'Splitto le varianti per campione.'
			for vcf_dir in callers:
				#print 'Splitto ' + vcf_dir
				split_vcf(vcf_dir,samples)
			print'\nSplitto le varianti per campione:Done'

		if opts.feat_extraction:
			print '\nFEATURES EXTRACTION.'
			varianti_total = dict()
			for dir_sample in os.listdir(opts.out_path):
				varianti = dict()
				vcf_path = opts.out_path +'/' + dir_sample
				if os.path.isdir(vcf_path):
					print "\nAnalizzo le varianti da: " + vcf_path
					for vcf_name in os.listdir(vcf_path) :
						print vcf_name
						if 'FreeB' in vcf_name:
							index = 0
						elif 'GATK' in vcf_name:
							index = 2
						elif 'VarScan' in vcf_name:
							index = 1
						
						in_file = open(vcf_path + '/' + vcf_name)
						vcfreader = read(in_file,index,varianti)
					
					set_features(varianti)
					print_var(varianti,opts.out_path,dir_sample)
				
				for var in varianti.keys():
					#varianti_total[var] = var.split('\t')[0]+'\t'+var.split('\t')[1]+'\t.\t'+var.split('\t')[2]+'\t'+var.split('\t')[3]+'\t.\t.\t.\t.\t.'
					varianti_total[var] = ''
			print_vcf(varianti_total,opts.out_path)
main()
