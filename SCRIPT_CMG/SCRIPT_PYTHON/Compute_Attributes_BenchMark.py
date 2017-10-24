import argparse
import re
import string
import sys

#------------------------------------------FUNCTION---------------------------------

def extract_FREE(New_Head,tag):

	# Tolgo le tag che leggo nelle if dalla nuova header
	if 'FREE-FORMAT-DPR-TOT' == tag or 'FREE-FORMAT-DPR-ALT' == tag or 'FREE-FORMAT-AD-ALT' == tag \
	or 'FREE-FORMAT-AD-REF' == tag or 'FREE-INFO-DPB' == tag or 'FREE-INFO-NUMALT' == tag or 'FREE-INFO-RPR' == tag \
	or 'FREE-INFO-RPL' == tag or 'FREE-INFO-RUN' == tag or 'FREE-INFO-technology.ILLUMINA' == tag \
	or 'FREE-INFO-SAF' == tag or 'FREE-INFO-SAR' == tag or 'FREE-INFO-SRF' == tag or 'FREE-INFO-SRR' == tag:
		pass

	#Aggioungo all'header le tag che devo calcolare
	elif 'FREE-FORMAT-RO' == tag:
		New_Head += ['FREE-FORMAT-RO']
		New_Head += ['FREE-FORMAT-AD-FREQ']

	elif 'FREE-INFO-AO' == tag:
		New_Head += ['FREE-INFO-AO']
		New_Head += ['FREE-INFO-AD-FREQ']

	else:
		New_Head += [tag]

	return New_Head

def extract_GATK(New_Head,tag):

 	if 'GATK-FORMAT-AD-ALT' == tag:
		New_Head += ['GATK-FORMAT-AD-ALT']
		New_Head += ['GATK-FORMAT-AD-FREQ']

	else:
		New_Head += [tag]

 	return New_Head

def extract_iEVA(New_Head,tag):

	if 'iEVA-FORMAT-iAD-REF' == tag:
		New_Head += ['iEVA-FORMAT-iAD-REF']

	elif 'iEVA-FORMAT-iAD-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iAD-ALT']
		New_Head += ['iEVA-FORMAT-iAD-FREQ']

	elif 'iEVA-FORMAT-iADup-REF' == tag:
		New_Head += ['iEVA-FORMAT-iADup-REF']
		New_Head += ['iEVA-FORMAT-iADup-REF-FREQ']

	elif 'iEVA-FORMAT-iADup-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iADup-ALT']
		New_Head += ['iEVA-FORMAT-iADup-ALT-FREQ']
		New_Head += ['iEVA-FORMAT-iADup-DELTA']

	elif 'iEVA-FORMAT-iAS-REF' == tag:
		New_Head += ['iEVA-FORMAT-iAS-REF']

	elif 'iEVA-FORMAT-iAS-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iAS-ALT']
		New_Head += ['iEVA-FORMAT-iAS-DELTA']

	elif 'iEVA-FORMAT-iXS-REF' == tag:
		New_Head += ['iEVA-FORMAT-iXS-REF']

	elif 'iEVA-FORMAT-iXS-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iXS-ALT']
		New_Head += ['iEVA-FORMAT-iXS-DELTA']

	elif 'iEVA-FORMAT-iXS0-REF' == tag:
		New_Head += ['iEVA-FORMAT-iXS0-REF']
		New_Head += ['iEVA-FORMAT-iXS0-REF-FREQ']

	elif 'iEVA-FORMAT-iXS0-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iXS0-ALT']
		New_Head += ['iEVA-FORMAT-iXS0-ALT-FREQ']
		New_Head += ['iEVA-FORMAT-iXS0-DELTA']

	elif 'iEVA-FORMAT-iUnMap-REF' == tag:
		New_Head += ['iEVA-FORMAT-iUnMap-REF']
		New_Head += ['iEVA-FORMAT-iUnMap-REF-FREQ']

	elif 'iEVA-FORMAT-iUnMap-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iUnMap-ALT']
		New_Head += ['iEVA-FORMAT-iUnMap-ALT-FREQ']
		New_Head += ['iEVA-FORMAT-iUnMap-DELTA']

	elif 'iEVA-FORMAT-iMQ-REF' == tag:
		New_Head += ['iEVA-FORMAT-iMQ-REF']

	elif 'iEVA-FORMAT-iMQ-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iMQ-ALT']
		New_Head += ['iEVA-FORMAT-iMQ-DELTA']

	elif 'iEVA-FORMAT-iMQ0-REF' == tag:
		New_Head += ['iEVA-FORMAT-iMQ0-REF']
		New_Head += ['iEVA-FORMAT-iMQ0-REF-FREQ']

	elif 'iEVA-FORMAT-iMQ0-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iMQ0-ALT']
		New_Head += ['iEVA-FORMAT-iMQ0-ALT-FREQ']
		New_Head += ['iEVA-FORMAT-iMQ0-DELTA']

	elif 'iEVA-FORMAT-iNPA-REF' == tag:
		New_Head += ['iEVA-FORMAT-iNPA-REF']
		New_Head += ['iEVA-FORMAT-iNPA-REF-FREQ']

	elif 'iEVA-FORMAT-iNPA-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iNPA-ALT']
		New_Head += ['iEVA-FORMAT-iNPA-ALT-FREQ']
		New_Head += ['iEVA-FORMAT-iNPA-DELTA']

	elif 'iEVA-FORMAT-iQual-REF' == tag:
		New_Head += ['iEVA-FORMAT-iQual-REF']

	elif 'iEVA-FORMAT-iQual-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iQual-ALT']
		New_Head += ['iEVA-FORMAT-iQual-DELTA']

	elif 'iEVA-FORMAT-iSA-REF' == tag:
		New_Head += ['iEVA-FORMAT-iSA-REF']
		New_Head += ['iEVA-FORMAT-iSA-REF-FREQ']

	elif 'iEVA-FORMAT-iSA-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iSA-ALT']
		New_Head += ['iEVA-FORMAT-iSA-ALT-FREQ']
		New_Head += ['iEVA-FORMAT-iSA-DELTA']

	elif 'iEVA-FORMAT-iNP-REF' == tag:
		New_Head += ['iEVA-FORMAT-iNP-REF']
		New_Head += ['iEVA-FORMAT-iNP-REF-FREQ']

	elif 'iEVA-FORMAT-iNP-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iNP-ALT']
		New_Head += ['iEVA-FORMAT-iNP-ALT-FREQ']
		New_Head += ['iEVA-FORMAT-iNP-DELTA']

	elif 'iEVA-FORMAT-iNPP-REF' == tag:
		New_Head += ['iEVA-FORMAT-iNPP-REF']
		New_Head += ['iEVA-FORMAT-iNPP-REF-FREQ']

	elif 'iEVA-FORMAT-iNPP-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iNPP-ALT']
		New_Head += ['iEVA-FORMAT-iNPP-ALT-FREQ']
		New_Head += ['iEVA-FORMAT-iNPP-DELTA']

	elif 'iEVA-FORMAT-iCR-REF' == tag:
		New_Head += ['iEVA-FORMAT-iCR-REF']
		New_Head += ['iEVA-FORMAT-iCR-REF-FREQ']

	elif 'iEVA-FORMAT-iCR-ALT' == tag:
		New_Head += ['iEVA-FORMAT-iCR-ALT']
		New_Head += ['iEVA-FORMAT-iCR-ALT-FREQ']
		New_Head += ['iEVA-FORMAT-iCR-DELTA']

	else:
		New_Head += [tag]

	return New_Head


def Extract_Feature_FreeBayes(line,Header_input,New_Head):

	Variant_Info = {}

	for elem in New_Head:

		if 'FREE-FORMAT-AD-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('FREE-FORMAT-AO')])/(float(line[Header_input.index('FREE-FORMAT-AO')])+float(line[Header_input.index('FREE-FORMAT-RO')])),4))
			except:
				Variant_Info[elem] = '?'

		elif 'FREE-INFO-AD-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('FREE-INFO-AO')])/(float(line[Header_input.index('FREE-INFO-DP')])),4))
			except:
				Variant_Info[elem] = '?'

		elif 'iEVA' in elem:
			continue

		else:
			try:
				Variant_Info[elem] = line[Header_input.index(elem)]
			except:
				continue

	return Variant_Info


def Extract_Feature_GATK(line,Header_input,New_Head):

	Variant_Info = {}

	for elem in New_Head:

		if 'GATK-FORMAT-AD-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('GATK-FORMAT-AD-ALT')])/(float(line[Header_input.index('GATK-FORMAT-AD-REF')])+float(line[Header_input.index('GATK-FORMAT-AD-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

		elif 'iEVA' in elem:
			continue

		else:
			try:
				Variant_Info[elem] = line[Header_input.index(elem)]
			except:
				continue

	return Variant_Info


def Extract_Feature_iEVA(Variant_Info,line,Header_input,New_Head):

	for elem in New_Head:

#------FORMAT-iAD-FREQ

		if 'iEVA-FORMAT-iAD-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iAD-ALT')])/(float(line[Header_input.index('iEVA-FORMAT-iDP')])),4))
			except:
				Variant_Info[elem] = '?'

#------FORMAT-iADup

		if 'iEVA-FORMAT-iADup-REF-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iADup-REF')])/(float(line[Header_input.index('iEVA-FORMAT-iADup-REF')])+float(line[Header_input.index('iEVA-FORMAT-iAD-REF')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iADup-ALT-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iADup-ALT')])/(float(line[Header_input.index('iEVA-FORMAT-iADup-ALT')])+float(line[Header_input.index('iEVA-FORMAT-iAD-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

#------FORMAT-iXS0

		if 'iEVA-FORMAT-iXS0-REF-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iXS0-REF')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-REF')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iXS0-ALT-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iXS0-ALT')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

#------FORMAT-iUnMap

		if 'iEVA-FORMAT-iUnMap-REF-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iUnMap-REF')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-REF')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iUnMap-ALT-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iUnMap-ALT')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

#-----FORMAT-iMQ0

		if 'iEVA-FORMAT-iMQ0-REF-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iMQ0-REF')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-REF')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iMQ0-ALT-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iMQ0-ALT')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

#-----FORMAT-iNPA

		if 'iEVA-FORMAT-iNPA-REF-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iNPA-REF')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-REF')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iNPA-ALT-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iNPA-ALT')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

#-----FORMAT-iSA

		if 'iEVA-FORMAT-iSA-REF-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iSA-REF')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-REF')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iSA-ALT-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iSA-ALT')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

#-----FORMAT-iNP

		if 'iEVA-FORMAT-iNP-REF-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iNP-REF')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-REF')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iNP-ALT-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iNP-ALT')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

#-----FORMAT-iNPP

		if 'iEVA-FORMAT-iNPP-REF-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iNPP-REF')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-REF')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iNPP-ALT-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iNPP-ALT')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

#-----FORMAT-iCR

		if 'iEVA-FORMAT-iCR-REF-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iCR-REF')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-REF')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iCR-ALT-FREQ' == elem:

			try:
				Variant_Info[elem] = str(round(float(line[Header_input.index('iEVA-FORMAT-iCR-ALT')])/(float(line[Header_input.index('iEVA-FORMAT-iAD-ALT')])),4))
			except:
				Variant_Info[elem] = '?'


		else:
			if 'DELTA' in elem:
				continue
			else:
				try:
					Variant_Info[elem] = line[Header_input.index(elem)]
				except:
					continue



	# Qui calcolo invece le features dei DELTA. le metto fuori perche mi servono prima i FREQ calcolati
	for elem in New_Head:

		if 'iEVA-FORMAT-iADup-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iADup-REF-FREQ')])-float(line[Header_input.index('iEVA-FORMAT-iADup-ALT-FREQ')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iADup-REF-FREQ')])+float(line[Header_input.index('iEVA-FORMAT-iADup-ALT-FREQ')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iAS-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iAS-REF')])-float(line[Header_input.index('iEVA-FORMAT-iAS-ALT')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iAS-REF')])+float(line[Header_input.index('iEVA-FORMAT-iAS-ALT')])),4))
			except:
				Variant_Info[elem] = '?'
		
		if 'iEVA-FORMAT-iXS-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iXS-REF')])-float(line[Header_input.index('iEVA-FORMAT-iXS-ALT')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iXS-REF')])+float(line[Header_input.index('iEVA-FORMAT-iXS-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iXS0-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iXS0-REF-FREQ')])-float(line[Header_input.index('iEVA-FORMAT-iXS0-ALT-FREQ')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iXS0-REF-FREQ')])+float(line[Header_input.index('iEVA-FORMAT-iXS0-ALT-FREQ')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iUnMap-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iUnMap-REF-FREQ')])-float(line[Header_input.index('iEVA-FORMAT-iUnMap-ALT-FREQ')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iUnMap-REF-FREQ')])+float(line[Header_input.index('iEVA-FORMAT-iUnMap-ALT-FREQ')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iMQ-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iMQ0-REF-FREQ')])-float(line[Header_input.index('iEVA-FORMAT-iMQ0-ALT-FREQ')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iMQ0-REF-FREQ')])+float(line[Header_input.index('iEVA-FORMAT-iMQ0-ALT-FREQ')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iMQ0-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iMQ0-REF-FREQ')])-float(line[Header_input.index('iEVA-FORMAT-iMQ0-ALT-FREQ')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iMQ0-REF-FREQ')])+float(line[Header_input.index('iEVA-FORMAT-iMQ0-ALT-FREQ')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iNPA-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iNPA-REF-FREQ')])-float(line[Header_input.index('iEVA-FORMAT-iNPA-ALT-FREQ')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iNPA-REF-FREQ')])+float(line[Header_input.index('iEVA-FORMAT-iNPA-ALT-FREQ')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iQual-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iQual-REF')])-float(line[Header_input.index('iEVA-FORMAT-iQual-ALT')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iQual-REF')])+float(line[Header_input.index('iEVA-FORMAT-iQual-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iSA-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iSA-REF')])-float(line[Header_input.index('iEVA-FORMAT-iSA-ALT')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iSA-REF')])+float(line[Header_input.index('iEVA-FORMAT-iSA-ALT')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iNP-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iNP-REF-FREQ')])-float(line[Header_input.index('iEVA-FORMAT-iNP-ALT-FREQ')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iNP-REF-FREQ')])+float(line[Header_input.index('iEVA-FORMAT-iNP-ALT-FREQ')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iNPP-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iNPP-REF-FREQ')])-float(line[Header_input.index('iEVA-FORMAT-iNPP-ALT-FREQ')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iNPP-REF-FREQ')])+float(line[Header_input.index('iEVA-FORMAT-iNPP-ALT-FREQ')])),4))
			except:
				Variant_Info[elem] = '?'

		if 'iEVA-FORMAT-iCR-DELTA' == elem:

			try:
				Variant_Info[elem] = str(round((float(line[Header_input.index('iEVA-FORMAT-iCR-REF-FREQ')])-float(line[Header_input.index('iEVA-FORMAT-iCR-ALT-FREQ')])) \
					/(float(line[Header_input.index('iEVA-FORMAT-iCR-REF-FREQ')])+float(line[Header_input.index('iEVA-FORMAT-iCR-ALT-FREQ')])),4))
			except:
				Variant_Info[elem] = '?'


	return Variant_Info




#-----------------------------------------MAIN----------------------------------


def main():
	
	parser = argparse.ArgumentParser('\n\nQuesto tool prende il formato csv in uscita dal tool UNIFORM-BenchMark e BenchMark-Extractor per calcolare le nuove feature di iEVA e dei variant caller')
	parser.add_argument('-I','--input',help="File csv contenente le varianti estratte da Benchmark-Extractor. Per funzionare, il tool DEVE AVERE IL NOME DEL VARIANT CALLER (VarScan, FreeBayes, GATK, indipendentemente da lettere maiuscole o minuscole) nel nome del file. Indicare come PATH/TO/NOME_FILE")
	parser.add_argument('-O','--outfile',help="File di output in csv format. Indicare come path/to/nome_file.csv")

	global opts

	opts = parser.parse_args()
	out = open(opts.outfile,'w')
	New_Head = []
	stamp = []
	new_line = {}

	with open(opts.input,'r') as csv:
		for line in csv:

			#Verifico che header contenga chr pos ref alt e id:
			if 'CHROM' and 'POS' and 'REF' and 'ALT' and 'ID' in line:

				#Salvo l'header
				line = line.rstrip()
				Header_input = line.split(',')

				for tag in Header_input:

					#Aggiungo le tag che splitto in FREE nella head
					if 'FreeBayes'.lower() in opts.input.lower():

						if 'FREE' in tag:
							New_Head = extract_FREE(New_Head,tag)

						elif 'iEVA' in tag:
							New_Head = extract_iEVA(New_Head,tag)

						else:
							New_Head += [tag]


					#Aggiungo le tag che splitto in VARSCAN nella head
					if 'VarScan'.lower() in opts.input.lower():

						if 'iEVA' in tag:
							New_Head = extract_iEVA(New_Head,tag)

						else:
							New_Head += [tag]


					#Aggiungo le tag che splitto in GATK nella head
					if 'GATK'.lower() in opts.input.lower():

						if 'GATK' in tag:
							New_Head = extract_GATK(New_Head,tag)

						elif 'iEVA' in tag:
							New_Head = extract_iEVA(New_Head,tag)

						else:
							New_Head += [tag]
	
				out.write(','.join(New_Head))


			#Altrimenti vado a modificare i dati
			else:

				def_line = []
				line = line.rstrip()
				line = line.split(',')


				if 'FreeBayes'.lower() in opts.input.lower():

					new_line = Extract_Feature_FreeBayes(line,Header_input,New_Head)

					new_line = Extract_Feature_iEVA(new_line,line,Header_input,New_Head)


				if 'VarScan'.lower() in opts.input.lower():

					new_line = Extract_Feature_iEVA(new_line,line,Header_input,New_Head)


				if 'GATK'.lower() in opts.input.lower():

					new_line = Extract_Feature_GATK(line,Header_input,New_Head)

					new_line = Extract_Feature_iEVA(new_line,line,Header_input,New_Head)

					#new_line = [item for val in line for item in val.split(',')]

				for attr in New_Head:

					stamp += [new_line[attr]]

			out.write(','.join(stamp) + '\n')
			stamp = []



	out.close()


main()