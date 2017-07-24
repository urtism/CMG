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
	QUAL=''
	Call='0'

class Freebayes(Caller):
	pass
	
class Samtools(Caller):
	pass

class Varscan(Caller):
	pass

class Features():

	GT_Varscan='.'
	GT_Freebayes='.'
	GT_Samtools='.'

	QB_Samtools='.'
	QB_Varscan='.'
	QB_Freebayes='.'
	MBQ_media='.'
	MBQ_mediana='.'
	
	CallSamtools='0'
	CallVarscan='0'
	CallFreebayes='0'

	FILTER_Samtools='.'
	FILTER_Freebayes='.'
	FILTER_Varscan='.'


''' //////////// FUNZIONI ////////////'''

def get_info_Freebayes(chrom,pos,ref,alt,qual,filter,info,format,sample,freebayes):
	'''estrae le informazioni dal vcf di freebayes'''
	
	freebayes.GT=sample[format.index('GT')]
	if freebayes.GT=='.' :
		freebayes.GT='./.'
	else:
		freebayes.QB=qual
		freebayes.FILTER=filter
	freebayes.Call=1

def get_info_Samtools(chrom,pos,ref,alt,qual,filter,info,format,sample,samtools):
	'''estrae le informazioni dal vcf di Samtools'''
	
	samtools.GT=sample[format.index('GT')]
	if samtools.GT=='./.':
		pass
	else:
		samtools.QB=qual		
		samtools.FILTER=filter
	samtools.Call=1
	
def get_info_varscan(chrom,pos,ref,alt,qual,filter,info,format,sample,varscan):
	'''estrae le informazioni dal vcf di varscan'''
	
	varscan.GT=sample[format.index('GT')]
	if varscan.GT== './.':
		pass
	else:
		varscan.Call=1
		varscan.QB_R=float(sample[format.index('RBQ')])
		varscan.QB_A=float(sample[format.index('ABQ')])
		try:
			varscan.QB=varscan.QB_A
		except:
			varscan.QB='.'
		Varscan.FILTER=filter
	Varscan.Call=1

def set_features(dictionary):
	'''setta i valori delle features in base alle info estratte dai vcf'''
	for variante in dictionary.keys():
		features=Features()
		varc_array=dictionary.get(variante)
		
		vett_MBQ=[]
		index=0
		#print variante
		#print varc_array
		for varcall in varc_array:
			if varcall is not "":
				vett_MBQ=vett_MBQ+[varcall.QB]
							
				if index == 0:
			
					features.GT_Freebayes=varc_array[0].GT
					features.QB_Freebayes=varc_array[0].QB
					features.CallFreebayes=varc_array[0].Call
					features.FILTER_Freebayes=varc_array[0].FILTER

				if index == 2:

					features.GT_Samtools=varc_array[2].GT
					features.QB_Samtools=varc_array[2].QB
					features.CallSamtools=varc_array[2].Call
					features.FILTER_Samtools=varc_array[2].FILTER

				elif index == 1:
					
					features.GT_Varscan=varc_array[1].GT
					features.QB_Varscan=varc_array[1].QB
					features.CallVarscan=varc_array[1].Call
					features.FILTER_Varscan=varc_array[1].FILTER

			index = index + 1	
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
		dictionary[variante]=varc_array + [features]

def switch(dictionary,ID,index,chrom,pos,ref,alt,qual,filter,info,format,sample):
	'''tramite index richiama la funzione di estrazione delle informazioni del variant caller associato all'indice'''
	if dictionary.has_key(ID):
		vettore=dictionary[ID]
	else:
		vettore=['','','']

	if index==0:
		# print 'Freebayes'
		freebayes=Freebayes()
		get_info_Freebayes(chrom,pos,ref,alt,qual,filter,info,format,sample,freebayes)
		if freebayes.GT != './.':
			vettore[0]=freebayes

	elif index==2:
		# print 'Samtools'
		samtools=Samtools()
		get_info_Samtools(chrom,pos,ref,alt,qual,filter,info,format,sample,samtools)
		if samtools.GT != './.':
			vettore[2]=samtools
	elif index==1 :
		# print 'varscan'
		varscan=Varscan()
		get_info_varscan(chrom,pos,ref,alt,qual,filter,info,format,sample,varscan)
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
			qual=var[5]
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
				switch(dictionary,ID,index,chrom,pos,ref,alt,qual,filter,info,format,sample)
				
def control(dictionary):
	''' esegue un controllo sulle varianti, se non hanno variant caller che le chiama vengono eliminate'''
	for variante in dictionary.keys():
		if dictionary[variante][:3] == ['','','']:
			#print "sto cancellando:",variante
			del dictionary[variante]
				
def print_var(dictionary,out,sample_name):

	dataset_varianti=open(out + '/' + sample_name + '.tsv','w')
	print out + '/' + sample_name + '.tsv'
	header=[]
	features_variante=[]

	header=["CHROM","POS","ID","REF","ALT","CallSamtools","CallFreebayes","CallVarscan","GT_Samtools","GT_Freebayes","GT_Varscan","QUAL_media"]
	#print features_variante
	dataset_varianti.write('\t'.join(header)+ '\n')
	for variante in dictionary.keys():
		features = dictionary.get(variante)[-1]
		features_variante=[str(features.CallSamtools),str(features.CallFreebayes),str(features.CallVarscan),str(features.GT_Samtools),str(features.GT_Freebayes),str(features.GT_Varscan),str(features.MBQ_media)]
		var=variante.split('\t')[0]+'\t'+variante.split('\t')[1]+'\t' +sample_name +'\t'+variante.split('\t')[2]+'\t'+variante.split('\t')[3]+ '\t' + '\t'.join(features_variante)
		
		for gt in [features.GT_Samtools,features.GT_Varscan,features.GT_Freebayes]:
			if gt != './.' and gt != '.' and gt != '0/0':
				dataset_varianti.write(var+ '\n')
				break
	dataset_varianti.close()

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
	elif 'Samtools' in vcf_name:
			variant_caller = 'Samtools'
	elif 'VarScan' in vcf_name:
			variant_caller = 'VarScan'
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
		
		for variante in varianti:
			
			variante_split = variante.split('\t')
			chrom = variante_split[0]
			pos = variante_split[1]
			id = variante_split[2]
			ref = variante_split[3]
			alt = variante_split[4]
			variante_common = variante_split[0:9]
			format_sample = variante_split[header_chrom.index(sample)   ]
			sample_vcf.write('\t'.join(variante_common + [format_sample]) + '\n')
			
		sample_vcf.close()

def main():

	parser = argparse.ArgumentParser('Parse Sanger VCF output from Variant callers to output a variant_dataset.txt.  Output is to stdout.')
	parser.add_argument('-f', '--freebayes', help="Freebayes vcf output file name",default=None)
	parser.add_argument('-g', '--samtools', help="Samtools vcf output file name",default=None)
	parser.add_argument('-v', '--varscan', help="Varscan vcf output file name",default=None)
	parser.add_argument('-s', '--split', help="Split vcf per samples", action='store_true')
	parser.add_argument('-F', '--feat_extraction', help="Enable features extraction", action='store_true')
	parser.add_argument('-o', '--out_path',help="path di output")

	global opts 
	opts = parser.parse_args()
	callers = [opts.samtools,opts.varscan,opts.freebayes]
	samples = samples_name_extract(open(opts.freebayes,'r'))
	#print callers
	try:
		os.mkdir(opts.out_path)
	except:
		pass
	
	if opts.split:
		print 'Splitto le varianti per campione.'
		for vcf_dir in callers:
			#print 'Splitto ' + vcf_dir
			split_vcf(vcf_dir,samples)
		print'\nSplitto le varianti per campione:Done'

	if opts.feat_extraction:
		print '\nFEATURES EXTRACTION.'
		for sample in samples:
			varianti_total = dict()
			varianti = dict()
			for vcf_name in callers :
				print vcf_name
				if 'FreeB' in vcf_name:
					index = 0
				elif 'Samtools' in vcf_name:
					index = 2
				elif 'VarScan' in vcf_name:
					index = 1
				
				in_file = open(vcf_name,'r')
				vcfreader = read(in_file,index,varianti)
			
			set_features(varianti)
			print_var(varianti,opts.out_path,sample)
			
			#for var in varianti.keys():
				#varianti_total[var] = var.split('\t')[0]+'\t'+var.split('\t')[1]+'\t.\t'+var.split('\t')[2]+'\t'+var.split('\t')[3]+'\t.\t.\t.\t.\t.'
			#	varianti_total[var] = ''
		#print_vcf(varianti_total,opts.out_path)
main()
