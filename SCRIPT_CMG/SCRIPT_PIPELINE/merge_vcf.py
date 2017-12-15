import re
import string
import argparse
import sys
import statistics
import scipy.stats as stats

def parse_header(vcf,i):
	h_filter=''
	h_info=''
	h_format=''
	var_format='.'

	if i==1:
		add='G'
	elif i==2:
		add='F'
	elif i==3:
		add='V'

	#print i,add
	for line in vcf:
		line=line.rstrip()
		line_split=line.split(',')

		#print line

		if line.startswith('##FILTER'):
			line_split[0]=line_split[0]+'_'+ add
			h_filter=h_filter+','.join(line_split) + '\n'
		
		elif line.startswith('##INFO'):
			line_split[0]=line_split[0]+'_'+ add
			h_info=h_info+','.join(line_split) + '\n'

		elif line.startswith('##FORMAT'):
			line_split[0]=line_split[0]+'_'+ add
			h_format=h_format+','.join(line_split) + '\n'

		elif line.startswith('chr'):
			line_split=line.split('\t')
			#print line_split
			#print line_split[8]
			var_format=modifica_format(line_split[8],i)
		#elif line.startswith('#CHROM') and i==1:
		

	return h_filter,h_info,h_format,var_format

def samples(vcf):
	sample_list=''
	for line in vcf:
		line=line.rstrip()
		line_split=line.split('\t')
		
		if line.startswith('#CHROM'):
			sample_list=line_split[9:]
	return sample_list

def estrai_varianti(varianti,vcf,i,sample_list):
	for line in vcf:
		line=line.rstrip()
		line_split=line.split('\t')
		if line.startswith('##'):
			continue
		elif line.startswith('#CHROM'):
			h_sample_list=line_split[9:]
		else:
			var='\t'.join(line_split[0:5])
			qual=line_split[5]
			filter=modifica_filter(line_split[6],i)
			info=modifica_info(line_split[7],i)
			format_samples=riordina_samples(line_split[9:],sample_list,h_sample_list)
			if var in varianti.keys():
				line_var=varianti.get(var)
				#print line_var[2]
			else:
				line_var=[['.','.','.','.'],['.','.','.','.'],['.','.','.','.']]
				
			(line_var[i-1])[0]=qual
			(line_var[i-1])[1]=filter
			(line_var[i-1])[2]=info
			(line_var[i-1])[3]=format_samples
			varianti[var]=line_var


def riordina_samples(format_list,sample_list,h_sample_list):
	new_format_list=[]
	for sample in sample_list:
		i=h_sample_list.index(sample)
		new_format_list+=[format_list[i]]
	return new_format_list

def modifica_filter(filter,i):
	if i==1:
		add='G'
	elif i==2:
		add='F'
	elif i==3:
		add='V'
	if filter != '.':
		filter_split=filter.split(';')
		filter=';'.join(a+'_'+add for a in filter_split)
	return filter

def modifica_info(info,i):
	if i==1:
		add='G'
	elif i==2:
		add='F'
	elif i==3:
		add='V'
	info_split=info.split(';')
	for el in info_split:
		el_split=el.split('=')
		info_split[info_split.index(el)]=el_split[0]+'_'+add+'='+el_split[1]

	info=';'.join(info_split)
	return info

def modifica_format(format,i):
	if i==1:
		add='G'
	elif i==2:
		add='F'
	elif i==3:
		add='V'
	format_split=format.split(':')
	format=':'.join(a+'_'+add for a in format_split)
	return format

def fix_format_values(line_var,s_list,i,format):
	format_sample=[]
	null_format='./.'
	if (line_var[i])[-1]=='.':
		for a in format.split(':')[1:]:
			#print a
			null_format=null_format + ':.'
		for s in s_list:
			format_sample=format_sample + [null_format]
	else:
		for form in (line_var[i])[-1]:
			null_format='./.'
			if form.startswith('.:'):
				for a in format.split(':')[1:]:
					null_format=null_format + ':.'
				format_sample = format_sample + [null_format]
			else:
				format_sample =	format_sample + [form]

	return format_sample


def main():
	
	parser = argparse.ArgumentParser('Parse 3 VCF output from Variant callers to output a VCF whith merged h_info. Output is to --out.')
	parser.add_argument('-g', '--gatk', help="GATK vcf output file name")
	parser.add_argument('-f', '--freebayes', help="Freebayes vcf output file name")
	parser.add_argument('-v', '--varscan', help="Varscan vcf output file name")
	parser.add_argument('-o', '--out', help="file name in output. It returns file_name.vcf ")
	
	global opts 
	opts = parser.parse_args()

	gatk=open(opts.gatk,'r').readlines()
	freebayes=open(opts.freebayes,'r').readlines()
	varscan=open(opts.varscan,'r').readlines()
	out_vcf=open(opts.out,'w')
	varianti = dict()
	HEADER='##file_format=VCFv4.2\n'

	h_filter_gatk,h_info_gatk,h_format_gatk,format_gatk=parse_header(gatk,1)
	h_filter_freebayes,h_info_freebayes,h_format_freebayes,format_freebayes=parse_header(freebayes,2)
	h_filter_varscan,h_info_varscan,h_format_varscan,format_varscan=parse_header(varscan,3)

	sample_list_gatk=samples(gatk)
	sample_list_freebayes=samples(freebayes)
	sample_list_varscan=samples(varscan)

	estrai_varianti(varianti,gatk,1,sample_list_gatk)
	estrai_varianti(varianti,freebayes,2,sample_list_gatk)
	estrai_varianti(varianti,varscan,3,sample_list_gatk)

	if format_varscan == '':
		format_varscan='GT_V:GQ_V:SDP_V:DP_V:RD_V:AD_V:FREQ_V:PVAL_V:RBQ_V:ABQ_V:RDF_V:RDR_V:ADF_V:ADR_V'



	HEADER=HEADER+h_filter_gatk+h_filter_freebayes+h_filter_varscan+h_info_gatk+h_info_freebayes+h_info_varscan+h_format_gatk+h_format_freebayes+h_format_varscan+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(sample_list_gatk)
	out_vcf.write(HEADER+'\n')
	for var in varianti.keys():
		line_var=varianti.get(var)
		qual=','.join([(line_var[0])[0],(line_var[1])[0],(line_var[2])[0]])
		filter=','.join([(line_var[0])[1],(line_var[1])[1],(line_var[2])[1]])
		info=','.join([(line_var[0])[2],(line_var[1])[2],(line_var[2])[2]])
		format=':'.join([format_gatk,format_freebayes,format_varscan])


		#print var

		(line_var[0])[-1]=fix_format_values(line_var,sample_list_gatk,0,format_gatk)
		(line_var[1])[-1]=fix_format_values(line_var,sample_list_gatk,1,format_freebayes)
		(line_var[2])[-1]=fix_format_values(line_var,sample_list_gatk,2,format_varscan)


		#print (line_var[0])[-1],(line_var[1])[-1],(line_var[2])[-1]

		format_sample=[]
		for sample in sample_list_gatk:
			i=sample_list_gatk.index(sample)
			format_sample+=[':'.join([((line_var[0])[-1])[i],((line_var[1])[-1])[i],((line_var[2])[-1])[i]])]
		
		format_samples='\t'.join(format_sample)

		out_vcf.write('\t'.join([var,qual,filter,info,format,format_samples])+'\n')
	
main()
