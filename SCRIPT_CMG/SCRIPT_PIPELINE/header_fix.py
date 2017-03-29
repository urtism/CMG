import argparse
import sys
from tr import tr
import re


parser = argparse.ArgumentParser('Parse VCF HEADER from FreeBayes,VarScan,GATK to fix it.  Output is to stdout.')
parser.add_argument('-f', '--file', help="freebayes vcf output file name")
parser.add_argument('-v', '--variantcaller', help="variant caller: F = freebayes, G = GATK, V = Varscan")

opts = parser.parse_args()

read=open(opts.file)


for line in read:
	if opts.variantcaller == 'F':
		
		if line.startswith('##FORMAT=<ID=DPR'):
			riga=(line.split(','))
			riga[riga.index('Number=A')]='Number=G'
			line=','.join(riga)
			
		# if line.startswith('##FORMAT=<ID=GQ'):
		# 	riga=(line.split(','))
		# 	riga[riga.index('Type=Integer')]='Type=Float'
		# 	line=','.join(riga)
	# 	if line.startswith('chr'):
	# 		riga=line.split('\t')
	# 		replaced = re.sub('[(chr)]', '', riga[0])
	# 		riga[0]=replaced
	# 		line='\t'.join(riga)
	
	elif opts.variantcaller == 'G':
	
		if line.startswith('##FORMAT=<ID=AD'):
			riga=(line.split(','))
			try:
				riga[riga.index('Number=.')]='Number=G'
				line=','.join(riga)
			except:
				pass

	elif opts.variantcaller == 'V':
	
		if line.startswith('chr'):
			riga=line.split('\t')
			ref=riga[3]
			alt=riga[4]
			if '/+' in ref:
				riga[3]=alt
				riga[4]=alt+ref.split('/+')[1]
				riga[7]+=";FIX"
			elif '/' in ref:
				riga[3]=ref.split('/')[0]
				riga[7]+=";FIX"
				
			if '/-' in alt:
				riga[4]=ref
				riga[3]=alt.split('/-')[0]
				riga[7]+=";FIX"
			elif '/' in alt:
				riga[4]=alt.split('/')[0]
				riga[7]+=";FIX"

			line='\t'.join(riga)
	print line.rstrip()
