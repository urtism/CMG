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
		if line.startswith('##FORMAT=<ID=GQ'):
			riga=(line.split(','))
			riga[riga.index('Type=Integer')]='Type=Float'
			line=','.join(riga)
	# 	if line.startswith('chr'):
	# 		riga=line.split('\t')
	# 		replaced = re.sub('[(chr)]', '', riga[0])
	# 		riga[0]=replaced
	# 		line='\t'.join(riga)
	
	elif opts.variantcaller == 'G':
	
		if line.startswith('##FORMAT=<ID=AD'):
			riga=(line.split(','))
			riga[riga.index('Number=.')]='Number=G'
			line=','.join(riga)
	
	print line.rstrip()