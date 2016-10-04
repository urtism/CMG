import argparse
import sys
from tr import tr
import re


parser = argparse.ArgumentParser('Parse VCF output from FreeBayes to output valid VCF.  Output is to stdout.')
parser.add_argument('-f', '--file', help="freebayes vcf output file name")
opts = parser.parse_args()

read=open(opts.file)

for line in read:
	if line.startswith('##FORMAT=<ID=DPR'):
		riga=(line.split(','))
		riga[riga.index('Number=A')]='Number=G'
		line=','.join(riga)
	print line.rstrip()