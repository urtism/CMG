#!/usr/bin/python

import sys

#bed='/home/jarvis/NGS_ANALYSIS/TARGET/trusight_cardio_manifest_a.bed'
bed=sys.argv[1]

for line in open(bed,'r'):
	line=line.rstrip()
	chr=line.split('\t')[0]
	start=line.split('\t')[1]
	stop=line.split('\t')[2]

	start=int(start)-150
	stop=int(stop)+150

	print '\t'.join([chr,str(start),str(stop)])
