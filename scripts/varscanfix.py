import argparse
import sys
from tr import tr
import re


parser = argparse.ArgumentParser('Parse VCF output from Varscan to output valid VCF.  Output is to stdout.')
parser.add_argument('-f', '--file', help="varscan vcf output file name")
opts = parser.parse_args()

read=open(opts.file)

for line in read:
	parts=line.split('\t')
	if line.startswith('#'):
		print line.rstrip()
	else:
		ref=parts[3]
		alt=parts[4]
		format=parts[8].split(':')
		samples=parts[9:]
		

		if len(ref)==2:
			r1=re.compile('\+[ACTG]+')
			r2=re.compile('[ACTG]+')

			if r1.match(ref[1]):  

			# AG/+GG	A ----> A	GA,AGGG
				refer=ref[0]
				a=str(ref[0])+str(ref[1][1:])
				alt=alt+[a]
			elif r2.match(ref[1]):

			# AG/GG	A ----> AGG	A,AG
				if ref[0][1:] in ref[1]:
					alt=alt+[ref[0][0]+ref[1][len(ref[0][1:]):]]
					a=str(ref[0][0])+str(ref[1])
					refer=a
				else:
					print "le due ref non sono simili----> bisogna gestirlo???"

		elif len(alt)==2:
			r1=re.compile('\-[ACTG]+')
			r2=re.compile('[ACTG]+')
			r3=re.compile('[ACTG]')

			if r1.match(alt[1]):

			# A	AG/-G ----> AG	A,AGG
				a=str(ref[0])+str(alt[1][1:])
				refer=a
				alt=[str(alt[0])+str(alt[1][1:])]+ref

			elif r3.match(alt[1]):
				refer=ref[0]
				
			elif r2.match(alt[1]):

			# A	AG/GG ----> A	AG,AGG
				if alt[0][1:] in alt[1]:
					alt[1]=str(ref[0])+str(alt[1])
					refer=ref[0]
				else:
					print "le due alt non sono simili----> bisogna gestirlo???"

			


		normal[format.index('FREQ')]=str(float((normal[format.index('FREQ')]).split('%')[0])/100)
		tumor[format.index('FREQ')]=str(float((tumor[format.index('FREQ')]).split('%')[0])/100)

		parts[3]=refer
		parts[4]=','.join(alt)
		parts[9]=':'.join(normal)
		parts[10]=':'.join(tumor)

		print ('\t'.join(parts)).rstrip()

