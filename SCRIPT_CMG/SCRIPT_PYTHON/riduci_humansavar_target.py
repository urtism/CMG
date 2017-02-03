import argparse
import os

def main():

	parser = argparse.ArgumentParser('add chr to phylop and cons db. OUTPUT TO STDOUT')
	parser.add_argument('-f', '--infile', help="HUMANSAVAR path")
	parser.add_argument('-g', '--genes', help="list of genes")


	global opts 
	opts = parser.parse_args()

	genefile=open(opts.genes,'r')
	
	genelist = genefile.readlines()
	#print genelist


	infile=open(opts.infile,'r')
	for line in infile:
		line = line.rstrip()
		if line.startswith('Main_'):
				print line

		elif line.split('\t')[0]+'\r\n' in genelist :
			print line

	infile.close()



main ()