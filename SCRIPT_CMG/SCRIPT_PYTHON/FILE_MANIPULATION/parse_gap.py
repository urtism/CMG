import argparse
import re
import string
import os
import math


def intersect_gap(filename):
	
	file = open(opts.directory + '/' + filename,'r')
	gap = []


	for line in file:

		line = line.rstrip()
		line = line.split(',')

		if line[0].startswith('#'):
 	 		gap.append(line)
			continue

		else:		
			target = open(opts.file,'r') 	
			
			for tar in target:

		 		tar = tar.rstrip()
				tar = tar.split('\t')
					
		 	 	if tar[0].startswith('@'):
		 	 		continue

		 	 	elif line[0] == tar[0] and line[1] >= tar[1] and line[2] <= tar[2]:
		 	 		#print line[0]
		 	 		#print tar[0]
		 	 		#print tar[1]
		 	 		#print line[1]
		 	 		#print tar[2]
		 	 		#print line[2]
		 	 		print line[3]
		 	 		line[3] = tar[4]
		 	 		print line[3]
		 	 		gap.append(line)
		 	 		#print line
			target.close()

	return gap

	file.close()

def main():

	parser = argparse.ArgumentParser('Estrai solo i gap presenti nel file bed o list fornito')	
	parser.add_argument('-d','--directory',help="directory contenente i gap")
	parser.add_argument('-f','--file',help="file contenente il target da cui fare il parse")
	parser.add_argument('-o','--outdir',help="directory di output")

	global opts

	opts = parser.parse_args()

	for filename in os.listdir(opts.directory):

		out = open(opts.outdir + '/' + filename,'w')
		
		gep = intersect_gap(filename)

		for elem in gep:
			print elem
			out.write(','.join(elem) + '\n')

		out.close()


main()