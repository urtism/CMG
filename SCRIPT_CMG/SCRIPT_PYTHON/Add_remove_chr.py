import argparse


def main():

	parser = argparse.ArgumentParser('add or remove chr in vcf file or bedfile. OUT IN STDOUT')
	parser.add_argument('-f', '--infile', help="vcf")

	global opts 
	opts = parser.parse_args()


	infile=open(opts.infile,'r')

	for line in infile:
		line=line.rstrip()
		if line.startswith('#'):
			print line
		else:
			line_split=line.split('\t')
			chr=line_split[0]
			if chr.startswith('chr'):
				print '\t'.join([chr.split('chr')[1]] + line_split[1:])
			else:
				print "chr"+line
		


main ()
