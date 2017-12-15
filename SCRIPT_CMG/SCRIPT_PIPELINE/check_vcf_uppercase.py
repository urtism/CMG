import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser('This script checks if ref and alt are uppercase')
	parser.add_argument('-v','--vcf', help="vcf to check")
	parser.add_argument('-o','--out', help="out file")

	global opts
	opts = parser.parse_args()
	
	vcf=open(opts.vcf,'r')
	out=open(opts.out,'w')

	for line in vcf:
		if line.startswith('#'):
			out.write(line)
		else:
			line_split=line.split('\t')
			line_split[3]=line_split[3].upper()
			line_split[4]=line_split[4].upper()
			out.write('\t'.join(line_split))