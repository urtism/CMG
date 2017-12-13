import string
import argparse

def main():

	parser = argparse.ArgumentParser('\n\nQuesto tool converte solo REF e ALT nel vcf in UPPERCASE per il vcf di BCFTOOLS\n')
	parser.add_argument('-I','--input',help="Input file in formato .vcf")
	parser.add_argument('-O','--outfile',help="File di output in vcf format")

	global opts
	opts = parser.parse_args()
	out = open(opts.outfile, 'w')
	Variant_list = {}
	Variant_table = {}

	with open(opts.input,'r') as vcf:

		for line in vcf:

			if line.startswith('#'):

				out.write(line)

			else:

				line = line.rstrip()
				line = line.split('\t')
				line[3]=line[3].upper()
				line[4]=line[4].upper()

				out.write('\t'.join(line) + '\n')
				

	out.close()

main()