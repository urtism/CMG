import pysam
import argparse

if __name__ == '__main__':

	parser = argparse.ArgumentParser('Bed to vcf')
	parser.add_argument('-b','--bed', help="in.bed")
	parser.add_argument('-r','--ref', help="reference.fasta")
	parser.add_argument('-o','--out', help="out.vcf")

	global opts
	opts = parser.parse_args()
	fa = pysam.Fastafile(opts.ref)
	vcf=open(opts.out,'w')
	vcf.write('##fileformat=VCFv4.1'+'\n')
	vcf.write('##INFO=<ID=ID,Number=1,Type=String,Description="Variant ID">'+'\n')
	vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"'+'\n')
	with open(opts.bed,'r') as bed:
		for line in bed:
			if line.startswith('#') or line.startswith('@'):
				continue
			else:
				chrom,start,end,ref,alt=line.rstrip().split('\t')[:5]
				id=line.rstrip().split('\t')[-1]

				try:
					fa.fetch(chrom,start)
				except:
					chrom='chr'+chrom

				if ref == '-':
					ref = fa.fetch(chrom, int(start)-1, int(start)).upper()

					alt = ref.upper() + alt.upper()
				elif alt== '-':
					ref = fa.fetch(chrom, int(start)-1, int(end)).upper()
					alt = ref[0].upper()

				vcf.write('\t'.join([chrom, start, '.', ref, alt,'.','.','ID='+id,'.','.'])+'\n')