from Bio.SeqUtils import GC
import pysam
import argparse


if __name__ == '__main__':

	parser = argparse.ArgumentParser('Calculate GC content from intervals in bed file. Output is a new bed file with GC content in percentage. Default output is stdout, use -o or --out to choose outpath ')
	parser.add_argument('-b', '--bed', help="Bed file ", required=True)
	parser.add_argument('-r', '--ref', help="Reference.fasta", required=True)
	parser.add_argument('-o', '--out', help="Bed file in output",default=None)

	global opts
	opts = parser.parse_args()

	fasta=pysam.FastaFile(opts.ref)
	bed=open(opts.bed,'r')
	out=open(opts.out,'w')

	for line in bed:
		
		line=line.rstrip('\n').split('\t')
		
		chr,start,stop = line[:3]
		
		seq = fasta.fetch( reference=chr, start=int(start), end=int(stop), region=None)
		GCcontent = str(GC(seq))
		print '\t'.join(line).rstrip()
		out.write('\t'.join(['\t'.join(line).rstrip(),GCcontent]) +'\n')
	
	bed.close()
	out.close()