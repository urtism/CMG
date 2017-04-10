import argparse
import re
import string
from Bio.Seq import Seq
from Bio import SeqIO, Entrez
from Bio import SeqUtils
from pyfasta import Fasta
import pysam



def Estrai_feature_bam(CHR,POS,REF,ALT,bam1,bam2):

	bamfile = pysam.AlignmentFile(bam1, "rb")
	bamfile2 = pysam.AlignmentFile(bam2, "rb")

	pair = 0
	unpair = 0
	unmapped = 0
	tot = 0
	qc_fail = 0
	is_proper_pair_mapped = 0
	is_read1 = 0
	is_read1_reverse = 0
	is_read2 = 0
	is_read2_reverse = 0
	is_reverse = 0
	is_secondary = 0
	is_supplementary = 0
	is_read1_forward = 0
	is_read2_forward = 0
	is_unmapped = 0
	coverage = 0
	mate_is_unmapped = 0
	mapping_quality = 0


	for read in bamfile.fetch(CHR, int(POS)-1, int(POS), until_eof=True):
		print read
		tot += 1
		mapping_quality += read.mapping_quality
		if read.is_unmapped:
			unmapped += 1
			
		if read.is_paired:
			pair += 1

		if not read.is_paired:
			print read
			unpair += 1

		if read.is_proper_pair:
			is_proper_pair_mapped += 1
			
		if read.is_unmapped:
			unmapped += 1
			
		if read.is_qcfail:
			qc_fail += 1
			
		if read.is_read1:
		 	is_read1 += 1
		 	if read.is_reverse:
		 		is_read1_reverse += 1
		 	if not read.is_reverse:
		 		is_read1_forward += 1
		 	
		if read.is_read2:
		 	is_read2 += 1
		 	if read.is_reverse:
		 		is_read2_reverse += 1
		 	if not read.is_reverse:
		 		is_read2_forward += 1
		 	
		if read.is_reverse:
		 	is_reverse += 1
		 	
		if read.is_secondary:
		 	is_secondary += 1
		 	
		if read.is_supplementary:
		 	is_supplementary += 1
		 	
		if read.mate_is_unmapped:
		 	mate_is_unmapped += 1
	

	if tot != 0:
		mean_mapping_qual = mapping_quality/tot
	else:
		mean_mapping_qual = 0

	coverage = bamfile.count(CHR, int(POS)-1, int(POS),until_eof=True)

	print str(CHR) + ' ' + ' ' + str(POS) + ' ' + str(REF) + ' ' + str(ALT)
	print 'Pair: ' + str(pair)
	print 'Unpair: ' + str(unpair)
	print 'Unmapped: ' + str(unmapped)
	print 'qc_fail: ' + str(qc_fail)
	print 'is_read1: ' + str(is_read1)
	print 'is_read1_reverse: ' + str(is_read1_reverse)
	print 'is_read1_forward: ' + str(is_read1_forward)
	print 'is_read2: ' + str(is_read2)
	print 'is_read2_reverse: ' + str(is_read2_reverse)
	print 'is_read2_forward: ' + str(is_read2_forward)
	print 'is_reverse: ' + str(is_reverse)
	print 'is_secondary: ' + str(is_secondary)
	print 'is_supplementary: ' + str(is_supplementary)
	print 'is_proper_pair_mapped: ' + str(is_proper_pair_mapped)
	print 'Coverage: ' + str(coverage)
	print 'mate_is_unmapped: ' + str(mate_is_unmapped)
	print 'mean_mapping_qual: ' + str(mean_mapping_qual)
	print 'is_unmapped: ' + str(is_unmapped)
	print 'Mapping_qual: ' + str(mapping_quality)
	print 'Tot: ' + str(tot) + '\n\n'

	bamfile.close()

#-------------------------BAM2-------------------



	pair = 0
	unpair = 0
	unmapped = 0
	tot = 0
	qc_fail = 0
	is_proper_pair_mapped = 0
	is_read1 = 0
	is_read1_reverse = 0
	is_read2 = 0
	is_read2_reverse = 0
	is_reverse = 0
	is_secondary = 0
	is_supplementary = 0
	is_read1_forward = 0
	is_read2_forward = 0
	is_unmapped = 0
	coverage = 0
	mate_is_unmapped = 0
	mapping_quality = 0

	for read2 in bamfile2.fetch(CHR, int(POS)-1, int(POS), until_eof=True):
		print read2
		tot += 1
		mapping_quality += read.mapping_quality
		if read.is_unmapped:
			unmapped += 1
			
		if read.is_paired:
			pair += 1

		if not read.is_paired:
			print read
			unpair += 1

		if read.is_proper_pair:
			is_proper_pair_mapped += 1
			
		if read.is_unmapped:
			unmapped += 1
			
		if read.is_qcfail:
			qc_fail += 1
			
		if read.is_read1:
		 	is_read1 += 1
		 	if read.is_reverse:
		 		is_read1_reverse += 1
		 	if not read.is_reverse:
		 		is_read1_forward += 1
		 	
		if read.is_read2:
		 	is_read2 += 1
		 	if read.is_reverse:
		 		is_read2_reverse += 1
		 	if not read.is_reverse:
		 		is_read2_forward += 1
		 	
		if read.is_reverse:
		 	is_reverse += 1
		 	
		if read.is_secondary:
		 	is_secondary += 1
		 	
		if read.is_supplementary:
		 	is_supplementary += 1
		 	
		if read.mate_is_unmapped:
		 	mate_is_unmapped += 1
			

	if tot != 0:
		mean_mapping_qual = mapping_quality/tot
	else:
		mean_mapping_qual = 0

	coverage = bamfile2.count(CHR, int(POS)-1, int(POS),until_eof=True)

	print str(CHR) + ' ' + ' ' + str(POS) + ' ' + str(REF) + ' ' + str(ALT)
	print 'Pair: ' + str(pair)
	print 'Unpair: ' + str(unpair)
	print 'Unmapped: ' + str(unmapped)
	print 'qc_fail: ' + str(qc_fail)
	print 'is_read1: ' + str(is_read1)
	print 'is_read1_reverse: ' + str(is_read1_reverse)
	print 'is_read1_forward: ' + str(is_read1_forward)
	print 'is_read2: ' + str(is_read2)
	print 'is_read2_reverse: ' + str(is_read2_reverse)
	print 'is_read2_forward: ' + str(is_read2_forward)
	print 'is_reverse: ' + str(is_reverse)
	print 'is_secondary: ' + str(is_secondary)
	print 'is_supplementary: ' + str(is_supplementary)
	print 'is_proper_pair_mapped: ' + str(is_proper_pair_mapped)
	print 'Coverage: ' + str(coverage)
	print 'mate_is_unmapped: ' + str(mate_is_unmapped)
	print 'mean_mapping_qual: ' + str(mean_mapping_qual)
	print 'is_unmapped: ' + str(is_unmapped)
	print 'Mapping_qual: ' + str(mapping_quality)
	print 'Tot: ' + str(tot) + '\n\n'

	bamfile2.close()






def main():

	parser = argparse.ArgumentParser('Con questo script vado ad eliminare le righe che sono annotate male a causa del sito multiallelico. Usato su esperimento somic')
	parser.add_argument('-i','--input',help="bam file")
	parser.add_argument('-b','--bam',help="bam da confrontare")
#	parser.add_argument('-o','--outfile',help="file di out con i risultati")

	global opts
	opts = parser.parse_args()
#	out = open(opts.outfile, 'w')
	bam1 = opts.input
	bam2 = opts.bam

	CHR = 'chr1'
	POS = 3311101
	REF = 'C'
	ALT = 'G'


	Estrai_feature_bam(CHR,POS,REF,ALT,bam1,bam2)

main()