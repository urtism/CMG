import argparse
import re
import string
from Bio.Seq import Seq
from Bio import SeqIO, Entrez
from Bio import SeqUtils
from pyfasta import Fasta
import pysam



def Find_Repetiotions(CHR,POS,REF,ALT,fasta,bam):

	ref = Fasta(fasta)

	sorted(ref.keys())
	converted_POS = int(POS)-1
	start = int(converted_POS)-5
	stop = int(converted_POS)+6
	start_perc = int(converted_POS)-10
	stop_perc = int(converted_POS)+11
	
	sequence = ref[CHR][start:stop]

	seq = sequence.upper()

	nucleotides = ['A','T','C','G']

	Homopolymer = 'Homopolymer=false'
	Repeat = 'RepSeq=0'

	for letter in nucleotides:

		if letter*5 in seq:

			start_hom = seq.index(letter*5)

			if start_hom <= 5 or start_hom == 7:
				#Controllo se ci sia qualche lettere minuscola che indicano sequenze ripetute nel genoma
				if any(c.islower() for c in sequence):

					Homopolymer = 'Homopolymer=true'
					Repeat = 'RepSeq=1'
					#print CHR + ' ' + POS + ' ' + REF + ' ' + ALT + ' ' + ref[CHR][start:converted_POS] + ' ' + ref[CHR][converted_POS:stop] + ' ' + Repeat + ' ' + Homopolymer
					break

				elif sequence.isupper():

					Homopolymer = 'Homopolymer=true'
					Repeat = 'RepSeq=0'
					#print CHR + ' ' + POS + ' ' + REF + ' ' + ALT + ' ' + ref[CHR][start:converted_POS] + ' ' + ref[CHR][converted_POS:stop] + ' ' + Repeat + ' ' + Homopolymer
					break

				else:

					Homopolymer = 'Homopolymer=true'
					Repeat = 'RepSeq=.'

			else:
				continue

		else:
			continue

	#Ora estraggo i q-score
	bamfile = pysam.AlignmentFile(bam, "rb")

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

	#for read in bamfile.fetch():
	#	if read.is_unmapped:
	#	 	is_unmapped += 1
	#	 	print read
		 	#if read.is_reverse:
		 	#	is_read1_reverse += 1
		 	#else:
		 	#	non_rev_read1 += 0


	for read in bamfile.fetch(CHR, int(POS)-1, int(POS), until_eof=True):
		#print str(read.reference_id) + '\n' + str(read.next_reference_id)

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
			
	#unpaired = (is_read1+is_read2) - pair

	if tot != 0:
		mean_mapping_qual = mapping_quality/tot
	else:
		mean_mapping_qual = 0


	#for pileupcolumn in bamfile.pileup(CHR, int(POS)-1, int(POS)+1):
	#	coverage = pileupcolumn.n

	coverage = bamfile.count(CHR, int(POS)-1, int(POS),until_eof=True)
	
	for read in bamfile.fetch(CHR, 3311100, 3311101, until_eof=True):
			print read
	print CHR + ' ' + ' ' + POS + ' ' + REF + ' ' + ALT + ' ' + ref[CHR][start:converted_POS] + ' ' + ref[CHR][converted_POS:stop] + ' ' + Repeat + ' ' + Homopolymer
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
	print 'Mapping_qual' + str(mapping_quality)
	print 'Tot: ' + str(tot)


	#Perc_pair = float(pair/tot)
	#Perc_unpair = float(unpair/tot)
	#print CHR + ' ' + POS + ' ' + REF + ' ' + ALT + ' ' + ref[CHR][start:converted_POS] + ' ' + ref[CHR][converted_POS:stop] + ' ' + Repeat + ' ' + Homopolymer
	#print Perc_unpair
	#print Perc_pair

		#print read.query_qualities[]

	# for pileupcolumn in samfile.pileup(CHR, start_perc, stop_perc):
	# 	q = pileupcolumn.query_qualities
	# 	q.query_sequence = q.query_sequence[start_perc:stop_perc]
	# 	q.query_qualities = q[stop_perc:stop_perc]
	# 	print q



	bamfile.close()

	return Homopolymer
			#continue



			#print CHR + ' ' + POS + ' ' + REF + ' ' + ALT + ' ' + ref[CHR][start:converted_POS] + ' ' + ref[CHR][converted_POS] + ' ' + ref[CHR][converted_POS+1:stop]
			
			#if converted_POS = 
			#print converted_POS 






def main():

	parser = argparse.ArgumentParser('Con questo script vado ad eliminare le righe che sono annotate male a causa del sito multiallelico. Usato su esperimento somic')
	parser.add_argument('-i','--input',help="vcf contenente chrOM e POS")
	parser.add_argument('-r','--reference',help="reference formato fasta")
	parser.add_argument('-b','--bam',help="vcf contenente chrOM e POS")
	parser.add_argument('-o','--outfile',help="file vcf annotato in INFO con una tag del tipo Repeats=Tandem/Homopolymer/null")

	global opts
	opts = parser.parse_args()
	vcf = open(opts.input, 'r')
	out = open(opts.outfile, 'w')
	fasta = opts.reference
	bam = opts.bam
	save_head = []
	count=0

	for line in vcf:

		#Vado a lavorare sull'header per estrarre le informazioni dell'header e salvarle nel nuovo file
		if line.startswith('#'):
			line = line.rstrip()

			if line.startswith('##'):
				save_head.append(line)

			elif line.startswith('#chrOM'):
				save_head.append(line)
				out.write('\n'.join(save_head) + '\n')

			else:
				save_head.append(line)

		else:
			line = line.rstrip().split('\t')
			CHR=line[0]
			POS=line[1]
			REF=line[3]
			ALT=line[4]

			Homopolymer = Find_Repetiotions(CHR,POS,REF,ALT,fasta,bam)
		count += 1

	vcf.close()
	out.close()

main()