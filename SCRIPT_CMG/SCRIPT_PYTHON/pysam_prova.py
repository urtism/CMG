import pysam
import math
import statistics

bam = pysam.AlignmentFile("/home/minime/Scrivania/TEST/20161213_02_Conn.bam", "rb")
fasta = pysam.FastaFile("/home/minime/NGS_TOOLS/hg19/ucsc.hg19.fasta")



chrom='chr2'
start=21225013
stop=21225014


#print bam.count(reference=chrom, start=start, end=stop, until_eof=False, read_callback='nofilter')
print bam.count_coverage(reference=chrom, start=start, end=stop)
#print bam.parse_region(reference=chrom, start=start, end=stop, tid=None)


# for pc in bam.pileup(reference=chrom, start=start, end=stop):
# 	for reads in pc.pileups:
# 		print reads

QB=[]
MQ=[]
BQ=[]

for pileupcolumn in bam.pileup(reference=chrom, start=start, end=stop):
	if pileupcolumn.reference_pos>=start and pileupcolumn.reference_pos<stop:
		for pileupread in pileupcolumn.pileups:
			QB+=[pileupread.alignment.query_qualities[pileupread.query_position]]
			MQ+=[pileupread.alignment.mapping_quality]
			BQ+=[pileupread.alignment.qual[pileupread.query_position]]
			pos=pileupread.alignment.get_reference_positions()[pileupread.query_position]
			base=pileupread.alignment.query_sequence[pileupread.query_position]
			dp=pileupcolumn.n


		print ("pos=%s, base=%s, dp=%s, BQ=%s, MQ=%s" % (pos, base, dp, str(statistics.mean(QB)),MQ[pileupread.query_position]))
		print BQ
			# print ("pos=%s, base=%s, dp=%s, BQ=%s" % (pileupread.alignment.get_reference_positions()[pileupread.query_position],pileupread.alignment.query_sequence[pileupread.query_position],
			# 		pileupcolumn.n, pileupread.alignment.query_qualities[pileupread.query_position]))


bam.close


samfile = pysam.AlignmentFile("/home/minime/Scrivania/TEST/20161213_02_Conn.bam", "rb" )
for pileupcolumn in samfile.pileup(chrom, start, stop):
    print ("\ncoverage at base %s = %s" %
           (pileupcolumn.pos, pileupcolumn.n))
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            # query position is None if is_del or is_refskip is set.
            print ('\tbase in read %s = %s' %
                  (pileupread.alignment.query_name,
                   pileupread.alignment.query_sequence[pileupread.query_position]))

samfile.close()
