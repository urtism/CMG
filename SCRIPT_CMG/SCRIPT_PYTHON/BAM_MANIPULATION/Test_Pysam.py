import pysam

files = ['/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_01_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_02_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_03_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_04_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_05_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_06_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_07_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_08_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_09_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_10_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_11_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_12_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_13_Conn.bam',
		'/home/matteo/Scrivania/TEST_DATA/20171120_Run_109_Cardio_2.1_VCF_IEVA/20171120_14_Conn.bam']

# for sample in files:
# 	print sample
# 	bampile = pysam.AlignmentFile(sample, "rb")
# 	for pileupcolumn in bampile.pileup(stepper='nofilter'):
# 		for pileupread in pileupcolumn.pileups:
			#print pileupread
			#if pileupread.alignment.is_duplicate:
			#	print 'DUPLICATE:\t'
			#if pileupread.alignment.is_unmapped:
			#	print 'UNMAPPED:\t'
			#if pileupread.alignment.mapping_quality == 0:
			#	print 'MQ0:\t'


			#if pileupread.alignment.is_paired == False:
			#	print 'NOT PAIRED:\t'
			#if pileupread.alignment.is_proper_pair == False:
			#	print 'NOT PROPER PAIRED:\t'

			# if pileupread.alignment.is_supplementary:
			# 	print 'SUPPLEMENTARY:\t'

			# if pileupread.alignment.is_secondary:
			# 	print 'SECONDARY:\t'