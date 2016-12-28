import pysam

bamfile = pysam.AlignmentFile("/home/jarvis/git/CMG/TEST_INPUT/SOURCE_FILE/20161213_01_Conn.bam","rb")

for read in bamfile.fetch('chr1'):
     print read
bamfile.close()