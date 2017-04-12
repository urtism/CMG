import argparse
import re
import string
import sys
from Bio.Seq import Seq
from Bio import SeqIO, Entrez
from Bio import SeqUtils
from pyfasta import Fasta
import pysam


#--------------------------------------------------- -FUNCTIONS-----------------------------------------------

def SimpleRepeats_Finder(Reference, sample_list, variant,H_CHR):

	#Leggo il file fasta reference
	Ref = Fasta(Reference)

	#Estraggo Chrom, Pos, Ref ed Alt
	sorted(Ref.keys())
	CHR = variant[H_CHR.index('#CHROM')]
	POS = int(variant[H_CHR.index('POS')])-1
	REF = variant[H_CHR.index('REF')]
	ALT = variant[H_CHR.index('ALT')]
	
	#Estraggo la sequenza dalla reference corrispondente alla variante trovata
	sequence = Ref[CHR][POS-20:POS+21]

	#Trasformo la sequenza in lettere maiuscole
	seq = sequence.upper()

	#Definisco i pattern da andare a trovare nella sequenza
	Homopolymer_pattern = re.compile(r'(A){5,}|(G){5,}|(C){5,}|(T){5,}')

	#Estraggo gli omopolimeri trovati in una lista contenente tuple. Es: [('AAAAAA',4,6),('TTTTTTTTTT',6,10)] in cui il primo elemento indica 
	#il pattern di omopolimeri, il secondo l'indice in cui inizia ed il terzo la lunghezza del pattern
	Homopolymer_list = [(m.group(), m.start(), len(m.group())) for m in re.finditer(Homopolymer_pattern, seq)]

	Simple_Repeat = '.'
	Simple_Repeat_Lenght = '.'

	#Vado a vedere tutti gli omopolimeri identificati nella sequenza. La variante deve cadere nell'omopolimero o iniziare in posizione
	#immediatamente precedente a quella dell'omopolimero. Per cui prendo l'indice di start dell'hom, lo sommo alla lunghezza ed applico il
	#seguente condizionale
	for poly in Homopolymer_list:

		Start_Hom = poly[1]
		End_Hom = poly[1]+poly[2]-1

		if poly[1] == 21:

			Simple_Repeat = 'Homopolymer'
			Simple_Repeat_Lenght = poly[2]
			break
		
		elif poly[1] <= 20 and End_Hom >= 20:

			Simple_Repeat = 'Homopolymer'
			Simple_Repeat_Lenght = poly[2]
			break

		else:
			continue

	#Passo ora all'identificazione delle Simple Repeats come come poly(AT|ATT|ATAA|etc...)
	if '.' in Simple_Repeat:

		Simple_Repeat_pattern = re.compile(r'(AC){4,}|(AG){4,}|(AT){4,}|(CA){4,}|(CG){4,}|(CT){4,}|(GA){4,}|(GC){4,}|'
											'(GT){4,}|(TA){4,}|(TC){4,}|(TG){4,}|(ACC){3,}|(AGG){3,}|(ATT){3,}|(CAA){3,}|'
											'(CGG){3,}|(CTT){3,}|(GAA){3,}|(GCC){3,}|(GTT){3,}|(TAA){3,}|(TCC){3,}|(TGG){3,}|'
											'(CCA){3,}|(GGA){3,}|(TTA){3,}|(AAC){3,}|(GGC){3,}|(TTC){3,}|(AAG){3,}|(CCG){3,}|'
											'(TTG){3,}|(AAT){3,}|(CCT){3,}|(GGT){3,}|(CAC){3,}|(GAG){3,}|(TAT){3,}|(ACA){3,}|'
											'(GCG){3,}|(TCT){3,}|(AGA){3,}|(CGC){3,}|(TGT){3,}|(ATA){3,}|(CTC){3,}|(GTG){3,}|'
											'(ACCC){3,}|(AGGG){3,}|(ATTT){3,}|(CAAA){3,}|(CGGG){3,}|(CTTT){3,}|(GAAA){3,}|(GCCC){3,}|'
											'(GTTT){3,}|(TAAA){3,}|(TCCC){3,}|(TGGG){3,}|(CACC){3,}|(GAGG){3,}|(TATT){3,}|(ACAA){3,}|'
											'(GCGG){3,}|(TCTT){3,}|(AGAA){3,}|(CGCC){3,}|(TGTT){3,}|(ATAA){3,}|(CTCC){3,}|(GTGG){3,}|'
											'(CCAC){3,}|(GGAG){3,}|(TTAT){3,}|(AACA){3,}|(GGCG){3,}|(TTCT){3,}|(AAGA){3,}|(CCGC){3,}|'
											'(TTGT){3,}|(AATA){3,}|(CCTC){3,}|(GGTG){3,}|(CCCA){3,}|(GGGA){3,}|(TTTA){3,}|(AAAC){3,}|'
											'(GGGC){3,}|(TTTC){3,}|(AAAG){3,}|(CCCG){3,}|(TTTG){3,}|(AAAT){3,}|(CCCT){3,}|(GGGT){3,}')

		Simple_Repeat_list = [(m.group(), m.start(), len(m.group())) for m in re.finditer(Simple_Repeat_pattern, seq)]

		for pattern in Simple_Repeat_list:

			Start_Hom = poly[1]
			End_Hom = poly[1]+poly[2]-1

			if poly[1] == 21:

				Simple_Repeat = 'Homopolymer'
				Simple_Repeat_Lenght = poly[2]
				break
			
			elif poly[1] <= 20 and End_Hom >= 20:

				Simple_Repeat = 'Homopolymer'
				Simple_Repeat_Lenght = poly[2]
				break

			else:
				continue


		#print str(Simple_Repeat_list) + ' ' + str(sequence) + ' ' + str(CHR) + ' ' + str(POS+1) + ' ' + str(REF) + ' ' + str(ALT)

	#print str(Homopolymer_list) + ' ' + sequence + ' ' + str(CHR) + ' ' + str(POS+1) + ' ' + str(REF) + ' ' + str(ALT) + ' ' + str(Simple_Repeat) + ' ' + str(Simple_Repeat_Lenght)#+ str(Start_Hom) + ' ' + str(End_Hom)

	#print str(CHR) + ' ' + str(POS+1) + ' ' + str(REF) + ' ' + str(ALT) + ' ' + Ref[CHR][POS-15:POS] + ' ' + Ref[CHR][POS] + ' ' + Ref[CHR][POS+1:POS+16] + ' ' + str(Homopolymer) + '\n'









		#Simple_Repeat_pattern = re.compile(r'(AC){4,}|(AG){4,}|(AT){4,}|(CA){4,}|(CG){4,}|(CT){4,}|(GA){4,}|(GC){4,}|'
		#									'(GT){4,}|(TA){4,}|(TC){4,}|(TG){4,}')


		# Simple_Repeat_pattern = re.compile(r'(ACC){3,}|(AGG){3,}|(ATT){3,}|(CAA){3,}|(CGG){3,}|(CTT){3,}|(GAA){3,}|'
		# 									'(GCC){3,}|(GTT){3,}|(TAA){3,}|(TCC){3,}|(TGG){3,}|(CCA){3,}|(GGA){3,}|'
		# 									'(TTA){3,}|(AAC){3,}|(GGC){3,}|(TTC){3,}|(AAG){3,}|(CCG){3,}|(TTG){3,}|'
		# 									'(AAT){3,}|(CCT){3,}|(GGT){3,}|(CAC){3,}|(GAG){3,}|(TAT){3,}|(ACA){3,}|'
		# 									'(GCG){3,}|(TCT){3,}|(AGA){3,}|(CGC){3,}|(TGT){3,}|(ATA){3,}|(CTC){3,}|'
		# 									'(GTG){3,}')

		# Simple_Repeat_pattern = re.compile(r'(ACCC){3,}|(AGGG){3,}|(ATTT){3,}|(CAAA){3,}|(CGGG){3,}|(CTTT){3,}|(GAAA){3,}|(GCCC){3,}|'
		# 									'(GTTT){3,}|(TAAA){3,}|(TCCC){3,}|(TGGG){3,}|(CACC){3,}|(GAGG){3,}|(TATT){3,}|(ACAA){3,}|'
		# 									'(GCGG){3,}|(TCTT){3,}|(AGAA){3,}|(CGCC){3,}|(TGTT){3,}|(ATAA){3,}|(CTCC){3,}|(GTGG){3,}|'
		# 									'(CCAC){3,}|(GGAG){3,}|(TTAT){3,}|(AACA){3,}|(GGCG){3,}|(TTCT){3,}|(AAGA){3,}|(CCGC){3,}|'
		# 									'(TTGT){3,}|(AATA){3,}|(CCTC){3,}|(GGTG){3,}|(CCCA){3,}|(GGGA){3,}|(TTTA){3,}|(AAAC){3,}|'
		# 									'(GGGC){3,}|(TTTC){3,}|(AAAG){3,}|(CCCG){3,}|(TTTG){3,}|(AAAT){3,}|(CCCT){3,}|(GGGT){3,}')









#-------------------------------------------------------MAIN--------------------------------------------------


def main():

	parser = argparse.ArgumentParser('Tool to extract sequence features from fasta file and bam file')
	parser.add_argument('-I','--input',help="input file in vcf format")
	parser.add_argument('-R','--reference',help="Reference file in fasta format")
	parser.add_argument('-L','--list',help="file containing bam list. Include path")
	parser.add_argument('-O','--outfile',help="Output vcf file")
	parser.add_argument('-RSeq','--RepeatSeq',action="store_true",help="Enable Simple Repeat finder. Check if a variant fall in repeated sequence like Homopolymer, poly(ATA). Output in INFO field")
	parser.add_argument('-RLen','--RepeatLenght',action="store_true",help="Report length of repeated sequence like Homopolymer, poly(ATA) if variants fall in. Cannot be calculated without enabling -SRep option. Output in INFO field")
	parser.add_argument('-RM','--RepeatMasker',action="store_true",help="Check if variant is a repeated sequence found with Repeat Masker. Tag in INFO field")
	parser.add_argument('-GC','--gcContent',action="store_true",help="Enable GC content. Tag in INFO field")
	parser.add_argument('-SBR','--StrandBiasReads',action="store_true",help="Strand Baias based on read orientation. Tag in FORMAT field")
	parser.add_argument('-AS','--AlignmentScore',action="store_true",help="Mean Alignment Score of reads mapping position. Tag in FORMAT field")
	parser.add_argument('-PDR','--DuplicateReference',action="store_true",help="Percentage of REF reads marked as duplicates. Tag in FORMAT field")
	parser.add_argument('-PDA','--DuplicateAlternate',action="store_true",help="Percentage of ALT reads marked as duplicates. Tag in FORMAT field")
	parser.add_argument('-DDup','--DeltaDuplicate',action="store_true",help="Difference for percentage duplicates in REF and ALT (REF-ALT). Cannot be calculated without enabling -PDR and -PDA options. Tag in FORMAT field")
	parser.add_argument('-QR','--MeanRefQscore',action="store_true",help="Mean Q-score for REF reads. Tag in FORMAT field")
	parser.add_argument('-QA','--MeanAltQscore',action="store_true",help="Mean Q-score for ALT reads. Tag in FORMAT field")
	parser.add_argument('-Q30','--PercentageReadsQscoreFilter',action="store_true",help="Percentage of reads with mean Q-score < 30. Tag in FORMAT field")
	parser.add_argument('-NM','--PercentageNumberOfMismatch',action="store_true",help="Percentage of reads with edit distance > 4. Tag in FORMAT field")

	global opts

	opts = parser.parse_args()

	vcf = open(opts.input,'r')
	out = open(opts.outfile,'w')
 	Reference = opts.reference
	bam_list = opts.list
	save_head = []
	file_format = []
	H_FILTER = []
	H_INFO = []
	H_FORMAT = []
	H_CHR = []

	#Estraggo dall'heder i campi FILTER, INFO e FORMAT
	for line in vcf:

		line = line.rstrip()

		if line.startswith('##fileformat'):
			file_format += [line]

		elif line.startswith('##FILTER'):
			H_FILTER += [line]

		elif line.startswith('##INFO'):
			H_INFO += [line]

		elif line.startswith('##FORMAT'):
			H_FORMAT += [line]

		elif line.startswith('#CHROM'):
			H_CHR += [line]
			break

	vcf.close()

	#Da sys.argv estraggo le opzioni date nella command line. Per ognuna di queste, scrivo il campo da inserire nell'header
	for option in sys.argv:

		if '-RSeq' in option:
			H_INFO += ['##INFO=<ID=RSeq,Number=1,Type=String,Description="Variant falls into repeated sequence like Homopolymer or SimpleRepeats in interval +-20bp">']
		if '-RLen' in option:
			H_INFO += ['##INFO=<ID=RLen,Number=1,Type=Integer,Description="Length of Repeated Sequence found in RSeq">']
		elif '-RM' in option:
			H_INFO += ['##INFO=<ID=RM,Number=.,Type=Falg,Description="Variant falls into repeated sequence from RepeatMasker">']
		elif '-GC' in option:
			H_INFO += ['##INFO=<ID=GC,Number=1,Type=Float,Description="GC content in sequence (+-25 bp)">']
		elif '-SBR' in option:
			H_FORMAT += ['##FORMAT=<ID=SBR,Number=4,Type=Float,Description="Strand Baias based on Read orientation: F1,R1,F2,R2">']
		elif '-AS' in option:
			H_FORMAT += ['##FORMAT=<ID=AS,Number=1,Type=Float,Description="Mean Alignment Score of reads mapping position">']
		elif '-PDR' in option:
			H_FORMAT += ['##FORMAT=<ID=PDR,Number=1,Type=Float,Description="Percentage of REF reads marked as duplicates">']
		elif '-PDA' in option:
			H_FORMAT += ['##FORMAT=<ID=PercentageDuplicateAlternate,Number=1,Type=Float,Description="Percentage of ALT reads marked as duplicates">']
		elif '-DDup' in option:
			H_FORMAT += ['##FORMAT=<ID=DeltaDuplicate,Number=1,Type=Float,Description="Difference for percentage duplicates in REF and ALT (REF-ALT)">']
		elif '-QR' in option:
			H_FORMAT += ['##FORMAT=<ID=MeanRefQscore,Number=1,Type=Float,Description="Mean Q-score for REF reads">']
		elif '-QA' in option:
			H_FORMAT += ['##FORMAT=<ID=MeanAltQscore,Number=1,Type=Float,Description="Mean Q-score for ALT reads">']
		elif '-Q30' in option:
			H_FORMAT += ['##FORMAT=<ID=PercentageReadsQscoreFilter,Number=1,Type=Float,Description="Percentage of reads with mean Q-score < 30">']
		elif '-NM' in option:
			H_FORMAT += ['##FORMAT=<ID=PercentageNumberOfMismatch,Number=1,Type=Float,Description="Percentage of reads with edit distance > 4">']
		else:
			continue

	#Scrivo in output la nuova header
	out.write('\n'.join(file_format) + '\n' + '\n'.join(H_FILTER) + '\n' + '\n'.join(H_INFO) + '\n' + '\n'.join(H_FORMAT) + '\n' + '\n'.join(H_CHR) + '\n')

	#Inserisco il controllo nelle opzioni -DDup e -SLen. Per essere calcolate bisogna avere -PDA -PDR e -SRep
	if '-RLen' in sys.argv and '-RSeq' not in sys.argv:
		sys.exit('\n' + "Cannot report --SimpleRepeatLenght without enabling --SimpleRepeat" + '\n')
	if '-DDup' in sys.argv and '-PDA' not in sys.argv:
		sys.exit('\n' + "Cannot report --DeltaDuplicate without enabling --PercentageDuplicateAlternate and --PercentageDuplicateReference" + '\n')
	if '-DDup' in sys.argv and '-PDR' not in sys.argv:
		sys.exit('\n' + "Cannot report --DeltaDuplicate without enabling --PercentageDuplicateAlternate and --PercentageDuplicateReference" + '\n')


	vcf = open(opts.input,'r')

	#Splitto l'header del vcf
	H_CHR = H_CHR[0].split('\t')

	#Estraggo dall'header la lista dei campioni:
	sample_list = H_CHR[0].split('\t')[9:]

	#Ora passiamo ad annotare le varianti nel vcf:
	for variant in vcf:

		if variant.startswith('#'):

			continue

		else:

			variant = variant.rstrip()
			variant = variant.split('\t')

			if opts.SimpleRepeat:
				variant = SimpleRepeats_Finder(Reference, sample_list, variant, H_CHR)
				continue

			if opts.RepeatedSequence:
			#Repeat = 'RepSeq=0'
			#if any(c.islower() for c in sequence):
			#	Repeat = 'RepSeq=1'
				continue

			if opts.gcContent:
			#gc_content = SeqUtils.GC(seq)
			#print gc_content
				continue

			if opts.StrandBiasReads:
				continue

			if opts.AlignmentScore:
				continue

			if opts.DuplicateReference:
				continue

			if opts.DuplicateAlternate:
				continue

			if opts.DeltaDuplicate:
				continue

			if opts.MeanRefQscore:
				continue

			if opts.MeanAltQscore:
				continue

			if opts.PercentageReadsQscoreFilter:
				continue

			if opts.PercentageNumberOfMismatch:
				continue

			out.write('\t'.join(variant) + '\n')

main()
