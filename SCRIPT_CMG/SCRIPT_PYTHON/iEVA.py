import argparse
import re
import string
import sys
from Bio.Seq import Seq
from Bio import SeqIO, Entrez
from Bio import SeqUtils
from pyfasta import Fasta
from repDNA.nac import Kmer
import pysam
import scipy.stats as stats


#-----------------------------------------------------FUNCTIONS-----------------------------------------------



#--------------------------------------------------CHECK_PARAMETER--------------------------------------------



def Check_Parameters(option,argv,min,max):

	if argv > max or argv < min:
		sys.exit("Argument " + str(option) + ':' + ' expecated value in range ' + str(min) + '-' + str(max) + '\n')



#-------------------------------------------------CHECK_ZERO_NUMBER-------------------------------------------



def Check_Zero(number):

	if number == 0.0:
		number = int(number)
	else:
		pass

	return number


#-----------------------------------------------EXTRACT_VARIANT_TYPE------------------------------------------



def Extract_Variant_Type(variant, H_CHR):

	CHR = variant[H_CHR.index('#CHROM')]
	REF = variant[H_CHR.index('REF')]
	ALT = variant[H_CHR.index('ALT')]

	#Per prima cosa vado a definire la classe della variante, se SNV, IN/DEL, INDEL o Sequence Alteration
	if re.match("[A|T|C|G]*$",REF) and re.match("[A|T|C|G]*$",ALT):

		if len(REF) > len(ALT) and len(ALT) == 1:
			Variant_Class = 'Del'

		elif len(ALT) > len(REF)  and len(REF) == 1:
			Variant_Class = 'Ins'

		elif len(ALT) == len(REF) or len(REF) == len(ALT) and len(REF) == 1 and len(ALT) == 1:
			Variant_Class = 'snv'

		elif len(ALT) == len(REF) or len(REF) == len(ALT) and len(REF) != 1 and len(ALT) != 1:
			Variant_Class = 'Alt'

		else:
			Variant_Class = '.'

	else:

		Variant_Class = '.'

	
	return Variant_Class



#--------------------------------------------------EXTRACT_READS_INFO-----------------------------------------



def Extract_Reds_Info(read, Reads_Info):

	Reads_Info['Total_Reads_Unfilter'] += 1

	if read.alignment.mapping_quality == 0 and read.alignment.is_duplicate == False:
		Reads_Info['Mapping_Quality_Zero'] += 1

	if read.alignment.is_unmapped and read.alignment.is_duplicate == False:
		Reads_Info['Unmapped_reads'] += 1	

	if read.alignment.is_duplicate:
		Reads_Info['Duplicate_reads'] += 1

	if read.alignment.is_duplicate == False:
		Reads_Info['Total_Reads_No_Dup'] += 1		

	if read.alignment.is_paired == False and read.alignment.is_duplicate == False:
		Reads_Info['Not_Paired_Reads'] += 1

	if read.alignment.is_proper_pair == False and read.alignment.is_duplicate == False:
		Reads_Info['Not_Proper_Paired_Reads'] += 1

	if read.alignment.is_supplementary and read.alignment.is_duplicate == False:
		Reads_Info['Supplementary_Align'] += 1

	if read.alignment.is_secondary and read.alignment.is_duplicate == False:
		Reads_Info['Not_Primary_Align'] += 1

	if read.alignment.is_paired and read.alignment.is_duplicate == False and read.alignment.is_unmapped == False \
	and read.alignment.is_secondary == False and read.alignment.is_supplementary == False:

		Reads_Info['Total_Reads_Filtered'] += 1
		Reads_Info['Alignment_Score'] += read.alignment.get_tag('AS')

		if read.alignment.is_read1:
		 	Reads_Info['is_read1'] += 1
		 	if read.alignment.mate_is_reverse:
		 		Reads_Info['is_read1_forward'] += 1
			elif read.alignment.is_reverse:
				Reads_Info['is_read1_reverse'] += 1
		 	
		elif read.alignment.is_read2:
			Reads_Info['is_read2'] += 1
			if read.alignment.is_reverse:
		 		Reads_Info['is_read2_reverse'] += 1
		 	if read.alignment.mate_is_reverse:
		 		Reads_Info['is_read2_forward'] += 1

	return Reads_Info



#-----------------------------------------------------CHECK_BAM-----------------------------------------------



def Check_Bam(Sample_list, bam_list):

	#Ora vado a controllare se il Sample name nel bam nella tag 'SM' e' lo stesso che trovo nella colonna sample del vcf.
	#Per farlo, il comando bamfile.header prende l'header del bam come una dictionary i cui valori sono composti da liste con
	#al loro interno annidate altre dictionary. Per estrarre il nome uso allora il comando bam_header.get('RG')[0].get('SM')

	Sample_dict = {}

	path_list = open(bam_list, 'r')

	#Apro i bam presenti in bam_list per processare le reads. Per ogni bam estraggo la feature
	for path in path_list:
		
		path = path.rstrip()
		
		bamfile = pysam.AlignmentFile(path, "rb")

		bam_header = bamfile.header

		Header_Sample_Name = bam_header.get('RG')[0].get('SM')

		#Se il campo 'SM' nell'header del bam corrisponde a quello del vcf, allora assegno ad ogni campione il suo path e rimuovo dalla
		#Sample_list il nome trovato.
		if Header_Sample_Name in Sample_list:

			Sample_dict[Header_Sample_Name] = path

			Sample_list.remove(Header_Sample_Name)

		else:

			continue

	path_list.close()

	#Se la sample list ha elementi al suo interno alllora non tutti i path del bam sono stati inseriti ed il tool da un Warn
	if len(Sample_list) != 0:
		for elem in Sample_list:
			print '\n' + "Sample " + str(elem) + " will not be annotated. Missing bam file or malformed sample field in vcf file." + '\n'
			message = ''
	return Sample_dict



#-----------------------------------------------SIMPLE_REPEATS_FINDER-----------------------------------------



def SimpleRepeats_Finder(Reference, variant, H_CHR, Variant_Class, WindowSize):

	WindowSize = WindowSize / 2

	#Leggo il file fasta reference
	Ref = Fasta(Reference)

	#Estraggo Chrom, Pos, Ref ed Alt
	sorted(Ref.keys())
	CHR = variant[H_CHR.index('#CHROM')]
	POS = int(variant[H_CHR.index('POS')])-1
	REF = variant[H_CHR.index('REF')]
	ALT = variant[H_CHR.index('ALT')]
	
	#Estraggo la sequenza dalla reference corrispondente alla variante trovata
	sequence = Ref[CHR][POS-WindowSize:POS+WindowSize+1]

	#Trasformo la sequenza in lettere maiuscole
	seq = sequence.upper()

	#CALCOLO COMPOSIZIONE IN PSEUDONUCLEOTIDI:
	kmer = Kmer(k=2,normalize=True)
	Dinuc = kmer.make_kmer_vec([seq])
	Dinucloetides = Dinuc[0]
	Combo_Seq = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
	#Psuedo_Nucleotide = dict(zip(Combo_Seq,Dinucloetides))

	#Definisco i pattern da andare a trovare nella sequenza
	Homopolymer_pattern = re.compile(r'(A){5,}|(G){5,}|(C){5,}|(T){5,}')

	#Estraggo gli omopolimeri trovati in una lista contenente tuple. Es: [('AAAAAA',4,6),('TTTTTTTTTT',6,10)] in cui il primo elemento indica 
	#il pattern di omopolimeri, il secondo l'indice in cui inizia ed il terzo la lunghezza del pattern
	Homopolymer_list = [(m.group(), m.start(), len(m.group())) for m in re.finditer(Homopolymer_pattern, seq)]

	RepeatSeq = '0'
	RepeatSeq_Lenght = '.'

	# La variante (se SNV) deve cadere nell'omopolimero mentre puo' iniziare in posizione POS del vcf se INDEL dal momento che gli INDEl sono
	# left-align.
	for poly in Homopolymer_list:

		End_Hom = poly[1]+poly[2]-1

		if poly[1] == WindowSize+1 and Variant_Class is 'Deletion' or poly[1] == WindowSize+1 and Variant_Class is 'Insertion':
			RepeatSeq = '2'
			RepeatSeq_Lenght = str(poly[2])
			break
		
		elif poly[1] <= WindowSize and End_Hom >= WindowSize and Variant_Class is 'SNV':
			RepeatSeq = '2'
			RepeatSeq_Lenght = str(poly[2])
			break

		else:
			continue

	#Passo ora all'identificazione delle Simple Repeats come come poly(AT|ATT|ATAA|etc...). Il condizionale e' lo stesso di prima
	if '0' in RepeatSeq:

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

			End_Rep = pattern[1]+pattern[2]-1

			if pattern[1] == WindowSize+1 and Variant_Class is 'Deletion' or pattern[1] == WindowSize+1 and Variant_Class is 'Insertion':
				RepeatSeq = '1'
				RepeatSeq_Lenght = (pattern[2])
				break
			
			elif pattern[1] <= WindowSize and End_Rep >= WindowSize and Variant_Class is 'SNV':
				RepeatSeq = '1'
				RepeatSeq_Lenght = str(pattern[2])
				break

			else:
				continue

	#print str(sequence) + ' ' + str(CHR) + ' ' + str(POS+1) + ' ' + str(REF) + ' ' + str(ALT) + ' ' + str(RepeatSeq) + ' ' + str(RepeatSeq_Lenght)
	
	return RepeatSeq, RepeatSeq_Lenght, Dinucloetides



#---------------------------------------------------REPEAT_MASKER---------------------------------------------



def Check_Repeat_Masker(Reference, variant, H_CHR):

	#Leggo il file fasta reference
	Ref = Fasta(Reference)

	#Estraggo Chrom, Pos, Ref ed Alt
	sorted(Ref.keys())
	CHR = variant[H_CHR.index('#CHROM')]
	POS = int(variant[H_CHR.index('POS')])-1

	#Estraggo la sequenza dalla reference corrispondente alla variante trovata
	sequence = Ref[CHR][POS]

	#Controllo se la variante cade su sequenza in lettera minuscola (trovata dal repeat masker) oppure se maiuscola (sequenza non ripetuta)
	if sequence.isupper():
		RM = '0'
	else:
		RM = '1'

	return RM



#----------------------------------------------------GC_CONTENT-----------------------------------------------



def GC_content(Reference, variant, H_CHR, WindowSize):

	WindowSize = WindowSize / 2

	#Leggo il file fasta reference
	Ref = Fasta(Reference)

	#Estraggo Chrom, Pos, Ref ed Alt
	sorted(Ref.keys())
	CHR = variant[H_CHR.index('#CHROM')]
	POS = int(variant[H_CHR.index('POS')])-1
	
	#Estraggo la sequenza dalla reference corrispondente alla variante trovata
	sequence = Ref[CHR][POS-WindowSize:POS+WindowSize+1]

	#Trasformo la sequenza in lettere maiuscole
	seq = sequence.upper()

	GC = str(round(SeqUtils.GC(seq), 3))

	return GC



#--------------------------------------------------ALIGNMENT-SCORE-------------------------------------------



def Sequence_Annotator(variant, Sample_dict, H_CHR, Reference, Variant_Class, SNVMinBaseQuality, SNVMinMappingQuality, IndelMinBaseQuality, IndelMinMappingQuality):

	CHR = variant[H_CHR.index('#CHROM')]
	POS = int(variant[H_CHR.index('POS')])-1
	REF = variant[H_CHR.index('REF')]
	ALT = variant[H_CHR.index('ALT')]

	variant_lenght = max(len(REF), len(ALT))

	STOP = int(variant[H_CHR.index('POS')]) + variant_lenght - 1

	Sample_Stat = {}

	#Estraggo per tutti i campioni l'AS score
	for sample in Sample_dict.keys():

		#Skippo le varianti che hanno * in alternativo, le trovo solitamente solo in GATK
		if '*' in ALT or ',' in ALT:

			Variant_Stat = {}

			Variant_Stat['Mean_Alignment_Score'] = '.'
			Variant_Stat['Percentage_Unmapped_Reads'] = '.'
			Variant_Stat['Strand_Bias_Reads'] = '.'
			Variant_Stat['VariantClass'] = '.'
			Variant_Stat['Read_Ref'] = '.'
			Variant_Stat['Read_Alt'] = '.'
			Variant_Stat['Qual_Alt'] = '.'
			Variant_Stat['Qual_Ref'] = '.'
			Variant_Stat['Percentage_Dup_Alt'] = '.'
			Variant_Stat['Percentage_Dup_Ref'] = '.'
			Variant_Stat['Number_Read_Dup_Alt'] = '.'
			Variant_Stat['Number_Read_Dup_Ref'] = '.'
			Variant_Stat['Delta_Dupicate'] = '.'
			Variant_Stat['Total_Duplicate_Reads'] = '.'
			Variant_Stat['Total_Reads_Unfilter'] = '.'
			Variant_Stat['Clipped_Reads_Ref'] = '.'
			Variant_Stat['Clipped_Reads_Alt'] = '.'
			Variant_Stat['Percentage_ClipRef'] = '.'
			Variant_Stat['Percentage_ClipAlt'] = '.'
			Variant_Stat['Coverage'] = '.'
			Variant_Stat['Percentage_Supplementary_Align'] = '.'
			Variant_Stat['Percentage_Not_Primary_Alignment_Reads'] = '.'
			Variant_Stat['Percentage_Not_Paired_Reads'] = '.'
			Variant_Stat['Percentage_Not_Proper_Paired_Reads'] = '.'
			Variant_Stat['Mapping_Quality_Zero'] = '.'

			Sample_Stat[sample] = Variant_Stat

		else:

			Variant_Stat = {}

			Reads_Info = {}

			Reads_Info['Total_Reads_No_Dup'] = 0
			Reads_Info['Total_Reads_Unfilter'] = 0
			Reads_Info['Total_Reads_Filtered'] = 0
			Reads_Info['Coverage'] = 0
			Reads_Info['Duplicate_reads'] = 0
			Reads_Info['Not_Paired_Reads'] = 0
			Reads_Info['Not_Proper_Paired_Reads'] = 0
			Reads_Info['Dup_Read_Alt'] = 0
			Reads_Info['Dup_Read_Ref'] = 0
			Reads_Info['Alignment_Score'] = 0
			Reads_Info['Unmapped_reads'] = 0
			Reads_Info['Mapping_Quality_Zero'] = 0
			Reads_Info['Read_Ref'] = 0
			Reads_Info['Read_Alt'] = 0
			Reads_Info['Base_Alt_Qual'] = 0
			Reads_Info['Base_Ref_Qual'] = 0
			Reads_Info['Clipped_Reads_Ref'] = 0
			Reads_Info['Clipped_Reads_Alt'] = 0
			Reads_Info['Supplementary_Align'] = 0
			Reads_Info['Not_Primary_Align'] = 0

			Reads_Info['is_read1'] = 0
			Reads_Info['is_read1_forward'] = 0
			Reads_Info['is_read1_reverse'] = 0

			Reads_Info['is_read2'] = 0
			Reads_Info['is_read2_forward'] = 0
			Reads_Info['is_read2_reverse'] = 0


			bampile = pysam.AlignmentFile(Sample_dict.get(sample), "rb")


			if Variant_Class == 'SNV':

				for pileupcolumn in bampile.pileup(CHR, POS, POS+variant_lenght, stepper='nofilter'):

					if pileupcolumn.pos == POS:

						for pileupread in pileupcolumn.pileups:

							#Estraggo le feature dal bam con la funzione Extract_Reds_Info
							Reads_Info = Extract_Reds_Info(pileupread, Reads_Info)

							if not pileupread.is_del and not pileupread.is_refskip:

								if pileupread.alignment.query_sequence[pileupread.query_position] == ALT \
								and pileupread.alignment.is_paired and not pileupread.alignment.is_unmapped \
								and not pileupread.alignment.is_secondary and not pileupread.alignment.is_supplementary:

									if pileupread.alignment.is_duplicate:
										Reads_Info['Dup_Read_Alt'] += 1

									elif pileupread.alignment.query_qualities[pileupread.query_position] < SNVMinBaseQuality \
									or pileupread.alignment.mapping_quality < SNVMinMappingQuality:
										continue

									else:
										Reads_Info['Coverage'] += 1
										Reads_Info['Read_Alt'] += 1
										Reads_Info['Base_Alt_Qual'] += pileupread.alignment.query_qualities[pileupread.query_position]

										if 'S' in pileupread.alignment.cigarstring or 'H' in pileupread.alignment.cigarstring:

											Reads_Info['Clipped_Reads_Alt'] += 1


								elif pileupread.alignment.query_sequence[pileupread.query_position] == REF \
								and pileupread.alignment.is_paired and not pileupread.alignment.is_unmapped \
								and not pileupread.alignment.is_secondary and not pileupread.alignment.is_supplementary:

									if pileupread.alignment.is_duplicate:
										Reads_Info['Dup_Read_Ref'] += 1

									elif pileupread.alignment.query_qualities[pileupread.query_position] < SNVMinBaseQuality \
									or pileupread.alignment.mapping_quality < SNVMinMappingQuality:
										continue

									else:
										Reads_Info['Coverage'] += 1
										Reads_Info['Read_Ref'] += 1
										Reads_Info['Base_Ref_Qual'] += pileupread.alignment.query_qualities[pileupread.query_position]

										if 'S' in pileupread.alignment.cigarstring or 'H' in pileupread.alignment.cigarstring:

											Reads_Info['Clipped_Reads_Ref'] += 1




			if Variant_Class == 'Deletion':

				for pileupcolumn in bampile.pileup(reference=CHR, start=POS, end=POS+variant_lenght-1, stepper='nofilter'):

					if pileupcolumn.pos == POS:

						for pileupread in pileupcolumn.pileups:

							Reads_Info = Extract_Reds_Info(pileupread, Reads_Info)

							Read_Tuple = pileupread.alignment.get_aligned_pairs(with_seq=True)

							# Includo nella variante la base in cui parte la delezione (la left-align come nel vcf) e quella
							# dopo la delezione per verificare in seguito che sia un match o comunque che non sia un'altra delezione 
							# Posso farlo tramite il seguente ciclo for che restituisce l'indice in cui trovo la POS data. Se non lo trova,
							# Exact_Match = [] --> Sara' vuoto

							Start_Index = 0
							Stop_Index = 0

							for index, val in enumerate(Read_Tuple):

								if val[1] == POS and val[0] is not None and val[2].isupper() and val[2] == REF[0]:
									Start_Index = index

								if val[1] == STOP+1 and val[0] is not None and val[2].isupper():
									Stop_Index = index
									break

							Exact_Match = Read_Tuple[Start_Index:Stop_Index]


							# -----> N.B Se la read termina prima di stop index, oppure la base subito dopo l'indel descritto in POS
							#        e' un'altra delezione, allora Indel_sequence sara' una stringa che indica TUTTA la seq della read.
							#        Questo accade perche' Exact_Match[-1][0] e' uguale a None
							try:
								Indel_Sequence = pileupread.alignment.query_sequence[Exact_Match[0][0]:(Exact_Match[-1][0])]
								#print Indel_Sequence
							except:
								Indel_Sequence = ''

							# Exact_Match sara' vuoto se la read e' piu' corta della variante (finisce prima delle delezione) e quindi posso
							# tagliare fuori le read piu' corte. Allo stesso modo, se trovo un'INSERZIONE allora la len(Extract_Match) sara' 
							# piu' grande e quindi la escludo:
							if Exact_Match == [] or len(Exact_Match) > variant_lenght+1:
								continue

							else:

								# 1) Con il comando pileupread.indel verifico la lunghezza dell'eventuale delezione che segue la pos == -(variant_lenght-1)
								# 2) Filtro le read seconodo le specifiche (es: pileupread.alignment.is_unmapped)
								# 3) Verifico che la lunghezza della seq estratta sia uguale a 2 (Base di POS e base di POS+1): len(Indel_Sequence) == ALT+1
								# 4) Gia' prima ho verificato che POS e STOP siano entrambi match esatti su REF
								# 5) Se sulla read trovo un Indel che corrisponde a POS, allora la seq della read avra' in quelle posizioni solo una base
								#    (corrispondente appunto a POS) e quindi len(Indel_Sequence) == 1

								if pileupread.indel == -(variant_lenght-1) and not pileupread.alignment.is_unmapped \
								and pileupread.alignment.is_proper_pair and pileupread.alignment.is_paired \
								and not pileupread.is_refskip and not pileupread.alignment.is_secondary \
								and not pileupread.alignment.is_supplementary and not pileupread.is_tail \
								and len(Indel_Sequence) == 1:

									if pileupread.alignment.is_duplicate:
										Reads_Info['Dup_Read_Alt'] += 1

									elif pileupread.alignment.query_qualities[pileupread.query_position] < IndelMinBaseQuality \
									or pileupread.alignment.mapping_quality < IndelMinMappingQuality:
										continue

									else:

										Reads_Info['Coverage'] += 1
										Reads_Info['Read_Alt'] += 1
										Alt_Qual = round(float(sum([pileupread.alignment.query_qualities[pileupread.query_position]] + [pileupread.alignment.query_qualities[pileupread.query_position+1]]))/float(2),2)
										Reads_Info['Base_Alt_Qual'] += Alt_Qual

										if 'S' in pileupread.alignment.cigarstring or 'H' in pileupread.alignment.cigarstring:
											Reads_Info['Clipped_Reads_Alt'] += 1
							

								# Per la REF della delezione e' un po' piu' complesso. Prima vado a verificare che la variante inizi in una pos
								# per cui ci sia un match (quindi non stia partendo gia' da una delezione) e che rispetti i filtraggi. Poi
								# passo a valutare se la prima base della variante espressa nella coordinata del vcf e' un match con la ref anche
								# in Exact_Match trovato prima (Exact_Match[-1][0] != None) e che la base seguente all'ultima espressa
								# nella delezione sia anch'essa un match (or Exact_Match[-1][1] != None)


								elif pileupread.alignment.is_paired and pileupread.alignment.is_proper_pair and not pileupread.alignment.is_unmapped \
								and not pileupread.is_refskip and not pileupread.query_position == None \
								and not pileupread.alignment.is_secondary and not pileupread.alignment.is_supplementary:

									Not_Ref = ''

									# Verifico che sia delezione e non inserzione ( if values[0] != None and values[1] ) e che la stessa REF non
									# contenga SNP al suo interno ( values[2].isupper() )
									for values in Exact_Match[1:-1]:
										if values[0] is not None and values[1] is not None and values[2].isupper():
											Not_Ref = 'False'
										
										else:
											Not_Ref = 'True'
											break
								
									# Filtro le read per cui la prima base nel vcf abbia un quality score minore di 8 oppure un mapping quality
									# minore di 15
									if pileupread.alignment.query_qualities[pileupread.query_position] < IndelMinBaseQuality \
									or pileupread.alignment.mapping_quality < IndelMinMappingQuality:
										Not_Ref = 'True'

									if Not_Ref is 'False':

										if pileupread.alignment.is_duplicate:
											Reads_Info['Dup_Read_Ref'] += 1

										else:
											Reads_Info['Coverage'] += 1
											Reads_Info['Read_Ref'] += 1
											Base_Qual = round(float(sum(pileupread.alignment.query_qualities[pileupread.query_position:pileupread.query_position+variant_lenght]))/float(variant_lenght),2)
											Reads_Info['Base_Ref_Qual'] += Base_Qual
											
											if 'S' in pileupread.alignment.cigarstring or 'H' in pileupread.alignment.cigarstring:
												Reads_Info['Clipped_Reads_Ref'] += 1




			if Variant_Class == 'Insertion':

				for pileupcolumn in bampile.pileup(reference=CHR, start=POS, end=POS+variant_lenght-1, stepper='nofilter'):

					if pileupcolumn.pos == POS:

						for pileupread in pileupcolumn.pileups:

							Reads_Info = Extract_Reds_Info(pileupread, Reads_Info)

							Read_Tuple = pileupread.alignment.get_aligned_pairs(with_seq=True)

							# Stesso procedimento utilizzato per gli Indel

							Start_Index = 0
							Stop_Index = 0

							for index, val in enumerate(Read_Tuple):

								#Estraggo Extract_Match
								if val[1] == POS and val[0] is not None and val[2].isupper():
									Start_Index = index

								if val[1] == (POS+2) and val[0] is not None and val[2].isupper():
									Stop_Index = index
									break

							Exact_Match = Read_Tuple[Start_Index:Stop_Index]


							# Con questo try verifico che la read non finisca subito dopo POS, altrimenti Indel_sequence da errore
							try:
								Indel_Sequence = pileupread.alignment.query_sequence[Exact_Match[0][0]:(Exact_Match[-1][0])]
							except:
								Indel_Sequence = ''


							# Se Exact_Match e' vuoto, allora salta
							if Exact_Match == []:
								continue
							
							# Per prima cosa si richiede che la grandezza della sequenza sia pari a quella riportata in vcf.
							# In questo modo escludo varianti in cui ad esempio l'ALT e' TAAA piuttosto che TAA
							elif len(Exact_Match) == variant_lenght+1:

								# 1) Con il comando pileupread.indel verifico la lunghezza dell'eventuale inserzione che segue la pos == (variant_lenght-1)
								# 2) Filtro le read seconodo le specifiche (es: pileupread.alignment.is_unmapped)
								# 3) Verifico che la sequenza della read estratta sia uguale all'inserzione: Indel_Sequence == ALT
								# 4) Gia' in precedenza quando generato Exact_Match, ho verificato che start e stop siano MATCH

								if pileupread.indel == (variant_lenght-1) and not pileupread.alignment.is_unmapped \
								and pileupread.alignment.is_proper_pair \
								and pileupread.alignment.is_paired and not pileupread.is_refskip \
								and not pileupread.alignment.is_secondary and not pileupread.alignment.is_supplementary \
								and not pileupread.is_tail and Indel_Sequence == ALT \
								and Exact_Match[-1][0] is not None and Exact_Match[-1][2].isupper() \
								and Exact_Match[0][0] is not None and Exact_Match[0][2].isupper():

									if pileupread.alignment.is_duplicate:
										Reads_Info['Dup_Read_Alt'] += 1

									elif float(sum(pileupread.alignment.query_qualities[Exact_Match[1][0]:Exact_Match[-1][0]]))/float(variant_lenght-1) < IndelMinBaseQuality \
									or pileupread.alignment.mapping_quality < IndelMinMappingQuality:
										continue

									else:
										Reads_Info['Coverage'] += 1
										Reads_Info['Read_Alt'] += 1
										Alt_Qual = round(float(sum(pileupread.alignment.query_qualities[Exact_Match[1][0]:Exact_Match[-1][0]]))/float(variant_lenght-1),2)
										Reads_Info['Base_Alt_Qual'] += Alt_Qual

										if 'S' in pileupread.alignment.cigarstring or 'H' in pileupread.alignment.cigarstring:
											Reads_Info['Clipped_Reads_Alt'] += 1
							

							# Passo ad analizzare le REF controllando per prima cosa che la lunghezza di Exact_Match sia pari a 2
							# ovvero POS e POS+1 (non contiene quindi Indel), come sopra
							elif len(Exact_Match) == 2:

								if pileupread.alignment.is_paired and pileupread.alignment.is_proper_pair \
								and not pileupread.alignment.is_unmapped \
								and not pileupread.is_refskip and not pileupread.alignment.is_secondary \
								and not pileupread.alignment.is_supplementary \
								and Exact_Match[0][0] is not None and Exact_Match[0][2].isupper() \
								and Exact_Match[-1][0] is not None and Exact_Match[-1][2].isupper():

									if pileupread.alignment.is_duplicate:
											Reads_Info['Dup_Read_Ref'] += 1

									elif float(sum(pileupread.alignment.query_qualities[Exact_Match[0][0]:(Exact_Match[-1][0])]))/float(2) < IndelMinBaseQuality \
									or pileupread.alignment.mapping_quality < IndelMinMappingQuality:
										continue

									else:
										Reads_Info['Coverage'] += 1
										Reads_Info['Read_Ref'] += 1
										Base_Qual = round(float(sum(pileupread.alignment.query_qualities[Exact_Match[0][0]:(Exact_Match[-1][0])+1]))/float(2),2)
										Reads_Info['Base_Ref_Qual'] += Base_Qual
										
										if 'S' in pileupread.alignment.cigarstring or 'H' in pileupread.alignment.cigarstring:
											Reads_Info['Clipped_Reads_Ref'] += 1

							else:
								continue

		#Tutte le features estratte le vado ad inserire in un dict con chiave il la feature e valore il corrispondente risultato.
		#Il dict verra' poi utilizzato come valore di un'altra dict le cui chiavi sono i campioni. Ad esempio avremo:
		#{'20161125_02_Cardio': {'SBR': 1.0, 'AS': 0}, '20161125_01_Cardio': {'SBR': 1.0, 'AS': '139.333333333'}}

		try:
			Variant_Stat['Mean_Alignment_Score'] = round(float(Reads_Info['Alignment_Score'])/float(Reads_Info['Total_Reads_No_Dup']),3)
		except:
			Variant_Stat['Mean_Alignment_Score'] = '.'

		try:
			Variant_Stat['Percentage_Unmapped_Reads'] = round(float(Reads_Info['Unmapped_reads'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Percentage_Unmapped_Reads'] = '.'

		#Calcolo lo strand bias con il Fisher's exact test
		try:
			oddsratio, SBR = stats.fisher_exact([Reads_Info['is_read1_forward'], Reads_Info['is_read1_reverse']], [Reads_Info['is_read2_forward'], Reads_Info['is_read2_reverse']])
		except:
			SBR = '.'

		try:
			Variant_Stat['Strand_Bias_Reads'] = round(SBR,4)
		except:
			Variant_Stat['Strand_Bias_Reads'] = '.'

		try:
			Variant_Stat['VariantClass'] = Reads_Info['Variant_Class']
		except:
			Variant_Stat['VariantClass'] = '.'

		try:
			Variant_Stat['Read_Ref'] = Reads_Info['Read_Ref']
		except:
			Variant_Stat['Read_Ref'] = '.'

		try:
			Variant_Stat['Read_Alt'] = Reads_Info['Read_Alt']
		except:
			Variant_Stat['Read_Alt'] = '.'

		try:
			Variant_Stat['Qual_Alt'] = round(float(Reads_Info['Base_Alt_Qual'])/float(Reads_Info['Read_Alt']),2)
		except:
			Variant_Stat['Qual_Alt'] = '.'

		try:
			Variant_Stat['Qual_Ref'] = round(float(Reads_Info['Base_Ref_Qual'])/float(Reads_Info['Read_Ref']),2)
		except:
			Variant_Stat['Qual_Ref'] = '.'

		try:
			Variant_Stat['Percentage_Dup_Alt'] = round(float(Reads_Info['Dup_Read_Alt'])/float(Reads_Info['Read_Alt']+Reads_Info['Dup_Read_Alt']),4)
		except:
			Variant_Stat['Percentage_Dup_Alt'] = '.'

		try:
			Variant_Stat['Percentage_Dup_Ref'] = round(float(Reads_Info['Dup_Read_Ref'])/float(Reads_Info['Read_Ref']+Reads_Info['Dup_Read_Ref']),4)
		except:
			Variant_Stat['Percentage_Dup_Ref'] = '.'

		try:
			Variant_Stat['Number_Read_Dup_Alt'] = Reads_Info['Dup_Read_Alt']
		except:
			Variant_Stat['Number_Read_Dup_Alt'] = '.'

		try:
			Variant_Stat['Number_Read_Dup_Ref'] = Reads_Info['Dup_Read_Ref']
		except:
			Variant_Stat['Number_Read_Dup_Ref'] = '.'

		try:
			Variant_Stat['Delta_Dupicate'] = round(float(Variant_Stat['PDR']-Variant_Stat['PDA']),4)
		except:
			Variant_Stat['Delta_Dupicate'] = '.'

		try:
			Variant_Stat['Total_Duplicate_Reads'] = round(float(Reads_Info['Duplicate_reads'])/float(Reads_Info['Total_Reads_Unfilter']),4)
		except:
			Variant_Stat['Total_Duplicate_Reads'] = '.'

		try:
			Variant_Stat['Total_Reads_Unfilter'] = Reads_Info['Total_Reads_Unfilter']
		except:
			Variant_Stat['Total_Reads_Unfilter'] = '.'

		try:
			Variant_Stat['Clipped_Reads_Ref'] = Reads_Info['Clipped_Reads_Ref']
		except:
			Variant_Stat['Clipped_Reads_Ref'] = '.'

		try:
			Variant_Stat['Clipped_Reads_Alt'] = Reads_Info['Clipped_Reads_Alt']
		except:
			Variant_Stat['Clipped_Reads_Alt'] = '.'

		try:
			Variant_Stat['Percentage_ClipRef'] = round(float(Reads_Info['Clipped_Reads_Ref'])/float(Reads_Info['Read_Ref']),4)
		except:
			Variant_Stat['Percentage_ClipRef'] = '.'

		try:
			Variant_Stat['Percentage_ClipAlt'] = round(float(Reads_Info['Clipped_Reads_Alt'])/float(Reads_Info['Read_Alt']),4)
		except:
			Variant_Stat['Percentage_ClipAlt'] = '.'

		try:
			Variant_Stat['Coverage'] = Reads_Info['Coverage']
		except:
			Variant_Stat['Coverage'] = '.'

		try:
			Variant_Stat['Percentage_Supplementary_Align'] = round(float(Reads_Info['Supplementary_Align'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Percentage_Supplementary_Align'] = '.'

		try:
			Variant_Stat['Percentage_Not_Primary_Alignment_Reads'] = round(float(Reads_Info['Not_Primary_Align'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Percentage_Not_Primary_Alignment_Reads'] = '.'

		try:
			Variant_Stat['Percentage_Not_Paired_Reads'] = round(float(Reads_Info['Not_Paired_Reads'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Percentage_Not_Paired_Reads'] = '.'

		try:
			Variant_Stat['Percentage_Not_Proper_Paired_Reads'] = round(float(Reads_Info['Not_Proper_Paired_Reads'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Percentage_Not_Proper_Paired_Reads'] = '.'

		try:
			Variant_Stat['Mapping_Quality_Zero'] = round(float(Reads_Info['Mapping_Quality_Zero'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Mapping_Quality_Zero'] = '.'


		#Inserisco la dict come valore della dict riferita ai campioni:
		Sample_Stat[sample] = Variant_Stat
		#print Reads_Info
		Reads_Info = {}


	return Sample_Stat



#-------------------------------------------------------MAIN--------------------------------------------------



def main():

	parser = argparse.ArgumentParser(description='iEVA tool: identification and Extraction of Variant Attributes.')
	
	Required = parser.add_argument_group('Required arguments')
	Required.add_argument('-I','--input',required=True,metavar='path/to/input',help="Input file in vcf format.")
	Required.add_argument('-Ref','--reference',required=True,metavar='path/to/fasta',help="Reference file in fasta format.")
	Required.add_argument('-O','--outfile',required=True,metavar='path/to/outfile',help="Output vcf file.")
	Required.add_argument('-L','--list',metavar='path/to/bamlist',help="REQUIRED only to extract features from bam file. Bam list file, including path, one for each raw. Sample name in SAMPLE column in vcf has to be the same of 'SM' header in bam file. Index file bai in same path of bam file.")

	
	Sequence = parser.add_argument_group('Sequence optional arguments')
	Sequence.add_argument('-WS','--WindowSize',default=40,type=int,metavar='20-600',help="Set window size for -SR -SRL -GC and -PNC command. Default=40, Min=20, Max=600")
	Sequence.add_argument('-SR','--SimpleRepeat',action="store_true",help="Enable Simple Repeat Finder. Check if a variant fall in a Repeated Sequence. 0 for None, 1 for Simple Repeat Sequence, 2 for Homopolymer sequence.")
	Sequence.add_argument('-SRL','--SimpleRepeatLength',action="store_true",help="Report length of repeated sequence")
	Sequence.add_argument('-PNC','--PseudoNucleotidesComposition',action="store_true",help="Describe nucleotide sequence using Pseudo Nucleotide Composition with Kmer size of 2. Values reported as: AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT")
	Sequence.add_argument('-RM','--RepeatMasker',action="store_true",help="Check if variant falls into a repeated sequence found with Repeat Masker tool.")
	Sequence.add_argument('-GC','--gcContent',action="store_true",help="Sequence GC content.")
	Sequence.add_argument('-VC','--VariantClass',action="store_true",help="Annotate variant class: SNV (snv), Insertion (Ins), Deletion (Del), SequenceAlteration (Alt)")	


	Bam = parser.add_argument_group('Bam optional arguments')
	Bam.add_argument('-SBR','--StrandBiasReads',action="store_true",help="Strand Baias based on read orientation using Fisher's exact test.")
	Bam.add_argument('-UnMap','--UnMappedReads',action="store_true",help="Percentage of unmapped reads for variant position.")
	Bam.add_argument('-MQ0','--MappingQualityZero',action="store_true",help="Percentage of reads having Mapping Quaility 0.")
	Bam.add_argument('-NPA','--NotPrimaryAlignment',action="store_true",help="Percentage of reads mapping position where alignment is not primary, i.e. they map also other genomic location. Flagged as NotPrimaryAlignment. Tag in FORMAT field.")
	Bam.add_argument('-SA','--SupplementaryAlignment',action="store_true",help="Percentage of reads flagged as supplementary alignment. These reads suggest structural a chimeric alignments (i.e. one read fragment map to another location).")
	Bam.add_argument('-NP','--NotPairedReads',action="store_true",help="Percentage of non paired reads.")
	Bam.add_argument('-NPP','--NotProperPairedReads',action="store_true",help="Percentage of non proper paired reads, like ambiguous pairing.")
	Bam.add_argument('-AS','--AlignmentScore',action="store_true",help="Mean Alignment Score of reads mapping position.")
	Bam.add_argument('-NDT','--NumberTotalDupReads',action="store_true",help="Total number of duplicate reads mapping position.")
	Bam.add_argument('-NDR','--NumberReadDupRef',action="store_true",help="Number of duplicate reads mapping REF allele.")
	Bam.add_argument('-NDA','--NumberReadDupAlt',action="store_true",help="Number of duplicate reads mapping ALT allele.")
	Bam.add_argument('-DR','--DuplicateReference',action="store_true",help="Percentage of REF reads marked as duplicates.")
	Bam.add_argument('-DA','--DuplicateAlternate',action="store_true",help="Percentage of ALT reads marked as duplicates.")
	Bam.add_argument('-DDup','--DeltaDuplicate',action="store_true",help="Difference for percentage duplicates in REF and ALT (REF-ALT). Negative values show preference in duplicate for REF allele, positive for ALT.")
	Bam.add_argument('-SNVmbq','--SNVMinBaseQuality',default=12,type=int,metavar='0-66',help="Minimum Base Quality treshold for base supporting SNV position. Default=12")
	Bam.add_argument('-SNVmmq','--SNVMinMappingQuality',default=30,type=int,metavar='0-60',help="Minimum Mapping Quality treshold for reads supporting SNV position. Default=30")
	Bam.add_argument('-INDELmbq','--IndelMinBaseQuality',default=10,type=int,metavar='0-66',help="Minimum Base Quality treshold for base supporting InDel position. Default=10")
	Bam.add_argument('-INDELmmq','--IndelMinMappingQuality',default=20,type=int,metavar='0-60',help="Minimum Mapping Quality treshold for reads supporting InDel position. Default=20")
	Bam.add_argument('-iDP','--iEvaDepth',action="store_true",help="iEVA Depth for variant position. Only proper paired and proper mapped reads will be included.")
	Bam.add_argument('-iAD','--iAlleleDepth',action="store_true",help="iEVA Allele Depth reported as: Ref,Alt")
	Bam.add_argument('-RR','--ReadRef',action="store_true",help="Number of reads mapping REF allele.")
	Bam.add_argument('-RA','--ReadAlt',action="store_true",help="Number of reads mapping ALT allele.")
	Bam.add_argument('-QR','--MeanRefQscore',action="store_true",help="Mean Q-score for REF reads.")
	Bam.add_argument('-QA','--MeanAltQscore',action="store_true",help="Mean Q-score for ALT reads.")
	Bam.add_argument('-TDP','--TotalDPUnfilter',action="store_true",help="Total number of reads mapping position. No filter is applied. Look at -iDP option for filtered DP.")
	Bam.add_argument('-NCR','--NumberClippedReadsRef',action="store_true",help="Number of clipped reads mapping REF allele.")
	Bam.add_argument('-NCA','--NumberClippedReadsAlt',action="store_true",help="Number of clipped reads mapping ALT allele.")
	Bam.add_argument('-CR','--ClippedReadsRef',action="store_true",help="Percentage of clipped reads mapping REF allele.")
	Bam.add_argument('-CA','--ClippedReadsAlt',action="store_true",help="Percentage of clipped reads mapping ALT allele.")

	global opts

	opts = parser.parse_args()

	#Prima controllo il corretto inserimento delle opzioni:

	if opts.WindowSize in sys.argv:
		argument = sys.argv.index('-WS')
		Check_Parameters(sys.argv[argument],opts.WindowSize,20,600)


	out = open(opts.outfile,'w')
 	Reference = opts.reference
	bam_list = opts.list
	Header_File = []
	H_CHR = []


	#Estraggo dall'heder i campi FILTER, INFO e FORMAT
	with open(opts.input) as vcf:

		for line in vcf:

			line = line.rstrip()

			if line.startswith('##'):
				Header_File += [line]

			elif line.startswith('#CHROM'):
				H_CHR += [line]
				break


	#Da sys.argv estraggo le opzioni date nella command line. Per ognuna di queste, scrivo il campo da inserire nell'header
	for option in sys.argv:

		if '-SR' in option:
			Header_File += ['##INFO=<ID=SR,Number=1,Type=String,Description="Variant falls into repeated sequence. None=0, SimpleRepeat=1, Homopolymer=2.">']
		
		elif '-SRL' in option:
			Header_File += ['##INFO=<ID=SRL,Number=1,Type=Integer,Description="Length of repeat sequence for SR tag">']

		elif '-PNC' in option:
			Header_File += ['##INFO=<ID=PNC,Number=1,Type=Float,Description="Pseudo Nucleotide Composition using Kmer size of 2. Reported as: AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT">']
		
		elif '-RM' in option:
			Header_File += ['##INFO=<ID=RM,Number=1,Type=Integer,Description="Variant falls into repeated sequence computed with RepeatMasker tool. True=1, False=0">']
		
		elif '-PNC' in option:
			Header_File += ['##INFO=<ID=PNC,Number=1,Type=Float,Description="Describe nucleotide sequence in a window size of 40 using Pseudo Nucleotide Composition. Values reported as: ">']

		elif '-GC' in option:
			Header_File += ['##INFO=<ID=GC,Number=1,Type=Float,Description="GC content in sequence">']

		elif '-VC' in option:
			Header_File += ['##FORMAT=<ID=VC,Number=1,Type=String,Description="Annotate variant class: SNV=snv, Insertion=Ins, Deletion=Del, SequenceAlteration=Alt">']
		
		elif '-SBR' in option:
			Header_File += ['##FORMAT=<ID=SBR,Number=1,Type=Float,Description="Strand bias based on read orientation (R1+,R1-,R2+,R2-)">']
		
		elif '-AS' in option:
			Header_File += ['##FORMAT=<ID=AS,Number=1,Type=Float,Description="Reads mean alignment score">']
		
		elif '-UnMap' in option:
			Header_File += ['##FORMAT=<ID=UnMap,Number=1,Type=Float,Description="Percentage of unmapped reads">']

		elif '-MQ0' in option:
			Header_File += ['##FORMAT=<ID=MQ0,Number=1,Type=Float,Description="Percentage of reads having Mapping Quaility=0">']

		elif '-NPA' in option:
			Header_File += ['##FORMAT=<ID=NPA,Number=1,Type=Float,Description="Percentage of reads flagged as not primary alignment">']
		
		elif '-SA' in option:
			Header_File += ['##FORMAT=<ID=SA,Number=1,Type=Float,Description="Percentage of reads flagged as supplementary alignment">']
		
		elif '-NP' in option:
			Header_File += ['##FORMAT=<ID=NP,Number=1,Type=Float,Description="Percentage of not paired reads">']
		
		elif '-NPP' in option:
			Header_File += ['##FORMAT=<ID=NPP,Number=1,Type=Float,Description="Percentage of not proper paired reads">']

		elif '-NDT' in option:
			Header_File += ['##FORMAT=<ID=NDT,Number=1,Type=Integer,Description="Total number of duplicate reads">']

		elif '-NDR' in option:
			Header_File += ['##FORMAT=<ID=NDR,Number=1,Type=Integer,Description="Number of duplicate reads mapping REF allele">']

		elif '-NDA' in option:
			Header_File += ['##FORMAT=<ID=NDA,Number=1,Type=Integer,Description="Number of duplicate reads mapping ALT allele">']	
		
		elif '-DR' in option:
			Header_File += ['##FORMAT=<ID=DR,Number=1,Type=Float,Description="Percentage of REF reads marked as duplicates">']
		
		elif '-DA' in option:
			Header_File += ['##FORMAT=<ID=DA,Number=1,Type=Float,Description="Percentage of ALT reads marked as duplicates">']
		
		elif '-DDup' in option:
			Header_File += ['##FORMAT=<ID=DDup,Number=1,Type=Float,Description="Difference for duplicate reads in REF and ALT (DupREF-DupALT)">']
		
		elif '-iDP' in option:
			Header_File += ['##FORMAT=<ID=iDP,Number=1,Type=Integer,Description="iEVA read depth. Only proper paired and proper mapped reads will be included.">']

		elif '-iAD' in option:
			Header_File += ['##FORMAT=<ID=iAD,Number=R,Type=Integer,Description="Allelic depths reported by iEVA as Ref,Alt">']

		elif '-RR' in option:
			Header_File += ['##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Number of reads mapping REF allele">']
		
		elif '-RA' in option:
			Header_File += ['##FORMAT=<ID=RA,Number=1,Type=Integer,Description="Number of reads mapping ALT allele">']

		elif '-QR' in option:
			Header_File += ['##FORMAT=<ID=QR,Number=1,Type=Float,Description="Mean Q-score for reads supporting REF allele">']
		
		elif '-QA' in option:
			Header_File += ['##FORMAT=<ID=QA,Number=1,Type=Float,Description="Mean Q-score for reads supporting ALT allele">']
		
		elif '-TDP' in option:
			Header_File += ['##FORMAT=<ID=TDP,Number=1,Type=Integer,Description="Total read depth. No filter applied">']
		
		elif '-NCR' in option:
			Header_File += ['##FORMAT=<ID=NCR,Number=1,Type=Integer,Description="Number of clipped reads mapping REF allele">']
		
		elif '-NCA' in option:
			Header_File += ['##FORMAT=<ID=NCA,Number=1,Type=Integer,Description="Number of clipped reads mapping ALT allele">']

		elif '-ClipRef' in option:
			Header_File += ['##FORMAT=<ID=ClipRef,Number=1,Type=Float,Description="Percentage of clipped reads supporting REF">']
		
		elif '-ClipAlt' in option:
			Header_File += ['##FORMAT=<ID=ClipAlt,Number=1,Type=Float,Description="Percentage of clipped reads supporting ALT">']

		else:
			continue

	if opts.AlignmentScore or opts.StrandBiasReads or opts.UnMappedReads or opts.NotPrimaryAlignment \
	or opts.SupplementaryAlignment or opts.NotPairedReads or opts.NotProperPairedReads \
	or opts.DuplicateReference or opts.DuplicateAlternate or opts.DeltaDuplicate \
	or opts.ClippedReadsRef or opts.ClippedReadsAlt or opts.MeanRefQscore \
	or opts.MeanAltQscore or opts.ReadRef or opts.ReadAlt or opts.NumberReadDupRef or opts.NumberReadDupAlt \
	or opts.NumberTotalDupReads or opts.TotalDPUnfilter or opts.NumberClippedReadsRef\
	or opts.NumberClippedReadsAlt or opts.iEvaDepth or opts.iAlleleDepth or opts.MappingQualityZero:
		if opts.list:
			pass
		else:
			sys.exit('To Enable bam features insert path to bam list with command -L (--list)')

	#Scrivo in output la nuova header
	out.write('\n'.join(Header_File) + '\n' + '\n'.join(H_CHR) + '\n')

	#Splitto l'header del vcf
	H_CHR = H_CHR[0].split('\t')

	#Estraggo dall'header la lista dei campioni:

	if opts.list:
		Sample_list = H_CHR[H_CHR.index('FORMAT')+1:]
		Sample_dict = Check_Bam(Sample_list, bam_list)

	#Ora passiamo ad annotare le varianti nel vcf:

	CHR_Counter = ''

	with open(opts.input) as vcf:

		for variant in vcf:

			if variant.startswith('#'):

				continue

			else:

				variant = variant.rstrip()
				variant = variant.split('\t')

				Variant_Class = Extract_Variant_Type(variant, H_CHR)


				if opts.SimpleRepeat or opts.SimpleRepeatLength or opts.PseudoNucleotidesComposition:
					RepeatSeq, RepeatSeq_Lenght, Psuedo_Nucleotide = SimpleRepeats_Finder(Reference, variant, H_CHR, Variant_Class, opts.WindowSize)

				if opts.SimpleRepeat:
					variant[H_CHR.index('INFO')] += ';' + 'SR=' + str(RepeatSeq)

				if opts.SimpleRepeatLength:
					variant[H_CHR.index('INFO')] += ';' + 'SRL=' + str(RepeatSeq_Lenght)

				if opts.PseudoNucleotidesComposition:
					variant[H_CHR.index('INFO')] += ';' + 'PNC=' + ','.join(str(Din) for Din in Psuedo_Nucleotide)

				if opts.RepeatMasker:
					RM = Check_Repeat_Masker(Reference, variant, H_CHR)
					variant[H_CHR.index('INFO')] += ';' + 'RM=' + RM

				if opts.gcContent:
					GC = GC_content(Reference, variant, H_CHR, WindowSize)
					variant[H_CHR.index('INFO')] += ';' + 'GC=' + GC

				if opts.VariantClass:
					variant[H_CHR.index('INFO')] += ';' + 'VC=' + str(Variant_Class)

				if opts.AlignmentScore or opts.StrandBiasReads or opts.UnMappedReads or opts.NotPrimaryAlignment \
				or opts.SupplementaryAlignment or opts.NotPairedReads or opts.NotProperPairedReads \
				or opts.DuplicateReference or opts.DuplicateAlternate or opts.DeltaDuplicate \
				or opts.ClippedReadsRef or opts.ClippedReadsAlt or opts.MeanRefQscore \
				or opts.MeanAltQscore or opts.ReadRef or opts.ReadAlt or opts.NumberReadDupRef or opts.NumberReadDupAlt \
				or opts.NumberTotalDupReads or opts.TotalDPUnfilter or opts.NumberClippedReadsRef\
				or opts.NumberClippedReadsAlt or opts.iEvaDepth or opts.iAlleleDepth or opts.MappingQualityZero:
					Sample_Stat = Sequence_Annotator(variant, Sample_dict, H_CHR, Reference, Variant_Class, opts.SNVMinBaseQuality, opts.SNVMinMappingQuality, opts.IndelMinBaseQuality, opts.IndelMinMappingQuality)

				if opts.StrandBiasReads:
					variant[H_CHR.index('FORMAT')] += ':' + 'SBR'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Strand_Bias_Reads')))

				if opts.UnMappedReads:
					variant[H_CHR.index('FORMAT')] += ':' + 'UnMap'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Percentage_Unmapped_Reads')))

				if opts.MappingQualityZero:
					variant[H_CHR.index('FORMAT')] += ':' + 'MQ0'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Mapping_Quality_Zero')))

				if opts.NotPrimaryAlignment:
					variant[H_CHR.index('FORMAT')] += ':' + 'NPA'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Percentage_Not_Primary_Alignment_Reads')))

				if opts.SupplementaryAlignment:
					variant[H_CHR.index('FORMAT')] += ':' + 'SA'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Percentage_Supplementary_Align')))

				if opts.NotPairedReads:
					variant[H_CHR.index('FORMAT')] += ':' + 'NP'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Percentage_Not_Paired_Reads')))

				if opts.NotProperPairedReads:
					variant[H_CHR.index('FORMAT')] += ':' + 'NPP'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Percentage_Not_Proper_Paired_Reads')))

				if opts.AlignmentScore:
					variant[H_CHR.index('FORMAT')] += ':' + 'AS'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Mean_Alignment_Score')))

				if opts.NumberTotalDupReads:
					variant[H_CHR.index('FORMAT')] += ':' + 'NDT'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Total_Duplicate_Reads'))

				if opts.NumberReadDupRef:
					variant[H_CHR.index('FORMAT')] += ':' + 'NDR'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Number_Read_Dup_Ref'))

				if opts.NumberReadDupAlt:
					variant[H_CHR.index('FORMAT')] += ':' + 'NDA'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Number_Read_Dup_Alt'))

				if opts.DuplicateReference:
					variant[H_CHR.index('FORMAT')] += ':' + 'DR'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Percentage_Dup_Ref')))

				if opts.DuplicateAlternate:
					variant[H_CHR.index('FORMAT')] += ':' + 'DA'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Percentage_Dup_Alt')))

				if opts.DeltaDuplicate:
					variant[H_CHR.index('FORMAT')] += ':' + 'DDup'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Delta_Dupicate')))

				if opts.iEvaDepth:
					variant[H_CHR.index('FORMAT')] += ':' + 'iDP'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Coverage'))

				if opts.iAlleleDepth:
					variant[H_CHR.index('FORMAT')] += ':' + 'iAD'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Read_Ref')) + ',' + str(Sample_Stat.get(sample).get('Read_Alt'))

				if opts.ReadRef:
					variant[H_CHR.index('FORMAT')] += ':' + 'RR'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Read_Ref'))

				if opts.ReadAlt:
					variant[H_CHR.index('FORMAT')] += ':' + 'RA'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Read_Alt'))

				if opts.MeanRefQscore:
					variant[H_CHR.index('FORMAT')] += ':' + 'QR'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Qual_Ref'))

				if opts.MeanAltQscore:
					variant[H_CHR.index('FORMAT')] += ':' + 'QA'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Qual_Alt'))

				if opts.TotalDPUnfilter:
					variant[H_CHR.index('FORMAT')] += ':' + 'TDP'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Total_Reads_Unfilter'))

				if opts.NumberClippedReadsRef:
					variant[H_CHR.index('FORMAT')] += ':' + 'NCR'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Clipped_Reads_Ref'))

				if opts.NumberClippedReadsAlt:
					variant[H_CHR.index('FORMAT')] += ':' + 'NCA'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Clipped_Reads_Alt'))

				if opts.ClippedReadsRef:
					variant[H_CHR.index('FORMAT')] += ':' + 'CR'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Percentage_ClipRef')))

				if opts.ClippedReadsAlt:
					variant[H_CHR.index('FORMAT')] += ':' + 'CA'
					for sample in Sample_dict.keys():
						variant[H_CHR.index(sample)] += ':' + str(Check_Zero(Sample_Stat.get(sample).get('Percentage_ClipAlt')))

				if variant[H_CHR.index('#CHROM')] != CHR_Counter and CHR_Counter == '':
					CHR_Counter = variant[H_CHR.index('#CHROM')]
					print '\n' + 'Extracting attributes on: ' + str(CHR_Counter)

				if variant[H_CHR.index('#CHROM')] != CHR_Counter and CHR_Counter != '':
					CHR_Counter = variant[H_CHR.index('#CHROM')]
					print '\n' + 'Extracting attributes on: ' + str(CHR_Counter)


				out.write('\t'.join(variant) + '\n')

	out.close()

	print '\n' + 'Exatraction: Done' + '\n'

main()
