import argparse

def extract_var_from_db(varianti,uniprot,srepeats):
	outfile=open(opts.out,'w')

	for var in varianti:
		var=var.rstrip()
		uniprot_type='-'
		uniprot_disease='-'
		repetition_period='-'
		repetition_num='-'
		repetition_perc_match='-'
		repetition_perc_indel='-'
		repetition_score='-'
		repetition_entropy='-'
		repetition_isHomo='-'
		if var.startswith('CHROM'):
			header=var.split('\t')
			uniprot_type='Uniprot_type'
			uniprot_disease='Uniprot_Disease'
			repetition_isHomo='HOMOPOLIMER'
			repetition_period='RepPeriod'
			repetition_num='RepNum'
			repetition_perc_match='RepPercMatch'
			repetition_perc_indel='RepPercIndel'
			repetition_score='RepScore'
			repetition_entropy='RepEntropy'
		else:
			chrom=var.split('\t')[header.index('CHROM')]
			pos=var.split('\t')[header.index('POS')]
			ref=var.split('\t')[header.index('REF')]
			alt=var.split('\t')[header.index('ALT')]
			HGVSp=var.split('\t')[header.index('HGVSp')]
			HGVSc=var.split('\t')[header.index('HGVSc')]
			gene=var.split('\t')[header.index('SYMBOL')]

			for line in uniprot:
				if gene in line and HGVSp in line:
					uniprot_type=line.split('\t')[4]
					uniprot_deasase=line.split('\t')[-1]

			for line in srepeats:
				line=line.strip()
				if line.startswith('#'):
					header_sr=line.split('\t')
				else:
					chrom_sr=line.split('\t')[header_sr.index('chrom')]
					chromStart_sr=line.split('\t')[header_sr.index('chromStart')]
					chromEnd_sr=line.split('\t')[header_sr.index('chromEnd')]
					name_sr=line.split('\t')[header_sr.index('name')]
					period_sr=line.split('\t')[header_sr.index('period')]
					copyNum_sr=line.split('\t')[header_sr.index('copyNum')]
					consensusSize_sr=line.split('\t')[header_sr.index('consensusSize')]
					perMatch_sr=line.split('\t')[header_sr.index('perMatch')]
					perIndel_sr=line.split('\t')[header_sr.index('perIndel')]
					score_sr=line.split('\t')[header_sr.index('score')]
					A_sr=line.split('\t')[header_sr.index('A')]
					C_sr=line.split('\t')[header_sr.index('C')]
					G_sr=line.split('\t')[header_sr.index('G')]
					T_sr=line.split('\t')[header_sr.index('T')]
					entropy_sr=line.split('\t')[header_sr.index('entropy')]
					sequence_sr=line.split('\t')[header_sr.index('sequence')]

					if chrom == chrom_sr and int(pos)>int(chromStart_sr) and  int(pos)<int(chromEnd_sr):
						repetition_period=period_sr
						repetition_num=copyNum_sr
						repetition_perc_match=perMatch_sr
						repetition_perc_indel=perIndel_sr
						repetition_score=score_sr
						repetition_entropy=entropy_sr
						repetition_isHomo='HOMOPOLIMER'
		
		info_db='\t'.join([uniprot_type,uniprot_disease,repetition_isHomo,repetition_period,
			repetition_num,repetition_perc_match,repetition_perc_indel,repetition_score,repetition_entropy])
		outfile.write(var+'\t'+info_db+'\n')		


def main():
	parser = argparse.ArgumentParser('Aggiunge informazioni su omopolimeri e polimorfismi da database')
	parser.add_argument('-i','--var',help="file di varianti annotate tab delimited")
	parser.add_argument('-u','--uniprot_db',help="database uniprot di polimorfismi")
	parser.add_argument('-s','--simple_repeats_db',help="database di simple repeats")
	parser.add_argument('-o','--out',help="file in output")
	
	global opts
	opts = parser.parse_args()
	varianti=open(opts.var,'r')
	uniprot=open(opts.uniprot_db,'r')
	srepeats=open(opts.simple_repeats_db,'r')
	extract_var_from_db(varianti,uniprot,srepeats)
	

main()