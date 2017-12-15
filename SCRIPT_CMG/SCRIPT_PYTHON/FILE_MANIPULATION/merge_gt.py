import argparse
from collections import Counter

if __name__ == '__main__':

	parser = argparse.ArgumentParser('Confronta i GT dei tre variant callers e restituisce un GT unico')
	parser.add_argument('-i','--tsv',help="file tsv in ingresso")
	parser.add_argument('-a','--analisi',help="Germline o Somatic")
	parser.add_argument('-o','--out',help="file di out")

	global opts
	opts = parser.parse_args()
	in_tsv=open(opts.tsv,'r')
	out_tsv=open(opts.out,'w')
	if opts.analisi == 'Germline':
		for line in in_tsv:
			line_split=line.rstrip().split('\t')
			if line.startswith('chr'):
				gts=[line_split[header.index('GT_GATK')],line_split[header.index('GT_Varscan')],line_split[header.index('GT_Freebayes')]]
				hom=0
				het=0
				for gt in gts:
					if gt=='1/1':
						hom+=1
					elif gt=='0/1' or gt=='1/0':
						het+=1
				if hom>=het:
					gt='1/1'
				else:
					gt='0/1'
				line_split[header.index('GT_Freebayes')+1:header.index('GT_Freebayes')+1]=[gt]
				
			else:
				header=line_split
				line_split[header.index('GT_Freebayes')+1:header.index('GT_Freebayes')+1]=['GT']

			out_tsv.write('\t'.join(line_split)+'\n')
