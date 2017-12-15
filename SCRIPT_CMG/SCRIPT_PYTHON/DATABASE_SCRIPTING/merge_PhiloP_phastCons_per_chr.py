import argparse
import os
import glob2


def main():
	parser = argparse.ArgumentParser('Adds number of exon for each interval')
	
	parser.add_argument('-p','--path',help="path dei tsv")
	parser.add_argument('-o','--out',help="path di output")

	global opts
	opts = parser.parse_args()

	
	

	for filename in glob2.glob(os.path.join(opts.path,'*placental*')):

		nome_no_path=filename.split('/')[-1]
		nome_no_ext=nome_no_path.split('.')[0]
		nome_db=nome_no_ext.split('_')[0]
		tipo=nome_no_ext.split('_')[1]
		chr=nome_no_ext.split('_')[2]
		pannello=nome_no_ext.split('_')[3]

		filename_primate= '_'.join([nome_db,'primate',chr,pannello]) + '.chr.tsv'
		filename_vertebrate= '_'.join([nome_db,'vertebrate',chr,pannello]) + '.chr.tsv'

		posizioni=dict()

		placental=open(filename,'r')
		print "apro",filename 
		for line in placental:
			line=line.rstrip()
			if line.startswith('CHROM'):
				continue
			else:
				chr=line.split('\t')[0]
				pos=line.split('\t')[1]
				score=line.split('\t')[2]

				posizioni[chr+'\t'+ pos]=[score,'.','.']

		placental.close()

		#print posizioni

		primate=open(opts.path + '/'+filename_primate,'r')
		print "apro",opts.path + '/'+filename_primate
		
		for line in primate:
			line=line.rstrip()
			if line.startswith('CHROM'):
				continue
			else:
				#print line
				chr=line.split('\t')[0]
				pos=line.split('\t')[1]
				score=line.split('\t')[2]
				print line
				if chr+'\t'+ pos in posizioni.keys():
					scores=posizioni.get(chr+'\t'+ pos)
					scores[1]=score
					posizioni[chr+'\t'+ pos]=scores
				else:
					posizioni[chr+'\t'+ pos]=['.',score,'.']

		primate.close()

		vertebrate=open(opts.path + '/'+filename_vertebrate,'r')
		print "apro",opts.path + '/'+filename_vertebrate
		for line in vertebrate:
			line=line.rstrip()
			if line.startswith('CHROM'):
				continue
			else:
				chr=line.split('\t')[0]
				pos=line.split('\t')[1]
				score=line.split('\t')[2]
				print line
				if chr+'\t'+ pos in posizioni.keys():
					scores=posizioni.get(chr+'\t'+ pos)
					scores[2]=score
					posizioni[chr+'\t'+ pos]=scores
				else:
					posizioni[chr+'\t'+ pos]=['.','.',score]
		vertebrate.close()

		out = open(opts.out + '/' + '_'.join([nome_db,chr,pannello]) + '.tsv','w')

		print 'scrivo',opts.out + '/' + '_'.join([nome_db,chr,pannello]) + '.tsv'
		out.write('CHROM'+'\t' + 'POS'+'\t' + nome_db + '_placental'+'\t' + nome_db + '_primate' +'\t' +nome_db + '_vertebrate'+'\n')
	
		for var in posizioni.keys():

			out.write(var +'\t' + '\t'.join(posizioni.get(var))+'\n')

		out.close()

main()