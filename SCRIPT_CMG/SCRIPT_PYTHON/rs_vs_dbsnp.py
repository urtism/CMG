
import argparse


if __name__ == '__main__':

	parser = argparse.ArgumentParser('incrocia dbsnp con la lista di rs per estrarre ref e alt')
	parser.add_argument('--dbsnp', help="path a dbsnp",default=None)
	parser.add_argument('--rs_list', help="lista di rs in vcf format",default=None)
	parser.add_argument('--out', help="file di out in vcf format",default=None)

	global opts
	opts = parser.parse_args()
	vcfout = open(opts.out,'w')
	rs=dict()
	with open(opts.rs_list,'r') as vcf:
		for line in vcf:
			line = line.rstrip()
			if line.startswith('#'):
				vcfout.write(line+'\n')
			else:
				#chr,pos,id,ref,alt,filter,info,format,sample = line.split('\t')
				id = line

				rs[id]=[]



	dbsnp = open(opts.dbsnp,'r')
	for row in dbsnp:
		if row.startswith('#'):
			continue
		else:
			if str(row.split('\t')[2]) in rs.keys():
				print row.split('\t')[2]
				chr = row.split('\t')[0]
				pos = row.split('\t')[1]
				ref = row.split('\t')[3]
				alt = row.split('\t')[4]
				vcfout.write('\t'.join([chr,pos,id,ref,alt,'.','.','.','.'])+'\n')
				
	dbsnp.close()

						
							





