import argparse



if __name__ == '__main__':

	parser = argparse.ArgumentParser('Incrocia varianti in formato tsv con il database Cosmic e aggiunge al file le informazioni estratte dal db')
	parser.add_argument('-i','--input',help="file tsv annotato con VEP")
	parser.add_argument('-c','--cosmic',help="path del database cosmic in formato vcf o tsv.")
	parser.add_argument('--organo',help="organo da ricercare in cosmic",default='liver')
	parser.add_argument('-o','--outfile',help="file di output in tsv format")
	global opts
	opts = parser.parse_args()

	cosmic=open(opts.cosmic,'r').readlines()
	vars=open(opts.input,'r')
	out=open(opts.outfile,'w')

	if cosmic[0].startswith('##fileformat=VCFv4.1'):
		cosm_vcf=True

	for line in vars:
		line=line.rstrip()
		if line.startswith('CHROM'):
			header=line.split()
			out.write('\t'.join(header+['ORGAN'])+'\n')
		else:
			chr,pos,id,ref,alt=line.split()[:5]
			id='.'
			cosm_id=line.split()[header.index('Existing_variation')].split('&')
			for i in cosm_id:
				if not i.startswith('COSM'):
					del cosm_id[cosm_id.index(i)]
			found=False
			for var in cosmic[2:]:
				if cosm_vcf:
					for i in cosm_id:
						if i +',' in var.split('\t')[7] and opts.organo in var.split('\t')[7]: 
							#out.write('\t'.join([chr,pos,id,ref,alt,var.split('\t')[7]])+'\n')
							out.write('\t'.join([line,opts.organo])+'\n')
							found=True
							break
						elif i +';' in var.split('\t')[7] and opts.organo in var.split('\t')[7]: 
							#out.write('\t'.join([chr,pos,id,ref,alt,var.split('\t')[7]])+'\n')
							out.write('\t'.join([line,opts.organo])+'\n')
							found=True
							break
						elif [chr,pos] == var.split('\t')[:1] and opts.organo in var.split('\t')[7]:
							#out.write('\t'.join([chr,pos,id,ref,alt,var.split('\t')[7]])+'\n')
							out.write('\t'.join([line,opts.organo])+'\n')
							found=True
							break
					if found:
						break
			if not found:
				out.write('\t'.join([line,'.'])+'\n')
					# if [chr,pos] == var.split('\t')[:1] and 'liver' in var.split('\t')[7]:
					# 	out.write('\t'.join([chr,pos,id,ref,alt,var.split('\t')[7]])+'\n')
					# else:
					# 	for i in cosm_id:
					# 		if i in var.split('\t')[7] and 'liver' in var.split('\t')[7]:
					# 			out.write('\t'.join([chr,pos,id,ref,alt,var.split('\t')[7]])+'\n')	



