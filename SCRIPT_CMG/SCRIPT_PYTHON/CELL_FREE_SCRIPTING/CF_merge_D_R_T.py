import argparse


def read_snp(varianti,somT):
	lista_snp = open(opts.snp_list,'r')
	for line in lista_snp:
		line=line.rstrip()
		if line.startswith('#') or  line.startswith('CHROM'):
			continue
		else:
			chrom = line.split('\t')[0]
			pos = line.split('\t')[1]
			var=['','','']
			for i in somT:
				var+=['']
			#print chrom + '\t'+ pos
			varianti[chrom + '\t'+ pos]=var
	lista_snp.close()

def read_donor(varianti,varianti_other,header,somT):
	lista_snp = open(opts.donor,'r')
	for line in lista_snp:
		line=line.rstrip()
		if line.startswith('CHROM'):
			HEADER=line.split('\t')
			header = header + '\t' +'\t'.join(['GT_GATK_D','GT_VARSCAN_D','GT_FREEBAYES_D','DP_D','AF_D','QB_D'])
		else:
			chrom = line.split('\t')[0]
			pos = line.split('\t')[1]
			ref = line.split('\t')[3]
			alt = line.split('\t')[4]
			gt = '\t'.join(line.split('\t')[HEADER.index('GT_GATK'):HEADER.index('GT_Freebayes') +1])
			#print gt
			#dp = '\t'.join(line.split('\t')[11:14])
			dp = line.split('\t')[HEADER.index('DP_media')]
			qb = line.split('\t')[HEADER.index('MBQ_media')]
			AF = line.split('\t')[HEADER.index('AF_media')]
			donor = '\t'.join([gt,dp,AF,qb])
			#print chrom + '\t'+ pos,donor, varianti.has_key(chrom + '\t'+ pos)
			#print chrom,pos
			if varianti.has_key(chrom + '\t'+ pos):
				#print "SI" 				
				var = varianti.get(chrom + '\t'+ pos)
				if var[0]=='':
					var[0]='.\t'+ref+'\t'+alt
				var[1] = donor
				varianti[chrom + '\t'+ pos]=var
			elif varianti_other.has_key(chrom + '\t'+ pos):
				donor = '\t'.join([gt,dp,AF,qb])
				var = varianti_other.get(chrom + '\t'+ pos)
				if var[0]=='':
					var[0]='.\t'+ref+'\t'+alt
				var[1] = donor
				varianti_other[chrom + '\t'+ pos]=var
			else:
				donor = '\t'.join([gt,dp,AF,qb])				
				var = ['.\t'+ref+'\t'+alt,donor,'']
				for i in somT:
					var+=['']
				varianti_other[chrom + '\t'+ pos]=var



	lista_snp.close()
	return header

def read_rec(varianti,varianti_other,header,somT):
	lista_snp = open(opts.recipient,'r')
	for line in lista_snp:
		line=line.rstrip()
		if line.startswith('CHROM'):
			HEADER=line.split('\t')
			header = header + '\t' + '\t'.join(['GT_GATK_R','GT_VARSCAN_R','GT_FREEBAYES_R','DP_R','AF_R','QB_R'])
		else:
			chrom = line.split('\t')[0]
			pos = line.split('\t')[1]
			ref = line.split('\t')[3]
			alt = line.split('\t')[4]
			gt = '\t'.join(line.split('\t')[HEADER.index('GT_GATK'):HEADER.index('GT_Freebayes') +1])
			#print gt
			#dp = '\t'.join(line.split('\t')[11:14])
			dp = line.split('\t')[HEADER.index('DP_media')]
			qb = line.split('\t')[HEADER.index('MBQ_media')]
			AF = line.split('\t')[HEADER.index('AF_media')]
			if varianti.has_key(chrom + '\t'+ pos):
				rec = '\t'.join([gt,dp,AF,qb])
				var = varianti.get(chrom + '\t'+ pos)
				if var[0]=='':
					var[0]='.\t'+ref+'\t'+alt
				var[2]= rec
				varianti[chrom + '\t'+ pos]=var
			elif varianti_other.has_key(chrom + '\t'+ pos):
				rec = '\t'.join([gt,dp,AF,qb])
				var = varianti_other.get(chrom + '\t'+ pos)
				if var[0]=='':
					var[0]='.\t'+ref+'\t'+alt
				var[2]= rec
				varianti_other[chrom + '\t'+ pos]=var
			else:
				rec = '\t'.join([gt,dp,AF,qb])				
				var = ['.\t'+ref+'\t'+alt,'',rec]
				for i in somT:
					var+=['']
				varianti_other[chrom + '\t'+ pos]=var

	lista_snp.close()
	return header

def read_sT(varianti,varianti_other,header,index,st,somT):
	lista_snp = open(st,'r')
	for line in lista_snp:
		line=line.rstrip()

		if line.startswith('CHROM'):
			HEADER=line.split('\t')
			header = header +  '\t' +'\t'.join(['SOMATIC_VARSCAN_T'+ str(index+1),'FILTER_MUTECT_T'+ str(index+1),'STATUS_VARDICT_T'+ str(index+1),'GT_CF_MUTECT_T'+ str(index+1),'GT_CF_VARSCAN_T'+ str(index+1),'GT_CF_VARDICT_T'+ str(index+1),
										'GT_G_MUTECT_T'+ str(index+1),'GT_G_VARSCAN_T'+ str(index+1),'GT_G_VARDICT_T'+ str(index+1),
										'DP_CF_T'+ str(index+1),'AF_CF_T'+ str(index+1),'QB_CF_T'+ str(index+1),'DP_G_T'+ str(index+1),'AF_G_T'+ str(index+1),'QB_G_T'+ str(index+1),'DELTA_MEDIA_T'+ str(index+1)])
		else:
			chrom = line.split('\t')[0]
			pos = line.split('\t')[1]
			ref = line.split('\t')[3]
			alt = line.split('\t')[4]
			somatic_varscan = line.split('\t')[HEADER.index('SomaticVarscan')]
			filter_mutect = line.split('\t')[HEADER.index('FILTER_Mutect')]
			status_vardict = line.split('\t')[HEADER.index('STATUS_Vardict')]

			gt_t_varscan = line.split('\t')[HEADER.index('GT_t_Varscan')]
			gt_t_vardict = line.split('\t')[HEADER.index('GT_t_Vardict')]
			gt_t_mutect = line.split('\t')[HEADER.index('GT_t_Mutect')]

			gt_n_varscan = line.split('\t')[HEADER.index('GT_n_Varscan')]
			gt_n_vardict = line.split('\t')[HEADER.index('GT_n_Vardict')]
			gt_n_mutect = line.split('\t')[HEADER.index('GT_n_Mutect')]
			dp_t = line.split('\t')[HEADER.index('DP_tum_media')]
			af_t = line.split('\t')[HEADER.index('AF_tum_media')]
			qb_t = line.split('\t')[HEADER.index('MBQT_media')]
			dp_n = line.split('\t')[HEADER.index('DP_norm_media')]
			af_n = line.split('\t')[HEADER.index('AF_norm_media')]
			qb_n = line.split('\t')[HEADER.index('MBQN_media')]
			try:
				delta = line.split('\t')[HEADER.index('Delta_media')]
			except:
				delta = line.split('\t')[HEADER.index('Delta_median')]
			if varianti.has_key(chrom + '\t'+ pos):
				t = '\t'.join([somatic_varscan,filter_mutect,status_vardict,gt_t_mutect,gt_t_varscan,gt_t_vardict,gt_n_mutect,gt_n_varscan,gt_n_vardict,dp_t,af_t,qb_t,dp_n,af_n,qb_n,delta])
				var = varianti.get(chrom + '\t'+ pos)
				if var[0]=='':
					var[0]='.\t'+ref+'\t'+alt
				var[3+index]= t
				varianti[chrom + '\t'+ pos]=var
			elif varianti_other.has_key(chrom + '\t'+ pos):
				t = '\t'.join([somatic_varscan,filter_mutect,status_vardict,gt_t_mutect,gt_t_varscan,gt_t_vardict,gt_n_mutect,gt_n_varscan,gt_n_vardict,dp_t,af_t,qb_t,dp_n,af_n,qb_n,delta])
				var = varianti_other.get(chrom + '\t'+ pos)
				if var[0]=='':
					var[0]='.\t'+ref+'\t'+alt
				var[3+index]= t
				varianti_other[chrom + '\t'+ pos]=var
			else:
				t = '\t'.join([somatic_varscan,filter_mutect,status_vardict,gt_t_mutect,gt_t_varscan,gt_t_vardict,gt_n_mutect,gt_n_varscan,gt_n_vardict,dp_t,af_t,qb_t,dp_n,af_n,qb_n,delta])			
				var = ['.\t'+ref+'\t'+alt,'','']
				for i in somT:
					var += ['']
				var[3+index]=t	
				varianti_other[chrom + '\t'+ pos]=var

	lista_snp.close()
	return header

def check(varianti):
	for variante in varianti.keys():
		var = varianti.get(variante)
		if var[0]=='':
			var[0] = '\t'.join(['.','.','.'])
		if var[1]=='':
			var[1] = '\t'.join(['0/0','0/0','0/0','.','0','.'])
		if var[2]=='':
			var[2] = '\t'.join(['0/0','0/0','0/0','.','0','.'])
		for el in var[3:]:
			if el=='':
				var[var.index(el)]='\t'.join(['.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.'])

		if not opts.recipient:
			del var[2]
		if not opts.donor:
			del var[1]



def print_varianti(varianti,out,header,vett):
	out.write(header+'\n')
	for variante in varianti.keys():
		var = varianti.get(variante)
		#print var
		if opts.ctrl_snp_list:
			if variante in vett:
				continue	
			else:
				out.write('\t'.join([variante] + var) + '\n')
		else:
			out.write('\t'.join([variante] + var) + '\n')




def main():

	parser = argparse.ArgumentParser('Parse VCF output from Variant callers to output a variant_dataset.txt.')
	parser.add_argument('-d', '--donor', help="donor.tsv",default=None)
	parser.add_argument('-r', '--recipient', help="recipient.tsv",default=None)
	parser.add_argument('--sT', help="lista di sample.tsv somatici Don_Rec divisi da ';'")
	parser.add_argument('-o', '--out',help="path di output")
	parser.add_argument('-l', '--snp_list',help="lista di snp da matchare",default=None)
	parser.add_argument('-i', '--ctrl_snp_list',help="lista di snp da escludere",default=None)


	global opts 
	opts = parser.parse_args()
	vett=[]

	if opts.ctrl_snp_list:
		ctrl=open(opts.ctrl_snp_list,'r')
		
		for line in ctrl:
			line=line.rstrip()
			if line.startswith('#') or line.startswith('CHROM'):
				continue
			else:
				vett+=['\t'.join(line.split('\t')[0:2])] 


	varianti = dict()
	varianti_other = dict()
	somT=opts.sT.split(';')

	header='CHROM\tPOS\tID\tREF\tALT'
	if opts.snp_list:
		read_snp(varianti,somT)
	if opts.donor:
		header=read_donor(varianti,varianti_other,header,somT)
	if opts.recipient:
		header=read_rec(varianti,varianti_other,header,somT)
	
	for st in somT:
		index=somT.index(st)
		header=read_sT(varianti,varianti_other,header,index,st,somT)
	check(varianti)
	check(varianti_other)

	out_snp_list=open(opts.out+'.list.tsv','w')
	#out_other=open(opts.out+'.other.tsv','w')
	print_varianti(varianti,out_snp_list,header,vett)
	#print_varianti(varianti_other,out_other,header,vett)

	
main()
