import argparse


if __name__ == '__main__':

	parser = argparse.ArgumentParser('Genera un report delle varianti trovate in un campione o in un gruppo di campioni')
	parser.add_argument('--tsvlist',help="lista di file tsv di varianti da ricercare nella finestra data",default=None)
	parser.add_argument('--vcflist',help="lista di file vcf di varianti da ricercare nella finestra data",default=None)
	parser.add_argument('-o','--out',help="path di out in cui stampare le statistiche")
	parser.add_argument('--varlist',help="file di varianti da usare per cercare le varianti nei tsv")
	parser.add_argument('-w','--window',help="finestra da applciare alle varianti di list per cercare le varianti nei tsv o vcf",default=500)
	
	global opts
	opts = parser.parse_args()

	varlist = open(opts.varlist,'r')
	out = open(opts.out,'w')
	VARLIST = dict()

	for varline in varlist:
		varline=varline.rstrip()
		if varline.startswith('Nome'):
			varheader=varline.split('\t')
		else:
			varnome,varid,varchr,varpos,varref,varalt=varline.split('\t')[:6]
			#gene = line.split('\t')[header.index('SYMBOL')]
			varHGVSc = varline.split('\t')[varheader.index('HGVSc')]
			varHGVSp = varline.split('\t')[varheader.index('HGVSp')]
			var = '-'.join([varid,varchr,varpos,varref,varalt,varHGVSc])
			if var in VARLIST.keys():
				print varnome,varid,varchr,varpos,'gia fatto'
				continue
			else:
				print varnome,varid,varchr,varpos
				VARLIST[var]=[]
				if opts.tsvlist is not None:
					tsvlist = open(opts.tsvlist,'r')
					for path in tsvlist:
						path=path.rstrip()
						if varid in path: 
							if 'GATK' in path:
								variant_caller = 'GATK' 
							elif 'FREEBAYES' in path:
								variant_caller = 'FREEBAYES'
							elif 'VarScan' in path:
								variant_caller = 'VarScan'

							variant_caller
							tsv=open(path,'r')
							for line in tsv:
								line=line.rstrip()
								if line.startswith('CHROM'):
									header=line.split('\t')
								else:
									chr,pos,id,ref,alt=line.split('\t')[:5]
									gene = line.split('\t')[header.index('SYMBOL')]
									HGVSc = line.split('\t')[header.index('HGVSc')]
									HGVSp = line.split('\t')[header.index('HGVSp')]
									if chr == varchr and int(pos) <= int(varpos) + opts.window and int(pos) >= int(varpos) - opts.window:
										VARLIST[var]+=['\t'.join([id,chr,pos,ref,alt,gene,HGVSc,HGVSp,variant_caller])]
							tsv.close()
					tsvlist.close()
				out.write(var+' :\n')
				for vars in VARLIST[var]:
					out.write(vars +'\n')
				out.write('\n')

	varlist.close()
	

	# for var in VARLIST.keys():
	# 	out.write(var+' :\n')
	# 	for vars in VARLIST[var]:
	# 		out.write(vars +'\n')
	# 	out.write('\n')



