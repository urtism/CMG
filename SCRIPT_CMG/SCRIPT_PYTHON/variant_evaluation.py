import argparse
import operator


def evaluate(e_dataset,c_dataset):
	tp=0.0
	fp=0.0
	tn=0.0
	fn=0.0
	gt_match=0.0
	gt_no_match=0.0
	enum=len(e_dataset.keys())
	cnum=len(c_dataset.keys())
	
	for var in e_dataset.keys():
		if var in c_dataset.keys():
			tp+=1
			if e_dataset[var] == c_dataset[var]:
				gt_match+=1
			else:
				gt_no_match+=1
			del(e_dataset[var])
			del(c_dataset[var])
	
	fp=len(e_dataset.keys())
	fn=len(c_dataset.keys())
	return tp,fp,tn,fn,enum,cnum,gt_match,c_dataset


def calc_stats(tp,fp,tn,fn,enum,cnum,gt_match):
	f_disc_rate=fp/float(enum)
	tp_rate=tp/(tp+fn)
	precision=tp/(tp+fp)
	detect_rate=tp/float(cnum)
	missing_rate=fn/float(cnum)
	gt_match_rate=gt_match/tp
	acc=(tp+tn)/(tp+fp+fn+tn)
	scoref1=2*tp/(2*tp+fp+fn)
	return f_disc_rate,tp_rate,precision,detect_rate,missing_rate,gt_match_rate,acc,scoref1

def print_stats(type,f_disc_rate,tp_rate,precision,detect_rate,missing_rate,gt_match_rate,acc,scoref1,outfile,dataset,enum,cnum):
	outfile.write('\t'.join([type,str(enum),str(cnum),str(detect_rate),str(f_disc_rate),str(tp_rate),str(precision),str(acc),str(scoref1),str(missing_rate),str(gt_match_rate)])+'\n')
	if dataset != None:
		outfile.write('\nVariants missed:'+'\n')
		sorted_dataset = sorted(dataset.items(), key=operator.itemgetter(0))
		#print sorted_dataset
		for var in sorted_dataset:
			outfile.write(var[0]+'\n')


def print_header(outfile):
	outfile.write('\t'.join(['VARIANTI','NUM CALLED VARIANT','NUM REAL VARIANTS','DETECTION RATE','FALSE DISCOVER RATE','SENSITIVITY (TP RATE)','PRECISION','ACCURACY','SCORE F1','MISSING RATE','GT MATCH RATE'])+'\n')

if __name__ == '__main__':

	parser = argparse.ArgumentParser('Confronta un file con il benchmark e calcola le statistiche per la valutazione del file in entrate')
	parser.add_argument('-e','--eval',help="file vcf o tsv da valutare")
	parser.add_argument('-c','--comp',help="file vcf o tsv da comparare per la valutazione del file in entrata")
	parser.add_argument('-o','--out',help="path di out in cui stampare le statistiche di valutazione")
	
	global opts
	opts = parser.parse_args()
	ineval = open(opts.eval,'r').readlines()
	
	csnv=dict()
	cindel=dict()
	esnv=dict()
	eindel=dict()
	
	if ineval[0].startswith("##fileformat=VCF"):
		start = [ineval.index(x) for x in ineval if x.startswith("#CHROM")][0]
		for line in ineval[start+1:]:
			#print line
			chr,pos,id,ref,alt,qual,filter,info,format,sample=line.rstrip().split('\t')
			id_var='\t'.join([chr,pos,ref,alt])
			#print format
			gt=sample.split(':')[format.split(':').index('GT')]
			if len(ref)==1 and len(alt)==1:
				esnv[id_var]=gt
			else:
				eindel[id_var]=gt
	else:
		header=ineval[0]
		for line in ineval[1:]:
			chr,pos,id,ref,alt=line.rstrip().split('\t')[:5]
			id_var='\t'.join([chr,pos,ref,alt])
			gt=line.rstrip().split('\t')[header.split('\t').index('GT')]
			if len(ref)==1 and len(alt)==1:
				esnv[id_var]=gt
			else:
				eindel[id_var]=gt

	incomp = open(opts.comp,'r').readlines()
	
	if incomp[0].startswith("##fileformat=VCF"):
		start = [incomp.index(x) for x in incomp if x.startswith("#CHROM")][0]
		for line in incomp[start+1:]:
			chr,pos,id,ref,alt,qual,filter,info,format,sample=line.rstrip().split('\t')
			id_var='\t'.join([chr,pos,ref,alt])
			gt=sample.split(':')[format.split(':').index('GT')]
			if len(ref)==1 and len(alt)==1:
				csnv[id_var]=gt
			else:
				cindel[id_var]=gt
	else:
		header=incomp[0]
		for line in incomp[1:]:
			chr,pos,id,ref,alt=line.rstrip().split('\t')[:5]
			id_var='\t'.join([chr,pos,ref,alt])
			gt=line.rstrip().split('\t')[header.split('\t').index('GT')]
			if len(ref)==1 and len(alt)==1:
				csnv[id_var]=gt
			else:
				cindel[id_var]=gt

	etot=dict(esnv, **eindel)
	ctot=dict(csnv, **cindel)
	stats_file=open(opts.out,'a+')
	print_header(stats_file)
	tp,fp,tn,fn,enum,cnum,gt_match,csnv=evaluate(esnv,csnv)
	if cnum > 0:
		f_disc_rate,tp_rate,precision,detect_rate,missing_rate,gt_match_rate,acc,scoref1=calc_stats(tp,fp,tn,fn,enum,cnum,gt_match)
		print_stats('SNV',f_disc_rate,tp_rate,precision,detect_rate,missing_rate,gt_match_rate,acc,scoref1,stats_file,None,enum,cnum)
	
	tp,fp,tn,fn,enum,cnum,gt_match,cindel=evaluate(eindel,cindel)
	if cnum > 0:
		f_disc_rate,tp_rate,precision,detect_rate,missing_rate,gt_match_rate,acc,scoref1=calc_stats(tp,fp,tn,fn,enum,cnum,gt_match)
		print_stats('INDEL',f_disc_rate,tp_rate,precision,detect_rate,missing_rate,gt_match_rate,acc,scoref1,stats_file,None,enum,cnum)

	tp,fp,tn,fn,enum,cnum,gt_match,ctot=evaluate(etot,ctot)
	f_disc_rate,tp_rate,precision,detect_rate,missing_rate,gt_match_rate,acc,scoref1=calc_stats(tp,fp,tn,fn,enum,cnum,gt_match)
	print_stats('TOTAL',f_disc_rate,tp_rate,precision,detect_rate,missing_rate,gt_match_rate,acc,scoref1,stats_file,ctot,enum,cnum)
