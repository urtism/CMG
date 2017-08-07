import argparse
import statistics as st
import numpy as np
import os


class SNV():
	variant='SNV'
	num_called_variant=[]
	num_real_variants=[]
	detection_rate=[]
	False_discover_rate=[]
	Sensitivity=[]
	precision=[]
	accuracy=[]
	score_f1=[]
	missing_rate=[]
	gt_match_rate=[]

class INDEL():
	variant='INDEL'
	num_called_variant=[]
	num_real_variants=[]
	detection_rate=[]
	False_discover_rate=[]
	Sensitivity=[]
	precision=[]
	accuracy=[]
	score_f1=[]
	missing_rate=[]
	gt_match_rate=[]

class TOTAL():
	variant='TOTAL'
	num_called_variant=[]
	num_real_variants=[]
	detection_rate=[]
	False_discover_rate=[]
	Sensitivity=[]
	precision=[]
	accuracy=[]
	score_f1=[]
	missing_rate=[]
	gt_match_rate=[]
		
if __name__ == '__main__':

	parser = argparse.ArgumentParser('Prende in ingresso il path contenente i file delle statistiche di validazione e calcola le statistiche medie su tutti i file')
	parser.add_argument('-i','--path',help="path dove sono presenti i file.stats da valutare")
	parser.add_argument('-o','--out',help="file di out con le statistiche generali")

	global opts
	opts = parser.parse_args()

	snv=SNV()
	indel=INDEL()
	total=TOTAL()
	out=open(opts.out,'w')

	for file in os.listdir(opts.path):
		sample=file.split('.')[0]
		
		for line in open(opts.path+"/"+ file,'r') :

			line=line.rstrip()
			#print line
			if line.startswith('VARIANTI'):

				header=line.split('\t')
				out.write('SAMPLE'+ '\t' +line +'\n')

			elif line.startswith('SNV'):
				out.write(sample+ '\t' +line +'\n')
				
				snv.variant='SNV'
				snv.num_called_variant+=[int(line.split('\t')[header.index('NUM CALLED VARIANT')])]
				snv.num_real_variants+=[int(line.split('\t')[header.index('NUM REAL VARIANTS')])]
				snv.detection_rate+=[float(line.split('\t')[header.index('DETECTION RATE')])]
				snv.False_discover_rate+=[float(line.split('\t')[header.index('FALSE DISCOVER RATE')])]
				snv.Sensitivity+=[float(line.split('\t')[header.index('SENSITIVITY (TP RATE)')])]
				snv.precision+=[float(line.split('\t')[header.index('PRECISION')])]
				snv.accuracy+=[float(line.split('\t')[header.index('ACCURACY')])]
				snv.score_f1+=[float(line.split('\t')[header.index('SCORE F1')])]
				snv.missing_rate+=[float(line.split('\t')[header.index('MISSING RATE')])]
				snv.gt_match_rate+=[float(line.split('\t')[header.index('GT MATCH RATE')])]
							
			elif line.startswith('INDEL'):
				#print 'indel'
				out.write(sample+ '\t' +line +'\n')
				indel.variant='INDEL'
				indel.num_called_variant+=[int(line.split('\t')[header.index('NUM CALLED VARIANT')])]
				indel.num_real_variants+=[int(line.split('\t')[header.index('NUM REAL VARIANTS')])]
				indel.detection_rate+=[float(line.split('\t')[header.index('DETECTION RATE')])]
				indel.False_discover_rate+=[float(line.split('\t')[header.index('FALSE DISCOVER RATE')])]
				indel.Sensitivity+=[float(line.split('\t')[header.index('SENSITIVITY (TP RATE)')])]
				indel.precision+=[float(line.split('\t')[header.index('PRECISION')])]
				indel.accuracy+=[float(line.split('\t')[header.index('ACCURACY')])]
				indel.score_f1+=[float(line.split('\t')[header.index('SCORE F1')])]
				indel.missing_rate+=[float(line.split('\t')[header.index('MISSING RATE')])]
				indel.gt_match_rate+=[float(line.split('\t')[header.index('GT MATCH RATE')])]

			elif line.startswith('TOTAL'):
				#print 'total'
				out.write(sample+ '\t' +line +'\n')
				total.variant='TOTAL'
				total.num_called_variant+=[int(line.split('\t')[header.index('NUM CALLED VARIANT')])]
				total.num_real_variants+=[int(line.split('\t')[header.index('NUM REAL VARIANTS')])]
				total.detection_rate+=[float(line.split('\t')[header.index('DETECTION RATE')])]
				total.False_discover_rate+=[float(line.split('\t')[header.index('FALSE DISCOVER RATE')])]
				total.Sensitivity+=[float(line.split('\t')[header.index('SENSITIVITY (TP RATE)')])]
				total.precision+=[float(line.split('\t')[header.index('PRECISION')])]
				total.accuracy+=[float(line.split('\t')[header.index('ACCURACY')])]
				total.score_f1+=[float(line.split('\t')[header.index('SCORE F1')])]
				total.missing_rate+=[float(line.split('\t')[header.index('MISSING RATE')])]
				total.gt_match_rate+=[float(line.split('\t')[header.index('GT MATCH RATE')])]

			
	
	#out.write('\t'.join(["VARIANTI","NUM CALLED VARIANT","NUM REAL VARIANTS","MEAN DETECTION RATE","STDV DETECTION RATE","MEAN FALSE DISCOVER RATE","STDV FALSE DISCOVER RATE","MEAN SENSITIVITY","STDV SENSITIVITY","MEAN PRECISION","STDV PRECISION","MEAN ACCURACY","STDV ACCURACY","MEAN SCORE F1","STDV SCORE F1","MEAN MISSING RATE","STDV MISSING RATE","MEAN GT MATCH RATE","STDV GT MATCH RATE"])+'\n')
	
	for vars in [snv,indel,total]:
		
		NUM_CALLED_VARIANT = str(int(round(np.sum(vars.num_called_variant),4)))
		NUM_REAL_VARIANTS = str(int(round(np.sum(vars.num_real_variants),4)))

		MEAN_DETECTION_RATE = str(round(st.mean(vars.detection_rate),4))
		STDV_DETECTION_RATE = str(round(st.stdev(vars.detection_rate),4))

		MEAN_FALSE_DISCOVER_RATE = str(round(st.mean(vars.False_discover_rate),4))
		STDV_FALSE_DISCOVER_RATE = str(round(st.stdev(vars.False_discover_rate),4))

		MEAN_SENSITIVITY = str(round(st.mean(vars.Sensitivity),4))
		STDV_SENSITIVITY = str(round(st.stdev(vars.Sensitivity),4))

		MEAN_PRECISION = str(round(st.mean(vars.precision),4))
		STDV_PRECISION = str(round(st.stdev(vars.precision),4))

		MEAN_ACCURACY = str(round(st.mean(vars.accuracy),4))
		STDV_ACCURACY = str(round(st.stdev(vars.accuracy),4))

		MEAN_SCORE_F1 = str(round(st.mean(vars.score_f1),4))
		STDV_SCORE_F1 = str(round(st.stdev(vars.score_f1),4))

		MEAN_MISSING_RATE = str(round(st.mean(vars.missing_rate),4))
		STDV_MISSING_RATE = str(round(st.stdev(vars.missing_rate),4))

		MEAN_GT_MATCH_RATE = str(round(st.mean(vars.gt_match_rate),4))
		STDV_GT_MATCH_RATE = str(round(st.stdev(vars.gt_match_rate),4))

		out.write('\t'.join(['GROUP',vars.variant,NUM_CALLED_VARIANT,NUM_REAL_VARIANTS,MEAN_DETECTION_RATE,MEAN_FALSE_DISCOVER_RATE,MEAN_SENSITIVITY,MEAN_PRECISION,MEAN_ACCURACY,MEAN_SCORE_F1,MEAN_MISSING_RATE,MEAN_GT_MATCH_RATE])+'\n')
	
		#out.write('\t'.join([vars.variant,NUM_CALLED_VARIANT,NUM_REAL_VARIANTS,MEAN_DETECTION_RATE,STDV_DETECTION_RATE,MEAN_FALSE_DISCOVER_RATE,STDV_FALSE_DISCOVER_RATE,MEAN_SENSITIVITY,STDV_SENSITIVITY,MEAN_PRECISION,STDV_PRECISION,MEAN_ACCURACY,STDV_ACCURACY,MEAN_SCORE_F1,STDV_SCORE_F1,MEAN_MISSING_RATE,STDV_MISSING_RATE,MEAN_GT_MATCH_RATE,STDV_GT_MATCH_RATE])+'\n')
	