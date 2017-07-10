bed='/home/jarvis/NGS_ANALYSIS/TARGET/ctDNA_2_113416_AmpliconsExport.bed'

for line in open(bed,'r'):
	line=line.rstrip()
	chr=line.split('\t')[0]
	start=line.split('\t')[1]
	stop=line.split('\t')[2]

	start=int(start)-150
	stop=int(stop)+150

	print '\t'.join([chr,str(start),str(stop)])
