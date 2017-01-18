import argparse
import os

def extract(lista):
	result=open(opts.out,'w')

	# for i in (1,2,3,4,5,6,7,8,9):
	# 	gatk[i]=open(opts.run + '/' + opts.data +'_0' + str(i) +'_'+ opts.tipo+'_GATK_Conferme.ods','r')
	# 	freebayes[i]=open(opts.run + '/' + opts.data +'_0' + str(i) +'_' + opts.tipo+'_FREEBAYES_Conferme.ods','r')
	# 	varscan[i]=open(opts.run + '/' + opts.data +'_0'+ str(i) +'_' + opts.tipo+'_varsc_Conferme.ods','r')
	# 	other[i]=open(opts.run + '/' + opts.data +'_0' + str(i) +'_'+ opts.tipo+'_GATK_Other_Transcripts_Conferme.ods','r')
	# for a in (10,11,12):
	# 	gatk[i]=open(opts.run + '/' + opts.data +'_' + str(a) +'_'+ opts.tipo+'_GATK_Conferme.ods','r')
	# 	freebayes[i]=open(opts.run + '/' + opts.data +'_' + str(a) +'_' + opts.tipo+'_FREEBAYES_Conferme.ods','r')
	# 	varscan[i]=open(opts.run + '/' + opts.data +'_'+ str(a) +'_' + opts.tipo+'_GATK_Conferme.ods','r')
	# 	other[i]=open(opts.run + '/' + opts.data +'_' + str(a) +'_'+ opts.tipo+'_GATK_Other_Transcripts_Conferme.ods','r')

	for line in lista:
		trovata = 0
		line= line.rstrip()
		num_paz= line.split('\t')[0]
		pos= line.split('\t')[1]

		if int(num_paz)<10:
			num_paz= '0'+num_paz

		os.system('mv '+ opts.run + '/' + opts.data +'_' + num_paz +'_'+ opts.tipo+'_GATK_Conferme.ods '+ opts.run + '/' + opts.data +'_' + num_paz +'_'+ opts.tipo+'_GATK_Conferme.tsv')	
		os.system('mv '+ opts.run + '/' + opts.data +'_' + num_paz +'_'+ opts.tipo+'_FREEBAYES_Conferme.ods ' + opts.run + '/' + opts.data +'_' + num_paz +'_'+ opts.tipo+'_FREEBAYES_Conferme.tsv')	
		os.system('mv '+ opts.run + '/' + opts.data +'_' + num_paz +'_'+ opts.tipo+'_VarScan_Conferme.ods ' +opts.run + '/' + opts.data +'_' + num_paz +'_'+ opts.tipo+'_VarScan_Conferme.tsv')	
		os.system('mv '+ opts.run + '/' + opts.data +'_' + num_paz +'_'+ opts.tipo+'_GATK_Other_Transcripts_Conferme.ods ' +opts.run + '/' + opts.data +'_' + num_paz +'_'+ opts.tipo+'_GATK_Other_Transcripts_Conferme.tsv')	
		
		gatk=open(opts.run + '/' + opts.data +'_' + num_paz +'_'+ opts.tipo+'_GATK_Conferme.tsv','r')
		freebayes=open(opts.run + '/' + opts.data +'_' + num_paz +'_' + opts.tipo+'_FREEBAYES_Conferme.tsv','r')
		varscan=open(opts.run + '/' + opts.data +'_'+ num_paz +'_' + opts.tipo+'_GATK_Conferme.tsv','r')
		other=open(opts.run + '/' + opts.data +'_' + num_paz +'_'+ opts.tipo+'_GATK_Other_Transcripts_Conferme.tsv','r')

		for var in gatk:
			var=var.rstrip()
			res='\t'.join(['GATK',var.split('\t')[0],var.split('\t')[1],var.split('\t')[2],var.split('\t')[3],var.split('\t')[4],var.split('\t')[24],var.split('\t')[25],var.split('\t')[26],var.split('\t')[27],var.split('\t')[28]])
			if var.split('\t')[1] == pos and trovata ==0 :
				result.write(res+'\n')
				trovata=1
		for var in freebayes:
			var=var.rstrip()
			res='\t'.join(['freebayes',var.split('\t')[0],var.split('\t')[1],var.split('\t')[2],var.split('\t')[3],var.split('\t')[4],var.split('\t')[24],var.split('\t')[25],var.split('\t')[26],var.split('\t')[27],var.split('\t')[28]])
			if var.split('\t')[1] == pos and trovata ==0 :
				result.write(res+'\n')
				trovata=1
		for var in varscan:
			var=var.rstrip()
			res='\t'.join(['varscan',var.split('\t')[0],var.split('\t')[1],var.split('\t')[2],var.split('\t')[3],var.split('\t')[4],var.split('\t')[24],var.split('\t')[25],var.split('\t')[26],var.split('\t')[27],var.split('\t')[28]])
			if var.split('\t')[1] == pos and trovata ==0 :
				result.write(res+'\n')
				trovata=1
		for var in other:
			var=var.rstrip()
			res='\t'.join(['GATK_OTHER',var.split('\t')[0],var.split('\t')[1],var.split('\t')[2],var.split('\t')[3],var.split('\t')[4],var.split('\t')[24],var.split('\t')[25],var.split('\t')[26],var.split('\t')[27],var.split('\t')[28]])
			if var.split('\t')[1] == pos and trovata ==0 :
				result.write(res+'\n')
				trovata=1
		if trovata==0:
				result.write('NON TROVATA'+'\n')

		gatk.close()
		freebayes.close()
		varscan.close()
		other.close()
		


def main():
	parser = argparse.ArgumentParser('Adds class to the opts.dataset of variants.  Output is to stdout.')
	
	parser.add_argument('-f','--file',help="opts.dataset di varianti tab delimited")
	parser.add_argument('-r','--run',help="true somatic vcf",default= None)
	parser.add_argument('-d','--data',help="true somatic vcf",default= None)
	parser.add_argument('-t','--tipo',help="true somatic vcf",default= None)
	parser.add_argument('-o','--out',help="true somatic vcf",default= None)

	
	global opts
	opts = parser.parse_args()
	
	extract(open(opts.file,'r'))
	
main()