import argparse
import numpy as np
import random as rm


def calcola_freq(gt,perc,err):
	if gt==2:
		perc=perc*2
	delta=err*perc
	perc_rand=rm.randrange(int((perc-delta)*10000),int((perc+delta)*10000))
	if perc_rand/10000.0 > 1.0:
		perc_rand=10000
	return str(float(perc_rand)/10000.0)



def main():
	parser = argparse.ArgumentParser('Genera il file delle varianti da simulare con BamSurgeon')
	
	parser.add_argument('-i','--file1',help="file dei pazienti con i genotipi degli rs simulati")
	parser.add_argument('-s','--file2',help="lista degli rs con le coordinate cromosomiche")
	parser.add_argument('-c','--coppie',help="numero di coppie donor recipient da simulare")
	parser.add_argument('-r','--rumore',help="frequenza di rumore da simulare",default=None)
	parser.add_argument('-o','--out',help="path di output")

	global opts
	opts = parser.parse_args()
	rs_list=dict()

	soglie_freq_rec=[0.25,0.5,0.75]
	soglie_freq_don=[0.008,0.02,0.05,0.1]
	err=0.1

	for line in open(opts.file2,'r'):
		line=line.rstrip()
		if line.startswith('Locus'):
			continue
		else:
			rs=line.split('\t')[0]
			chr=line.split('\t')[1]
			pos=line.split('\t')[2]
			alt=line.split('\t')[4]

			rs_list[rs]=['chr'+chr,pos,pos,alt]

	list_file=(open(opts.file1,'r')).readlines()
	#list_file2=(lambda list_file: list_file[list_file.index(el)]=el.rstrip() for el in list_file)
	#print list_file2
	#L = map(uova, prosciutto)
	for el in list_file:
		list_file[list_file.index(el)]=el.rstrip()


	list_rs=(list_file[0].rstrip()).split('\t')[1:]
	#print list_rs
	list_paz=list_file[1:-3]

	coppie=[]
	for i in opts.coppie:
		i_don=rm.randrange(1,len(list_paz))
		i_rec=rm.randrange(1,len(list_paz))
		while i_don == i_rec:
			i_rec=rm.randrange(1,len(list_paz))
		coppie+=[[list_paz[i_don],list_paz[i_rec]]]

	for rs in list_rs:
		print rs,rs_list[rs]

	for c in coppie:
		ind_coppia=str(coppie.index(c)+1)
		don_mut=[]
		rec_mut=[]
		don=c[0].split('\t')
		rec=c[1].split('\t')
		print len(rec)
		nome_don=don[0]
		nome_rec=rec[0]

		for fq in soglie_freq_rec:
			out_rec=open(opts.out+'/C'+ind_coppia+'.rec.'+str(fq)+'.list','w')	
			for rs in list_rs:
				#print list_rs.index(rs) +1
				if rec[list_rs.index(rs) +1] =='0':
					if opts.rumore:
						out_rec.write('\t'.join(rs_list[rs][:3])+'\t')
						out_rec.write(str(calcola_freq(1,float(opts.rumore),err))+'\t')
						out_rec.write('\t'.join(rs_list[rs][-1])+'\n')
					else:
						continue
				elif rec[list_rs.index(rs) +1] =='1':
					out_rec.write('\t'.join(rs_list[rs][:3])+'\t')
					out_rec.write(str(calcola_freq(1,fq,err))+'\t')
					out_rec.write('\t'.join(rs_list[rs][-1])+'\n')

				elif rec[list_rs.index(rs) +1] =='2':
					out_rec.write('\t'.join(rs_list[rs][:3])+'\t')
					out_rec.write(str(calcola_freq(2,fq,err))+'\t')
					out_rec.write('\t'.join(rs_list[rs][-1])+'\n')
			out_rec.close()

		for fq in soglie_freq_don:
			out_don=open(opts.out+'/C'+ind_coppia+'.don.'+str(fq)+'.list','w')
			for rs in list_rs:
				if rs == 'rs1028528':
					print rs+'\t'+'\t'.join(rs_list[rs]),don[list_rs.index(rs) +1]
				if don[list_rs.index(rs) +1] =='0':
					if opts.rumore:
						print rs+'\t'+'\t'.join(rs_list[rs][:3])+'\t'+str(calcola_freq(1,float(opts.rumore),err))
						out_don.write('\t'.join(rs_list[rs][:3])+'\t')
						out_don.write(str(calcola_freq(1,float(opts.rumore),err))+'\t')
						out_don.write('\t'.join(rs_list[rs][-1])+'\n')
					else:
						continue
				elif don[list_rs.index(rs) +1] =='1':
					print rs+'\t'+'\t'.join(rs_list[rs][:3])+'\t'+str(calcola_freq(1,fq,err)+'\t')
					out_don.write('\t'.join(rs_list[rs][:3])+'\t')
					out_don.write(str(calcola_freq(1,fq,err))+'\t')
					out_don.write('\t'.join(rs_list[rs][-1])+'\n')
				elif don[list_rs.index(rs) +1] =='2':
					print rs+'\t'+'\t'.join(rs_list[rs][:3])+'\t'+str(calcola_freq(2,fq,err))
					out_don.write('\t'.join(rs_list[rs][:3])+'\t')
					out_don.write(str(calcola_freq(2,fq,err))+'\t')
					out_don.write('\t'.join(rs_list[rs][-1])+'\n')
			out_don.close()

main()