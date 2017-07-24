import random as rm
import numpy as np

num_paz=10000
lista_paz=dict()

rs_out=open('/home/jarvis/Scrivania/rs_random.list','w')
header=[]
for paz in range(num_paz):
	paz=str(paz)
	#print paz
	lista_paz[paz]=[]
	rs_file=open('/home/jarvis/Scrivania/snp_trap.adjust.list','r')
	for rs in rs_file:
		if rs.startswith('Locus'):
			continue
		else:
			id_rs=rs.split('\t')[0]
			if id_rs not in header:
				header+=[id_rs]
			pAA=rs.split('\t')[6]
			pAB=rs.split('\t')[7]
			pBB=rs.split('\t')[8]

			randrange=float(rm.randrange(1, num_paz))
			if randrange <= float(pAA)*num_paz:
				 #print randrange,str(float(pAA)*num_paz)
				lista_paz[paz]+=['0']
				
			elif randrange > float(pAA)*num_paz and randrange <= (float(pAA)*num_paz+float(pAB)*num_paz):
				lista_paz[paz]+=['1']
			else:
				lista_paz[paz]+=['2']
	rs_file.close()
	#print lista_paz[paz]			
rs_out.write('\t'.join(['Locus']+header)+'\n')
for paz in lista_paz.keys():
	rs_out.write('\t'.join([paz]+lista_paz[paz])+'\n')

gt_00=np.zeros(94)
gt_01=np.zeros(94)
gt_11=np.zeros(94)
count=0
"""
for paz in lista_paz.keys():
	for paz2 in lista_paz.keys():
		if paz==paz2:
			continue
		else:
			count+=1
			i=0
			for gt in lista_paz[paz]:
				print gt
				if gt == '0' and lista_paz[paz2][i] == '0':
					gt_00[i]+=1
				elif gt == '0' and lista_paz[paz2][i] == '1':
					gt_01[i]+=1
				elif gt == '0' and lista_paz[paz2][i] == '2':
					gt_11[i]+=1
				i+=1


rs_out.write('0/0 0/0\t'+'\t'.join((str(el/count) for el in gt_00))+'\n')
rs_out.write('0/0 0/1\t'+'\t'.join((str(el/count) for el in gt_01))+'\n')
rs_out.write('0/0 0/2\t'+'\t'.join((str(el/count) for el in gt_11))+'\n')

"""
#rs_out.write('\t'.join((str(el) for el in ['0/0']+gt_00))+'\n')
#rs_out.write('\t'.join((str(el) for el in ['0/1']+gt_01))+'\n')
#rs_out.write('\t'.join((str(el) for el in ['1/1']+gt_02))+'\n')

print count