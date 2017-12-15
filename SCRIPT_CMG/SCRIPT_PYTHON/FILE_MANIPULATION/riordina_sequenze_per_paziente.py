import argparse 
import os
#import SeqConverter
from shutil import copyfile

def main():
	
	parser = argparse.ArgumentParser('Genera per ogni paziente una cartella che contiene le sequenze in sanger')
	parser.add_argument('-a', '--file_a', help="lista di pazienti e di sequenze sanger")
	parser.add_argument('--illumina_list', help="lista pazienti illumina",default=None)
	parser.add_argument('-o', '--out', help="dir di output")
	

	global opts 
	opts = parser.parse_args()


	illumina_paz=dict()

	if opts.illumina_list:
		illumina_list=open(opts.illumina_list,'r')
		for line in illumina_list:
			line=line.rstrip()
			if line.startswith('SC_ID'):
				header=line.split('\t')
			else:
				sc_id=line.split('\t')[header.index('SC_ID')]
				nome=line.split('\t')[header.index('NOME_PAZ')]
				id=line.split('\t')[header.index('ID')]
				num_run=line.split('\t')[header.index('NUM_RUN')]
				data_run=line.split('\t')[header.index('DATA_RUN')]
				mean_dp=line.split('\t')[header.index('MEAN_DEPTH')]
				illumina_paz[nome]=[sc_id,id]


	sanger=open(opts.file_a,'r')
	sanger_paz=dict()
	sanger_paz_phd=dict()

	for line in sanger:
		line=line.rstrip()
		if line.startswith('ID_SC'):
			header=line.split('\t')
		else:

			nome=line.split('\t')[header.index('Sample')]
			cod=line.split('\t')[header.index('Index')]
			path_split=(line.split('\t')[header.index('Path')]).split('/')
			path_split[0]='/home/jarvis'
			path='/'.join(path_split)
			nome_dir=path_split[-1]
			ab1_path=''
			phd_path=''
			seq_path=''
			#print '/'.join(path_split[:-1])
			try:			
				for direc in os.listdir('/'.join(path_split[:-1])):
					if nome_dir in direc:
						path='/'.join(['/'.join(path_split[:-1]),direc])
						break
			except:
				print line,'/'.join(path_split[:-1])

			for file in os.listdir(path):
				if cod in file:
					seq_path = path+'/'+file
					if nome in sanger_paz.keys():
						sanger_paz[nome]+=[seq_path]						
					else:						
						sanger_paz[nome]=[seq_path]
					
			if seq_path == '':
				print line

	for paz in sanger_paz.keys():
		i=0
		if paz in illumina_paz.keys():
			id=illumina_paz[paz][1]
			sc_id=illumina_paz[paz][0]
			if sc_id != '.' and sc_id != '':
				new_dir=opts.out+'/'+sc_id
			elif id != '.' and id != '':
				new_dir=opts.out+'/'+id
			else:
				new_dir=opts.out+'/UNKNOWN'

			try:
				os.mkdir(new_dir)
			except:
				pass
		else:
			new_dir=opts.out+'/'+'_'.join(paz.split(' '))
			try:
				os.mkdir(new_dir)
			except:
				pass
		


		#print sanger_paz[paz]

		for path in sanger_paz[paz]:
			i+=1			
			#new_path= new_dir+'/'+str(i)+'.ab1'
			if 'phd' in path:
				new_path= new_dir+'/'+str(i)+'.phd'
			elif 'ab' in path:
				new_path= new_dir+'/'+str(i)+'.ab1'
			#print path,os.path.isfile(path)

			if os.path.isfile(path):
				copyfile(path, new_path)
				print path,new_path
			else:
				print paz,path + ' NON TROVATO'


main()