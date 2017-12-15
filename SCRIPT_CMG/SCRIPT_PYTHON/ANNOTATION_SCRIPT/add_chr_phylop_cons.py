import argparse
import os


def main():

	parser = argparse.ArgumentParser('add chr to phylop and cons db. OUTPUT TO STDOUT')
	parser.add_argument('-p', '--path', help="PHYLOP AND CONS PATH")

	global opts 
	opts = parser.parse_args()


	for path in  os.listdir(opts.path):
		newpath=opts.path+'/'+path
		for filename in  os.listdir(newpath):
			#print filename
			infile=open(newpath+'/'+filename,'r')
			outfile=open(newpath+'/'+filename.split('.')[0]+'.chr.tsv','w')
			print newpath+'/'+filename

			for line in infile:
				line = line.rstrip()
				if line.startswith('#') or line.startswith('track') :
					continue

				elif line.startswith('variableStep'):
					chr=(line.split(' ')[1]).split('=')[1]
					outfile.write('\t'.join(['CHROM','POS','SCORE']) +'\n')
				else:
					outfile.write('\t'.join([chr]+[line]) +'\n')

			infile.close()
			outfile.close()



main ()