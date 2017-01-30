import argparse
import os


def main():

	parser = argparse.ArgumentParser('add chr to phylop and cons db. OUTPUT TO STDOUT')
	parser.add_argument('-p', '--path', help="PHYLOP AND CONS PATH")

	global opts 
	opts = parser.parse_args()


	
	for filename in os.listdir(opts.path):
		#print filename
		infile=open(opts.path+'/'+filename,'r')
		for line in infile:
			line = line.rstrip()
			if line.startswith('#') or line.startswith('track') :
				continue

			elif line.startswith('variableStep'):
				chr=(line.split(' ')[1]).split('=')[1]
			else:
				print '\t'.join([chr]+[line])

		infile.close()



main ()