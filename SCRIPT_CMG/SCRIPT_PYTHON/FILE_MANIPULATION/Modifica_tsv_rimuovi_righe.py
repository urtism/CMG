import argparse
import re
import string


def main():

	parser = argparse.ArgumentParser('Con questo script vado ad eliminare le righe che sono annotate male a causa del sito multiallelico. Usato su esperimento somic')
	
	parser.add_argument('-i','--input',help="file tab delimited da corregere")
	parser.add_argument('-o','--outfile',help="file di output tab delimited corretto")

	global opts
	
	opts = parser.parse_args()
	tab = open(opts.input,'r')
	out = open(opts.outfile,'w')

	set1=[]

	for line in tab:
		line = line.rstrip()
		line = line.split('\t')

		if len(line) == 85:
			set1.append(line)
				#out.write('\t'.join(line) + '\n')

		else:
			continue


	for elem in set1:
		out.write('\t'.join(elem) + '\n')

	# 	set1.append(line)
	# 	set2.append(line)


	# for elem in set1:
	# 	if len()

	# for elem in set1:

	# 	id_fam = str(elem[0]+';')

	# 	for val in set2:
	# 		if elem[0] == val[0]:
	# 			continue
	# 		if elem[1].split('|')[-1] in val[1].split('|')[-1] or elem[1].split('|')[-1] in val[2].split('|')[-1] or elem[2].split('|')[-1] in val[1].split('|')[-1] or elem[2].split('|')[-1] in val[2].split('|')[-1]:
	# 			id_fam = str(id_fam) + str(val[0]) + ';'
	# 	link.append(id_fam)

	# newlink = []

	# for value in link:
	# 	value = value.split(';')[:-1]
	# 	value.sort()
	# 	newlink.append(';'.join(value))
	
	# newlink = list(set(newlink))
	
	# for var in newlink:
	# 	out.write(str(var) + '\n')
	
main()