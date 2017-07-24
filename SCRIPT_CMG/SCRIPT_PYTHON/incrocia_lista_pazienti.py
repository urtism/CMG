import argparse


def main():
	
	parser = argparse.ArgumentParser('Incrocia due liste di pazienti illumina, controlla quali pazienti nel file_b sono presenti nel file_a. Output is to --out.')
	parser.add_argument('-a', '--file_a', help="lista di pazienti a")
	parser.add_argument('-b', '--file_b', help="lista di pazienti b")
	parser.add_argument('-o', '--out', help="file name in output")
	

	global opts 
	opts = parser.parse_args()

	file_a = open(opts.file_a,'r')
	file_b = open(opts.file_b,'r')
	out_a = open(opts.out+'.list','w')
	out_b = open(opts.out+'_in_comune.list','w')

	pazienti_a=dict()

	for line in file_a:
		line=line.rstrip()
		if line.startswith('NOME_PAZ'):
			continue
		else:
			line_split=line.split('\t')
			nome=line_split[0]
			pazienti_a[nome]=line

	for line in file_b:
		line=line.rstrip()
		if line.startswith('NOME_PAZ'):
			continue
		else:
			line_split=line.split('\t')
			nome=line_split[0]
			if nome in pazienti_a.keys():
				out_b.write(line + '\n')
			else:
				pazienti_a[nome]=line


	for p in pazienti_a.keys():
		out_a.write(pazienti_a[p] +'\n')

main()