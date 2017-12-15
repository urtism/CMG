from Bio import SeqIO
import argparse



def convert(file_a,tipo_a,file_b,tipo_b):
	SeqIO.convert(file_a, tipo_a,file_b, tipo_b)


def main():
	
	parser = argparse.ArgumentParser('converte sequenze in diversi tipi di formato')
	parser.add_argument('-a', '--file_a', help="file da convertire")
	parser.add_argument('-t_a', '--tipo_a', help="dipo di file da convertire")
	parser.add_argument('-b', '--file_b', help="file convertito")
	parser.add_argument('-t_b', '--tipo_b', help="dipo di file nel quale convertire")
	

	global opts 
	opts = parser.parse_args()

	convert(opts.file_a,opts.tipo_a,opts.file_b,opts.tipo_b)

main() 



