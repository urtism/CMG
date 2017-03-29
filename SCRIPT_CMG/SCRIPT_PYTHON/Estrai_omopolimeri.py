import argparse
import re
import string
















def main():

	parser = argparse.ArgumentParser('Con questo script vado ad eliminare le righe che sono annotate male a causa del sito multiallelico. Usato su esperimento somic')
	
	parser.add_argument('-i','--input',help="vcf contenente CHROM e POS")
	parser.add_argument('-o','--outfile',help="file vcf annotato in INFO con una tag del tipo Repeats=Tandem/Homopolymer/null")

