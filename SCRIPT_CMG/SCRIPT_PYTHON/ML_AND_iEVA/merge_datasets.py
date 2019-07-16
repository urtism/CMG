import os
import tarfile
import argparse
import pandas as pd
import hashlib
import seaborn
import re
import numpy as np
import matplotlib.pyplot as plt
from six.moves import urllib
from joblib import dump, load

from pandas.api.types import is_object_dtype
from pandas.plotting import scatter_matrix


def add_vc(dataset,tag):

	for f in dataset.columns.values.tolist():
		if f != 'ID':
			newf = tag + '_' + f
			dataset = dataset.rename(index=str, columns={f: newf})
	return dataset


if __name__ == '__main__':

	parser = argparse.ArgumentParser('Questo script contiene le funzioni base per la manipolazione del dataset.\n')
	parser.add_argument('-g', '--gatk', default=None, help="Dataset da analizzare")
	parser.add_argument('-f', '--freebayes', default=None, help="Split del dataset in TESTset e TRAINset, indicare la test_size in decimali (es. 0.3 se si vuole avere un 70/30)")
	parser.add_argument('-v', '--varscan', default=None, help="Label utilizzato per stratificare il dataset splittato")
	parser.add_argument('-p', '--platypus', default=None, help="Classificatore allenato da utilizzare per la classificazione di un nuovo dataset")
	parser.add_argument('-s', '--samtools', default=None, help="train set")
	parser.add_argument('-n', '--snver', default=None, help="train set")
	parser.add_argument('-l', '--scalpel', default=None, help="train set")
	parser.add_argument('-O', '--outpath',default=None, help="Path di output")

	global opts
	opts = parser.parse_args()

	MERGE = pd.read_csv(opts.gatk,sep='\t')
	#freebayes = pd.read_csv(opts.freebayes,sep='\t')
	varscan = pd.read_csv(opts.varscan,sep='\t')
	#platypus = pd.read_csv(opts.platypus,sep='\t') 
	#samtools = pd.read_csv(opts.samtools,sep='\t')
	#snver = pd.read_csv(opts.snver,sep='\t')
	#scalpel = pd.read_csv(opts.scalpel,sep='\t')
	
	#gatk = add_vc(gatk,'GATK')
	#freebayes = add_vc(freebayes,'FB')
	varscan = add_vc(varscan,'VS')
	#platypus = add_vc(platypus,'PT')
	#samtools = add_vc(samtools,'ST')
	#snver = add_vc(snver,'SV')
	#scalpel = add_vc(scalpel,'SC')

	#MERGE= pd.merge(gatk, freebayes, on='ID', how='outer')
	MERGE= pd.merge(MERGE, varscan, on='ID', how='outer')
	#MERGE= pd.merge(MERGE, platypus, on='ID', how='outer')
	#MERGE= pd.merge(MERGE, samtools, on='ID', how='outer')
	#MERGE= pd.merge(MERGE, snver, on='ID', how='outer')
	#MERGE= pd.merge(MERGE, scalpel, on='ID', how='outer')

	MERGE.set_index(['ID'], inplace =True)


	MERGE.to_csv(path_or_buf=opts.outpath,sep='\t')





