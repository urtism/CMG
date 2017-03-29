import argparse
import re
import string
import numpy
import pylab
import plotly
import matplotlib.pyplot as plt
import plotly.graph_objs
from collections import Counter


def estrai_stat(tab,out):

	id_var=[]
	def_var=[]
	var=[]
	STATUS_Vardict=[]
	DP_TUM=[]
	DP_NORM=[]
	AF_NORM=[]
	AF_TUM=[]
	AF_tum_media=[]
	DP_tum_media=[]
	Delta_mediana=[]
	Delta_perc_median=[]
	MBQ_TUM=[]
	MBQT_media=[]
	MBQ_NORM=[]
	VARIANT_CLASS=[]
	Consequence=[]
	SYMBOL=[]
	Amino_acids=[]
	Codons=[]
	HGVSc=[]
	HGVSp=[]
	Existing_variation=[]
	GMAF=[]
	CLNDBN=[]
	CLNORIGIN=[]
	COSMIC_OCCURENCE=[]

	Class = {}

	Class['splice_acceptor_variant&coding_sequence_variant&intron_variant'] = 'Splice & coding'
	Class['frameshift_variant&stop_lost'] = 'Frameshift'
	Class['splice_region_variant&5_prime_UTR_variant'] = 'Intron'
	Class['splice_donor_variant&intron_variant'] = 'Intron'
	Class['3_prime_UTR_variant'] = 'Intron'
	Class['splice_acceptor_variant'] = 'Splice'
	Class['missense_variant&splice_region_variant'] = 'Splice & coding'
	Class['inframe_deletion'] = 'Inframe Indel'
	Class['stop_gained&splice_region_variant'] = 'Nonsense'
	Class['stop_gained&frameshift_variant'] = 'Nonsense'
	Class['downstream_gene_variant'] = 'Intron'
	Class['5_prime_UTR_variant'] = 'Intron'
	Class['synonymous_variant'] = 'Synonymous'
	Class['frameshift_variant&splice_region_variant'] = 'Frameshift'
	Class['protein_altering_variant'] = 'Protein altering'
	Class['splice_donor_variant'] = 'Splice'
	Class['splice_acceptor_variant&intron_variant'] = 'Splice'
	Class['inframe_insertion'] = 'Inframe Indel'
	Class['splice_donor_variant&coding_sequence_variant&intron_variant'] = 'Splice & coding'
	Class['stop_retained_variant'] = 'Stop retained'
	Class['splice_acceptor_variant&coding_sequence_variant'] = 'Splice & coding'
	Class['stop_gained'] = 'Nonsense' 
	Class['upstream_gene_variant'] = 'Intron'
	Class['splice_region_variant&synonymous_variant'] = 'Splice & coding'
	Class['stop_lost'] = 'Nonsense'
	Class['start_lost'] = 'Nonsense'
	Class['frameshift_variant'] = 'Frameshift'
	Class['frameshift_variant&start_lost'] = 'Frameshift'
	Class['splice_region_variant&intron_variant'] = 'Intron'
	Class['missense_variant'] = 'Missense'
	Class['intron_variant'] = 'Intron'

	for line in tab:
		line = line.rstrip()
		line = line.split('\t')

		if line[0].startswith('CHROM'):
			header=line

		else:
			#if 'DP_TUM' in header and 'MBQ_TUM' in header and 'Delta_perc_median' in header and 'AF_TUM' in header:
			if float(line[header.index('DP_TUM')]) >= 40.0 and line[header.index('MBQ_TUM')] != '.' and float(line[header.index('MBQ_TUM')]) >= 30.0:
				if line[header.index('AF_TUM')] != '.' and float(line[header.index('AF_TUM')]) >= 0.05:
					if line[header.index('Delta_perc_median')] != '.' and float(line[header.index('Delta_perc_median')]) == 100:

						id_var=[]
						var=[]
						STATUS_Vardict=[]
						DP_TUM=[]
						DP_tum_media=[]
						DP_NORM=[]
						AF_NORM=[]
						AF_TUM=[]
						Delta_mediana=[]
						Delta_perc_median=[]
						MBQ_TUM=[]
						MBQT_media=[]
						MBQ_NORM=[]
						VARIANT_CLASS=[]
						Consequence=[]
						SYMBOL=[]
						Amino_acids=[]
						Codons=[]
						HGVSc=[]
						HGVSp=[]
						Existing_variation=[]
						GMAF=[]
						CLNDBN=[]
						CLNORIGIN=[]
						COSMIC_OCCURENCE=[]

						var += [ str(line[header.index('CHROM')]) + '-' + str(line[header.index('POS')]) + '-' + str(line[header.index('REF')]) + '-' + str(line[header.index('ALT')]) ]
						STATUS_Vardict += [line[header.index('STATUS_Vardict')]]
						DP_TUM += [line[header.index('DP_TUM')]]
						DP_NORM += [line[header.index('DP_NORM')]]
						AF_NORM += [line[header.index('AF_NORM')]]
						AF_TUM += [line[header.index('AF_TUM')]]
						Delta_mediana += [line[header.index('Delta_mediana')]]
						Delta_perc_median += [line[header.index('Delta_perc_median')]]
						MBQ_TUM += [line[header.index('MBQ_TUM')]]
						MBQ_NORM += [line[header.index('MBQ_NORM')]]
						Consequence += [Class.get(line[header.index('Consequence')])]
						SYMBOL += [line[header.index('SYMBOL')]]
						Amino_acids=[line[header.index('Amino_acids')]]
						Codons=[line[header.index('Codons')]]
						HGVSc=[line[header.index('HGVSc')]]
						HGVSp=[line[header.index('HGVSp')]]
						Existing_variation += [line[header.index('Existing_variation')]]
						GMAF=[line[header.index('GMAF')]]
						CLNDBN += [line[header.index('CLNDBN')]]
						CLNORIGIN += [line[header.index('CLNORIGIN')]]
						COSMIC_OCCURENCE += [line[header.index('COSMIC_OCCURENCE')]]
						id_var = var + DP_TUM + AF_TUM + MBQ_TUM + Consequence + SYMBOL + Amino_acids + Codons + HGVSc + HGVSp + Existing_variation + GMAF + CLNDBN + CLNORIGIN + COSMIC_OCCURENCE
						def_var = def_var + [id_var]

			else:
				continue

	#Con questo comando posso contare le ricorrenze nella lista delle (ad esempio) consequence. Con count conto le ripetizioni nella lista
	#mentre assegno alla consequence la chiave e al valore il count.

	Count_Consequence = {i:Consequence.count(i) for i in Consequence}

	#print Count_Consequence

	Chiave = []
	Valore = []

	for key in Count_Consequence.keys():
		Chiave += [key]
		Valore += [Count_Consequence.get(key)]

	#out.write('Varianti somatiche filtrate con: DP_TUM >= 40 & MBQ_TUM >= 30 & AF_TUM >= 0.03 & Delta_perc_median >= 100' + '\n\n')

	#for i in range(len(Chiave)):
	#	new_list=[]
	#	new_list = [Chiave[i]] + [str(Valore[i])]
	#	print new_list
	#	out.write('\t'.join(new_list) + '\n')

	#out.write('\n\nConteggio varianti' + '\n\n')

	out.write('ID_VAR' + '\t' + 'DP_TUM' + '\t' + 'AF_TUM' + '\t' + 'MBQ_TUM' + '\t' + 'Consequence' + '\t' + 'SYMBOL' + '\t' + 'Amino_acids' + '\t' + 'Codons' + '\t' + 'HGVSc' + '\t' + 'HGVSp' + '\t' + 'Existing_variation' + '\t' + 'GMAF' + '\t' + 'CLNDBN' + '\t' + 'CLNORIGIN' + '\t' + 'COSMIC_OCCURENCE' + '\n')

	for varianti in def_var:
		out.write('\t'.join(varianti) + '\n')


# L'output sara' di questo tipo:
# {'splice_donor_variant&intron_variant': 1, '3_prime_UTR_variant': 220, 'splice_acceptor_variant': 31}



						#Definisco l'explode dell'immagine (quanto una fetta deve essere staccata dalle altre)
						# exp = range(len(Chiave))
						# exp = [float(spe)/4 for spe in exp]
						# print exp


						# fig1, ax1 = plt.subplots()
						# ax1.pie(Valore, labels=Chiave, autopct='%1.1f%%', startangle=300, explode=exp)
						# #explode=(0.0,0.05,0.1,0.2,0.4,0.1,0.3)
						# ax1.axis('equal')
						# plt.legend(Chiave, loc="best")
						# plt.title('Mutazioni nel paz', bbox={'facecolor':'0.9', 'pad':5})
						# plt.savefig('/home/minime/Scrivania/SOMIC/20170221_06_T_Cardio_Mut-graph.png', bbox_inches='tight')




def main():

	parser = argparse.ArgumentParser('Estraggo le statistiche dal file di input')
	parser.add_argument('-i','--input',help="file tab delimited")
	parser.add_argument('-o','--outfile',help="file di output tab delimited")

	global opts

	opts = parser.parse_args()
	tab = open(opts.input,'r')
	out = open(opts.outfile,'w')

	print opts.input

	estrai_stat(tab,out)

	tab.close()
	out.close()

main()