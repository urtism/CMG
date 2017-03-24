import argparse
import re
import string
import numpy
import pylab
import plotly
import matplotlib.pyplot as plt
import plotly.graph_objs
from collections import Counter

tab = open('/home/minime/Scrivania/SOMIC/20170217_01_T_Cardio_Ann_Cosmic-clnv.tsv','r')

var=[]
STATUS_Vardict=[]
DP_TUM=[]
DP_NORM=[]
AF_NORM=[]
AF_TUM=[]
Delta_mediana=[]
Delta_perc_median=[]
MBQ_TUM=[]
MBQ_NORM=[]
Consequence=[]
Existing_variation=[]
CLNDBN=[]
CLNORIGIN=[]
COSMIC_OCCURENCE=[]


Class = {}


Class['splice_donor_variant&intron_variant'] = 'Intron' 
Class['3_prime_UTR_variant'] = 'Intron'
Class['splice_acceptor_variant'] = 'Splice' 
Class['missense_variant&splice_region_variant'] = 'Splice & coding'
Class['inframe_deletion'] = 'Inframe Indel' 
Class['stop_gained&splice_region_variant'] = 'Nonsense' 
Class['downstream_gene_variant'] = 'Intron'
Class['5_prime_UTR_variant'] = 'Intron'
Class['synonymous_variant'] = 'Synonymous' 
Class['frameshift_variant&splice_region_variant'] = 'Frameshift' 
Class['protein_altering_variant'] = 'Protein altering'
Class['splice_donor_variant'] = 'Splice'
Class['inframe_insertion'] = 'Inframe Indel'
Class['splice_donor_variant&coding_sequence_variant&intron_variant'] = 'Splice & coding'
Class['stop_retained_variant'] = 'Stop retained'
Class['splice_acceptor_variant&coding_sequence_variant'] = 'Splice & coding'
Class['stop_gained'] = 'Nonsense' 
Class['upstream_gene_variant'] = 'Intron'
Class['splice_region_variant&synonymous_variant'] = 'Splice & coding'
Class['start_lost'] = 'Nonsense'
Class['frameshift_variant'] = 'Frameshift'
Class['splice_region_variant&intron_variant'] = 'Intron'
Class['missense_variant'] = 'Missense'
Class['intron_variant'] = 'Intron'

for line in tab:
	line = line.rstrip()
	line = line.split('\t')

	if line[0].startswith('CHROM'):
		header=line

	else:
		if line[header.index('MBQ_TUM')] != '.' and float(line[header.index('MBQ_TUM')]) >= 30.0 :
			if line[header.index('AF_TUM')] != '.' and float(line[header.index('AF_TUM')]) >= 0.05:
				if line[header.index('Delta_perc_median')] != '.' and float(line[header.index('Delta_perc_median')]) == 100:

					var  +=  [ str(line[header.index('CHROM')]) + '-' + str(line[header.index('POS')]) + '-' + str(line[header.index('REF')]) + '-' + str(line[header.index('ALT')]) ]
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
					Existing_variation += [line[header.index('Existing_variation')]]
					CLNDBN += [line[header.index('CLNDBN')]]
					CLNORIGIN += [line[header.index('CLNORIGIN')]]
					COSMIC_OCCURENCE += [line[header.index('COSMIC_OCCURENCE')]]

		else:
			continue

#Con questo comando posso contare le ricorrenze nella lista delle (ad esempio) consequence. Con count conto le ripetizioni nella lista
#mentre assegno alla consequence la chiave e al valore il count.

Count_Consequence = {i:Consequence.count(i) for i in Consequence}

print Count_Consequence



Chiave = []
Valore = []

for key in Count_Consequence.keys():
	Chiave += [key]
	Valore += [Count_Consequence.get(key)]

print Valore
print Chiave

# L'output sara' di questo tipo:
# {'splice_donor_variant&intron_variant': 1, '3_prime_UTR_variant': 220, 'splice_acceptor_variant': 31}


			# pylab.figure(1, figsize=(6,6))

			# patch, texts = pylab.pie(Valore, autopct='%1.1f%%', startangle=90)

			# pylab.legend(patch, Chiave, loc="best")

			# pylab.tight_layout()

			# pylab.show()


fig1, ax1 = plt.subplots()
ax1.pie(Valore, labels=Chiave, autopct='%1.1f%%', startangle=300, explode=(0.0,0.05,0.1,0.2,0.4,0.1,0.3))
ax1.axis('equal')
#plt.legend(Chiave, loc="best")
#plt.title('Mutazioni nel paz', bbox={'facecolor':'0.9', 'pad':5})
plt.savefig('/home/minime/Scrivania/SOMIC/20170217_01_T_Cardio_Mut-graph.png', bbox_inches='tight')