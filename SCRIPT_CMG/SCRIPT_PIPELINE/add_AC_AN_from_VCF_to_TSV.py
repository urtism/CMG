import argparse

parser = argparse.ArgumentParser('Incrocia un file tsv di varianti con il file vcf mergiato con i campi AC e AN e le aggiunge al tsv')
parser.add_argument('-i', '--tsv', help="tsv di varianti a cui aggiungere AC e AN")
parser.add_argument('-v', '--vcf', help="file vcf mergiato con i campi AC e AN" )
parser.add_argument('-o', '--out', help="file tsv con aggiunti i campi AC e AN" )

opts = parser.parse_args()

tsv_in = open(opts.tsv,'r')



varianti=dict()

for line in tsv_in:
	line=line.rstrip()
	if line.startswith('CHROM'):
		header=line.split('\t')
	else:
		line_split=line.split('\t')
		chr,pos,id,ref,alt = line_split[:5]
		idvar='\t'.join([chr,pos,ref,alt])
		varianti[idvar]=line_split

tsv_in.close()

vcf = open(opts.vcf,'r')
for line in vcf:
	if line.startswith('#CHROM'):
		HEADER_VCF=line.split('\t')
	elif line.startswith('##'):
		continue
	else:
		line_split=line.split('\t')
		info_split=line_split[7].split(';')
		chr,pos,id,ref,alt = line_split[:5]
		idvar='\t'.join([chr,pos,ref,alt])
		if idvar in varianti.keys():
			temp=varianti[idvar]
			AC=info_split[-2].split('=')[1]
			AN=info_split[-1].split('=')[1]
			temp[12:12]=[AC,AN]
			varianti[idvar]=temp

vcf.close()

tsv_out = open(opts.out,'w')

header[12:12]=['AC','AN']

tsv_out.write('\t'.join(header)+'\n')
for var in varianti.keys():
	tsv_out.write('\t'.join(varianti[var])+'\n')
