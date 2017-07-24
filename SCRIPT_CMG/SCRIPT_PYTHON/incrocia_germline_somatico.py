import argparse

if __name__ == '__main__':

	parser = argparse.ArgumentParser('Confronta il file  delle varianti germinali con il file delle varianti somatiche ')
	parser.add_argument('-g','--germline',help="file vcf o tsv di varianti germinali")
	parser.add_argument('-s','--somatic',help="file vcf o tsv di varianti da tessuto somatico")
	parser.add_argument('-o','--out',help="path di out in cui stampare le statistiche di valutazione")
	
	global opts
	opts = parser.parse_args()

	somatic=open(opts.somatic,'r')
	germline=open(opts.germline,'r')
	out=open(opts.out,'w')
	germl=dict()
	somat=dict()

	perif=dict()
	cirr=dict()
	tum=dict()

	for line in germline:
		line=line.rstrip()
		if line.startswith('CHROM'):
			header_germ=line.split('\t')
		else:
			chr,pos,id,ref,alt=line.split('\t')[:5]
			gt=','.join([line.split('\t')[header_germ.index('GT_GATK')],line.split('\t')[header_germ.index('GT_Varscan')],line.split('\t')[header_germ.index('GT_Freebayes')]])

			germl['\t'.join([chr,pos,ref,alt])]=[line.split('\t'),gt]
	germline.close()

	for line in somatic:
		line=line.rstrip()
		if line.startswith('CHROM'):
			header_som=line.split('\t')
		else:
			chr,pos,id,ref,alt=line.split('\t')[:5]
			gt_tm=','.join([line.split('\t')[header_som.index('GT_t_Mutect')],line.split('\t')[header_som.index('GT_t_Varscan')],line.split('\t')[header_som.index('GT_t_Vardict')]])
			gt_cirr=','.join([line.split('\t')[header_som.index('GT_n_Mutect')],line.split('\t')[header_som.index('GT_n_Varscan')],line.split('\t')[header_som.index('GT_n_Vardict')]])
			somat['\t'.join([chr,pos,ref,alt])]=[line.split('\t'),gt_tm,gt_cirr]
	somatic.close()

	header_som[header_som.index('GT_n_Varscan')]='GT_cirr_Varscan'
	header_som[header_som.index('GT_n_Vardict')]='GT_cirr_Vardict'
	header_som[header_som.index('GT_n_Mutect')]='GT_cirr_Mutect'
	header_som[header_som.index('AF_n_Mutect')]='AF_cirr_Mutect'
	header_som[header_som.index('AF_n_Vardict')]='AF_cirr_Vardict'
	header_som[header_som.index('AF_n_Varscan')]='AF_cirr_Varscan'
	header_som[header_som.index('AF_norm_media')]='AF_cirr_media'
	header_som[header_som.index('DP_norm_media')]='DP_cirr_media'
	header_som[header_som.index('AO_norm_media')]='AO_cirr_media'
	header_som[header_som.index('RO_norm_media')]='RO_cirr_media'
	header_som[header_som.index('MBQN_media')]='MBQC_media'

	header=header_som
	header[5:7]=['TESSUTO','GT_PERIF']
	out.write('\t'.join(header) +'\n')

	for var in somat.keys():
		print var
		variante=somat[var][0]
		if var in germl.keys():
			variante[5:7]=['PERIFERICO',somat[var][1]]
			
		else:
			chiamata_som=[somat[var][0][header_som.index('SomaticVarscan')],somat[var][0][header_som.index('SomaticVardict')],somat[var][0][header_som.index('FILTER_Mutect')]]
			if '1' in chiamata_som or 'PASS' in chiamata_som:
				variante[5:7]=['TUMORE','0/0']
			else:
				variante[5:7]=['CIRR+TUMORE','0/0']
		out.write('\t'.join(variante) +'\n')



