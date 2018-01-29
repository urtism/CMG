import argparse

def main():

	parser = argparse.ArgumentParser('Parse merged VCF to calculate AC.')
	parser.add_argument('-m', '--merged', help="vcf merged from GATK,Freebayes and Varscan2")
	parser.add_argument('-o', '--out',help="output file.")

	global opts 
	opts = parser.parse_args()

	vcf_in = open(opts.merged,'r')
	vcf_out = open(opts.out,'w')

	for line in vcf_in:
		line=line.rstrip()
		if line.startswith('#'):
			vcf_out.write(line + '\n')
		else:
			chr,pos,id,ref,alt,qual,filter,info,format = line.split('\t')[:9]
			samples = line.split('\t')[9:]
			AC=0
			AN=0

			for sample in samples:
				gt_G = sample.split(':')[format.split(':').index('GT_G')]
				gt_F = sample.split(':')[format.split(':').index('GT_F')]
				gt_V = sample.split(':')[format.split(':').index('GT_V')]

				try:
					a_a_G = int(gt_G.split('/')[0]) + int(gt_G.split('/')[1]) 
				except:
					a_a_G = 0
				try:
					a_a_F = int(gt_F.split('/')[0]) + int(gt_F.split('/')[1]) 
				except:
					a_a_F = 0
				try:
					a_a_V = int(gt_V.split('/')[0]) + int(gt_V.split('/')[1]) 
				except:
					a_a_V = 0

				if gt_G!= './.'	or gt_F!= './.'	or gt_V!= './.':
					AN+=2			

				AC += max(a_a_G,a_a_F,a_a_V)
			info = info + ';AC=' + str(AC) +';AN=' +str(AN)

			vcf_out.write('\t'.join( [chr,pos,id,ref,alt,qual,filter,info,format]+samples)+'\n')

main()
