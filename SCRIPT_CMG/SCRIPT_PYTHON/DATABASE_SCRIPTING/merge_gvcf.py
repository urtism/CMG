import subprocess
import argparse
import os



def merge_gvcfs(gatk,reference,target,merged,gvcf,tomerge): 

	if not os.path.isfile(tomerge + '.idx'):
		args = [gatk,'--java-options', '-Xmx32g','IndexFeatureFile',"-F",tomerge]
		success = subprocess.call(args)

	if not os.path.isfile(gvcf+'.idx'):
		args = [gatk,'--java-options', '-Xmx32g','IndexFeatureFile',"-F",gvcf]
		success = subprocess.call(args)


	args = [gatk,'--java-options', '-Xmx32g','CombineGVCFs',"-R",reference,'-L',target,'-V',gvcf,'-V',tomerge,'-O',merged]
	success = subprocess.call(args)
	if not success:
		return merged


def main():
	parser = argparse.ArgumentParser('Merge gvcf using GATK CombineGVCFs')
	parser.add_argument('-V', '--gvcf_list', help="gvcflist containing gvcfs to merge")
	parser.add_argument('-g', '--group', help="number of gvcf to merge as subgroup",default=1)
	parser.add_argument('-R', '--reference', help="reference fasta file")
	parser.add_argument('-L', '--interval', help="interval file")
	parser.add_argument('-O','--out',help="out file")

	global opts 
	opts = parser.parse_args()

	gvcf_list = open(opts.gvcf_list,'r').readlines()
	gatk='/home/jarvis/NGS_TOOLS/GATK/v4.0.6.0/gatk'

	j=0
	while len(gvcf_list)!= 1:
		#print  "LEEEENG", len(gvcf_list)
		j+=1
		gvcf_sublist = []

		for i in range(0, len(gvcf_list),int(opts.group)):
			svcfsubgr = gvcf_list[i:i+int(opts.group)]

			

			subg_merged = '.'.join((opts.out).split('.')[:-1] + [str(j),str(i/10),'g.vcf'])
			gvcf_sublist += [subg_merged]

			
			for gvcf in svcfsubgr:
				merged = subg_merged

				if svcfsubgr.index(gvcf) == 0:
					gvcf=gvcf.rstrip()
					tomerge = gvcf
					print tomerge
				else:
					gvcf=gvcf.rstrip()
					merged = merge_gvcfs(gatk,opts.reference,opts.interval,merged,gvcf,tomerge)
					tomerge= '.'.join(merged.split('.')[:-2]) +'.tomerge.g.vcf'
					status = subprocess.call("mv " + merged + ' ' + tomerge, shell=True)
					status = subprocess.call("mv " + merged + '.idx ' + tomerge+'.idx', shell=True)

			status = subprocess.call("mv " + tomerge + ' ' + subg_merged, shell=True)		
			status = subprocess.call("mv " + tomerge + '.idx ' + subg_merged+'.idx', shell=True)

			#print len(gvcf_list[i:i+int(opts.group)])

		gvcf_list = gvcf_sublist

	status = subprocess.call("mv " + subg_merged + ' ' + opts.out, shell=True)		
	status = subprocess.call("mv " + subg_merged + '.idx ' + opts.out+'.idx', shell=True)


	exit()

	# for gvcf in gvcf_list:
	# 	merged = opts.out
	# 	print merged

	# 	if gvcf_list.index(gvcf) == 0:
	# 		gvcf=gvcf.rstrip()
	# 		tomerge = gvcf
	# 		print tomerge
	# 	else:
	# 		gvcf=gvcf.rstrip()
	# 		merged = merge_gvcfs(gatk,opts.reference,opts.interval,merged,gvcf,tomerge)
	# 		tomerge= '.'.join(merged.split('.')[:-2]) +'.tomerge.g.vcf'
	# 		status = subprocess.call("mv " + merged + ' ' + tomerge, shell=True)
	# 		status = subprocess.call("mv " + merged + '.idx ' + tomerge+'.idx', shell=True)

	# status = subprocess.call("mv " + tomerge + ' ' + merged, shell=True)		
	# status = subprocess.call("mv " + tomerge + '.idx ' + merged+'.idx', shell=True)
main()