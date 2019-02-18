import argparse
import subprocess
import json
import os
import textwrap
import random as rm
import datetime
import pysam
import random as r

def prRed(prt): print("\033[91m {}\033[00m" .format(prt))
def prGreen(prt): print("\033[92m {}\033[00m" .format(prt))

def makedirs(dirs):
    for d in dirs:
        if not os.path.exists(d):
			os.makedirs(d)

def vars_from_db(db,num_snv,num_indel,out):

	database=open(db,'r').readlines()
	vcf=open(out,'w')
	start = [database.index(x) for x in database if x.startswith("#CHROM")][0]
	varianti=database[start:]
	i=0
	simulate=[['chr1','0']]
	while i < int(num_snv):
		rand=r.randrange(len(varianti))
		var=varianti[rand].rstrip().split('\t')
		del(varianti[rand])
		alt=var[4].split(',')[0]
		chr=var[0]
		pos=var[1]
		if len(var[3])==1 and len(alt)==1:
			for v in simulate:
				if chr == v[0] and int(pos) < int(v[1]) + 150 and chr == v[0] and int(pos) > int(v[1]) - 150:
					#print v,var
					pass
				else:
					i+=1
					vcf.write('\t'.join(var[:4]+[alt,'.','.','.','.','.'])+'\n')
					#print var
					simulate+=[[chr,pos]]
					break
		else:	
			continue 
	varianti=database[start:]
	i=0
	while i < int(num_indel):
		rand=r.randrange(len(varianti))
		var=varianti[rand].rstrip().split('\t')
		del(varianti[rand])
		alt=var[4].split(',')[0]
		if len(var[3])>1 or len(alt)>1:
			for v in simulate:
				if chr == v[0] and int(pos) < int(v[1]) + 150 and int(pos) > int(v[1]) - 150:
					#print v,var
					pass
				else:
					i+=1
					vcf.write('\t'.join(var[:4]+[alt,'.','.','.','.','.'])+'\n')
					#print var
					simulate+=[[chr,pos]]
					break
		else:
			continue
	vcf.close()
	return out

def check_simul_vars(log_bs_dir,bam,vars,outpath):

	for dir in os.listdir(log_bs_dir):
		if bam.split('/')[-1] in dir:
			log_dir=dir
		else:
			continue
		bam_log=outpath+'/'+bam.split('/')[-1]+'.log'
		varianti= open(vars,'r')
		log= open(bam_log,'w')
		for line in varianti:
			line=line.rstrip()
			id_var='_'.join(line.split('\t')[:3])
			for file in os.listdir(log_bs_dir+'/'+log_dir):
				if id_var in file:
					try:
						line_to_check=(open(log_bs_dir+'/'+log_dir+'/'+file,'r').readlines())[-1]
					except:
						var=line.split('\t')[:3]
						log.write('\t'.join(var+['.','NOT SIMULATED']) +'\n')
						print ':'.join(line.split('\t')[:2]),'--> variant not simulated'
						continue
					if line_to_check.startswith('indel'):
						vaf=str(round(float(line_to_check.split('\t')[-2]),3))
						type=line.split('\t')[4]
						var=line.split('\t')[:3]
						log.write('\t'.join(var+[vaf,type,'SIMULATED']) +'\n')
					elif line_to_check.startswith('snv'):
						vaf=str(round(float(line_to_check.split('\t')[-2]),3))
						var=line.split('\t')[:3]
						log.write('\t'.join(var+[vaf,'SIMULATED']) +'\n')
					else:
						var=line.split('\t')[:3]
						log.write('\t'.join(var+['.','NOT SIMULATED']) +'\n')
						print ':'.join(line.split('\t')[:2]),'--> variant not simulated'
		varianti.close()
		log.close()

def print_header(vcf):
	header=True
	for line in vcf:
		if line.startswith('#'):
			header=False
			break
		else:
			header=True
	if header:
		vcf.write(textwrap.dedent("""\
		##fileformat=VCFv4.1
		##phasing=none
		##INDIVIDUAL=TRUTH
		##SAMPLE=<ID=TRUTH,Individual="TRUTH",Description="bamsurgeon spike-in">
		##INFO=<ID=Somatic,Number=0,Type=Flag,Description="Somatic mutation in primary">
		##INFO=<ID=Germline,Number=0,Type=Flag,Description="Germline mutation in primary">
		##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSPIKEIN""")+'\n')

def print_vcf(analisi,log_bs_dir,outpath,reference):
	vcf=open(outpath+'/'+analisi+'.vcf','a+')
	print_header(vcf)

	for pathname in os.listdir(log_bs_dir):
		if os.path.isdir(log_bs_dir+'/'+pathname):
			for filename in os.listdir(log_bs_dir+'/'+pathname):
				#print filename
				if filename.endswith('.log'):
					with open(log_bs_dir +'/'+pathname + '/' + filename, 'r') as infile:
						for line in infile:
							if line.startswith('snv'):
								#chrom, pos, mut = line.strip().split()
								c = line.strip().split()
								chrom = c[1].split(':')[0]
								pos = c[3]
								mut = c[4]
								dpr = c[6]
								vaf = str(round(float(line.split('\t')[-2]),3))
								if float(vaf) >= 0.75:
									gt='1/1'
								elif float(vaf) > 0.00 and float(vaf) < 0.75:
									gt='0/1'

								ref,alt = mut.split('-->')
								try:
									ref=ref.upper()
								except:
									pass
								try:
									alt=alt.upper()
								except:
									pass
								vcf.write( "\t".join((chrom,pos,'.',ref,alt,'.','PASS', analisi+';VAF=' + vaf,'GT',gt))+'\n')

							elif line.startswith('indel'):
								vaf=str(round(float(line.split('\t')[-2]),3))
								if float(vaf) >= 0.75:
									gt='1/1'
								elif float(vaf) > 0.00 and float(vaf) < 0.75:
									gt='0/1'
								fa = pysam.Fastafile(reference)
								indelinfo = line.strip().split()[1].split(':')
								if indelinfo[0] == 'DEL':
									chrom, start, end = indelinfo[1:4]
									ref = fa.fetch(chrom, int(start)-1, int(end)).upper()
									alt = ref[0].upper()

								if indelinfo[0] == 'INS':
									chrom, start, seq = indelinfo[1:4]
									ref = fa.fetch(chrom, int(start)-1, int(start)).upper()
									alt = ref.upper() + seq.upper()
									
								assert ref != '' and alt != '' and start != ''

								vcf.write('\t'.join((chrom, start, '.', ref, alt, '.', 'PASS', analisi+';VAF=' + vaf , 'GT', gt))+'\n')
					subprocess.call("mv " + log_bs_dir+'/'+pathname+ '/'+filename + ' ' + log_bs_dir+'/'+pathname+ '/'+filename +'.checked', shell=True)

def Vcf_to_bamsurgeon(vars,min,max,err):
	print vars
	snp=open('/'.join(vars.split('.')[:-1])+'.snp','w')
	indel=open('/'.join(vars.split('.')[:-1])+'.indel','w')
	vcf=open(vars,'r')
	for line in vcf:
		line=(line.rstrip()).split('\t')
		if line[0].startswith('#'):
			continue
		else:
			chr=line[0]
			pos=line[1]
			id=line[2]
			ref=line[3]
			alt=line[4]
			format=line[-2]
			freq=line[-1]
			if freq == '.' or format != 'vaf':
				freq=Freq_calc(min,max,err)
			if len(ref)==1 and len(alt)==1:
				if alt =='.':
					snp.write('\t'.join([chr,pos,pos,freq]) + '\n')
				else:
					snp.write('\t'.join([chr,pos,pos,freq,alt]) + '\n')
			elif len(ref)>1 and len(alt)==1:
				del_pos=str(int(pos) + len(ref) - len(alt))
				indel.write('\t'.join([chr,pos,del_pos,freq,'DEL']) + '\n')
			elif len(ref)==1 and len(alt)>1:
				INS=alt[1:]
				indel.write('\t'.join([chr,pos,str(int(pos)+1),freq,'INS',INS]) + '\n')
	snp.close()
	indel.close()
	vcf.close()
	return '/'.join(vars.split('.')[:-1])+'.snp','/'.join(vars.split('.')[:-1])+'.indel'

def Freq_calc(min,max,err):
	freq=float(rm.randrange(int(float(min)*100),int(float(max)*100)))/100.0
	if err != None:
		delta=err*freq
		freq_rand=rm.randrange(int((freq-delta)*10000),int((freq+delta)*10000))
		if freq_rand/10000.0 > 1.0:
			freq_rand=10000
		freq=freq_rand	
	return str(freq)

def Set_bamsurgeon(conf):
	ms_args=[]
	ms_args+=['--picardjar',conf["picard"]["path"]]
	if conf["bamsurgeon"]["aligner"]!="":
		ms_args+=['--aligner',conf["bamsurgeon"]["aligner"]]
	else:
		ms_args+=['--aligner','mem']
	if conf["bamsurgeon"]["snvfrac"]!="":
		ms_args+=['--snvfrac',conf["bamsurgeon"]["snvfrac"]]
	else:
		pass
	if conf["bamsurgeon"]["mutfrac"]!="":
		ms_args+=['--mutfrac',conf["bamsurgeon"]["mutfrac"]]
	else:
		pass
	if conf["bamsurgeon"]["numsnvs"]!="":
		ms_args+=['--numsnvs',+conf["bamsurgeon"]["numsnvs"]]
	else:
		pass
	if conf["bamsurgeon"]["coverdiff"]!="":
		ms_args+=['--coverdiff',conf["bamsurgeon"]["coverdiff"]]
	else:
		pass
	if conf["bamsurgeon"]["haplosize"]!="":
		ms_args+=['--haplosize',conf["bamsurgeon"]["haplosize"]]
	else:
		pass
	if conf["bamsurgeon"]["procs"]!="":
		ms_args+=['--procs', conf["bamsurgeon"]["procs"]]
	else:
		pass
	if conf["bamsurgeon"]["mindepth"]!="":
		ms_args+=['--mindepth',conf["bamsurgeon"]["mindepth"]]
	else:
		pass
	if conf["bamsurgeon"]["maxdepth"]!="":
		ms_args+=['--maxdepth',conf["bamsurgeon"]["maxdepth"]]
	else:
		pass
	if conf["bamsurgeon"]["minmutreads"]!="":
		ms_args+=['--minmutreads',conf["bamsurgeon"]["minmutreads"]]
	else:
		pass
	if conf["bamsurgeon"]["avoidreads"]!="":
		ms_args+=['--avoidreads',conf["bamsurgeon"]["avoidreads"]]
	else:
		pass
	if conf["bamsurgeon"]["maxopen"]!="":
		ms_args+=['--maxopen',conf["bamsurgeon"]["maxopen"]]
	else:
		pass
	if conf["bamsurgeon"]["seed"]!="":
		ms_args+=['--seed',conf["bamsurgeon"]["seed"]]
	else:
		pass

	ms_args+=conf["bamsurgeon"]["flags"].split(',')
	path_addsnv=conf["bamsurgeon"]["addsnv"]
	path_addindel=conf["bamsurgeon"]["addindel"]

	return path_addsnv,path_addindel,ms_args

def Set_art(conf):
	art_args=[]
	if conf["art"]["qprof1"]!="":
		art_args+= ['--qprof1', conf["art"]["qprof1"]]
	else:
		pass
	if conf["art"]["qprof2"]!="":
		art_args+= ['--qprof2', conf["art"]["qprof2"]]
	else:
		pass
	if conf["art"]["rcount"]!="":
		art_args+= ['--rcount', conf["art"]["rcount"]]
	else:
		pass
	if conf["art"]["id"]!="":
		art_args+= ['--id', conf["art"]["id"]]
	else:
		pass
	if conf["art"]["insRate"]!="":
		art_args+= ['--insRate', conf["art"]["insRate"]]
	else:
		pass
	if conf["art"]["insRate2"]!="":
		art_args+= ['--insRate2', conf["art"]["insRate2"]]
	else:
		pass
	if conf["art"]["delRate"]!="":
		art_args+= ['--delRate', conf["art"]["delRate"]]
	else:
		pass
	if conf["art"]["delRate2"]!="":
		art_args+= ['--delRate2', conf["art"]["delRate2"]]
	else:
		pass
	if conf["art"]["maxIndel"]!="":
		art_args+= ['--maxIndel', conf["art"]["maxIndel"]]
	else:
		pass
	if conf["art"]["len"]!="":
		art_args+= ['--len', conf["art"]["len"]]
	else:
		art_args+= ['--len','150']
	if conf["art"]["mflen"]!="":
		art_args+= ['--mflen', conf["art"]["mflen"]]
	else:
		art_args+= ['--mflen','300']
	if conf["art"]["maskN"]!="":
		art_args+= ['--maskN', conf["art"]["maskN"]]
	else:
		pass
	if conf["art"]["minQ"]!="":
		art_args+= ['--minQ', conf["art"]["minQ"]]
	else:
		pass
	if conf["art"]["maxQ"]!="":
		art_args+= ['--maxQ', conf["art"]["maxQ"]]
	else:
		pass
	if conf["art"]["qShift"]!="":
		art_args+= ['--qShift', conf["art"]["qShift"]]
	else:
		pass
	if conf["art"]["qShift2"]!="":
		art_args+= ['--qShift2', conf["art"]["qShift2"]]
	else:
		pass
	if conf["art"]["rndSeed"]!="":
		art_args+= ['--rndSeed', conf["art"]["rndSeed"]]
	else:
		pass
	if conf["art"]["fcov"]!="":
		art_args+= ['--fcov', conf["art"]["fcov"]]
	else:
		art_args+= ['--fcov', '1000']
	if conf["art"]["sdev"]!="":
		art_args+= ['--sdev', conf["art"]["sdev"]]
	else:
		art_args+= ['--sdev', '10']
	if conf["art"]["seqSys"]!="":
		art_args+= ['--seqSys', conf["art"]["seqSys"]]
	else:
		art_args+= ['--seqSys', 'MSv3']
	art_args+=conf["art"]["flags"].split(',')
	path_art=conf["art"]["path"]
	if conf["art"]["prefix"]!="":
		art_args+= ['--out', opts.out_path+'/FASTQ/' +conf["art"]["prefix"]]
	else:
		art_args+= ['--out', opts.out_path+'/FASTQ/art']
	return path_art,art_args

def Set_bwa(conf):
	bwa_args=[]
	if conf["bwa"]["aligner"]!="":
		bwa_args+=['mem']
	else:
		bwa_args+=[conf["bwa"]["aligner"]]
	if conf["bwa"]["threads"]!="":
		bwa_args+=['-t 8']
	else:
		bwa_args+=['-t '+ conf["bwa"]["threads"]]

	if conf["bwa"]["output-prefix"]!="":
		bwa_args+=[opts.out_path+'/BAM/' + conf["bwa"]["output-prefix"]]
	else:
		bwa_args+=[opts.out_path+'/BAM/simul.bam']

	path_bwa=conf["bwa"]["path"]

	return path_bwa,bwa_args

def Set_pirs_args(conf):

	pirs_arguments=[]

	if conf["pirs"]["converage"]!="":
		pirs_arguments+=['--coverage='+conf["pirs"]["converage"]]
	else:
		pirs_arguments+=['--coverage=1000']
	if conf["pirs"]["read-len"]!="":
		pirs_arguments+=["--read-len="+conf["pirs"]["read-len"]]
	else:
		pirs_arguments+=["--read-len=150"]
	if conf["pirs"]["insert-len-mean"]!="":
		pirs_arguments+=["--insert-len-mean="+ conf["pirs"]["insert-len-mean"]]
	else:
		pass
	if conf["pirs"]["insert-len-sd"]!="":
		pirs_arguments+=["--insert-len-sd="+conf["pirs"]["insert-len-sd"]]
	else:
		pass
	if conf["pirs"]["base-calling-profile"]!="":
		pirs_arguments+=["--base-calling-profile="+conf["pirs"]["base-calling-profile"]]
	else:
		pass
	if conf["pirs"]["indel-error-profile"]!="":
		pirs_arguments+=["--indel-error-profile="+conf["pirs"]["indel-error-profile"]]
	else:
		pass
	if conf["pirs"]["error-rate"]!="":
		pirs_arguments+=["--error-rate="+conf["pirs"]["error-rate"]]
	else:
		pass
	if conf["pirs"]["substitution-error-algorithm"]!="":
		pirs_arguments+=["--substitution-error-algorithm="+conf["pirs"]["substitution-error-algorithm"]]
	else:
		pass
	if conf["pirs"]["quality-shift"]!="":
		pirs_arguments+=["--quality-shift="+conf["pirs"]["quality-shift"]]
	else:
		pass
	if conf["pirs"]["output-prefix"]!="":
		pirs_arguments+=["--output-prefix="+opts.out_path+'/FASTQ/' +conf["pirs"]["output-prefix"]]
	else:
		pirs_arguments+=["--output-prefix=pirs"]
	if conf["pirs"]["threads"]!="":
		pirs_arguments+=["--threads="+conf["pirs"]["threads"]]
	else:
		pirs_arguments+=["--threads=4"]
	
	pirs_arguments+=conf["pirs"]["flags"].split(',')
	
	pirs_arguments+=["--output-file-type=gzip"]
	path_pirs=conf["pirs"]["path"]
	return path_pirs,pirs_arguments

def Set_picard(conf):
	picard_args=[]
	if conf["picard"]["ram"]!="":
		picard_args+=[conf["picard"]["ram"]]
	else:
		picard_args+=[]

	path_picard=conf["picard"]["path"]
	return path_picard,picard_args

def Set_samtools(conf):
	return conf["samtools"]["path"]

###################################################################################################################################################################################################

def Bam_surgeon(path_bs,opts_bs,inbam,outbam,vars,log):
	print 'Simulating variants using BamSurgeon in: '+ inbam
	tmpdir='/'.join(outbam.split('/')[:-1]+['tmp'])
	makedirs([tmpdir])
	args = ['python',path_bs,'-v',vars,'-f', inbam,'-r', opts.ref,'-o',outbam,'--tmpdir',tmpdir] + opts_bs
	success = subprocess.call(args,stdout=log,stderr=log)
	try:
		os.rmdir(tmpdir)
	except:
		pass
	
	if not success:
		print 'Simulation: DONE ---> New bam:'+outbam
	else:
		prRed('Error in Variant spikein. Check log file.')
		exit(1)

def Simulate_fastq_pirs(path_pirs,opts_pirs,fasta,log):
	print 'Fastq simulation from: '+fasta
	if '--diploid' in opts_pirs:
		opts_pirs+=[fasta,fasta]
	else:
		opts_pirs+=[fasta]
	#print opts_pirs
	args = [path_pirs,'simulate'] + opts_pirs
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		print 'Fastq simulation: DONE'
	else:
		prRed('Error in Fastq simulation. Check log file.')
		exit(1)

def Simulate_fastq_art(path_art,opts_art,fasta,log):
	print 'Fastq simulation from: '+fasta
	
	args = [path_art,'-sam','-i',fasta] + opts_art
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		print 'Fastq simulation: DONE'
	else:
		prRed('Error in Fastq simulation. Check log file.')
		exit(1)


def Alignment_bwa(path_bwa,opts_bwa,fastq1,fastq2,log):
	
	sam=opts_bwa[-1]+'.sam'
	sam_file=open(sam,'w+')
	if fastq1==None:
		print 'Alignment using BWA:'+fastq2
		args = [path_bwa,opts_bwa[0],opts.ref,fastq2,opts_bwa[1]]
	elif fastq2==None:
		print 'Alignment using BWA:'+fastq1
		args = [path_bwa,opts_bwa[0],opts.ref,fastq1,opts_bwa[1]]
	else:
		print 'Alignment using BWA:'+fastq1 + ' '+ fastq2
		args = [path_bwa,opts_bwa[0],opts.ref,fastq1,fastq2,opts_bwa[1]]
	
	success = subprocess.call(args, stdout=sam_file,stderr=log)

	if not success:
		print 'Alignment: DONE ---> Sam: '+sam
		return sam
	else:
		prRed('Error in Alignment. Check log file.')
		exit(1)
	

def Index_bam(path_picard,opts_picard,bam,log):
	print 'Indexing bam'
	bai=bam+'.bai'
	args = ['java',opts_picard[0],'-jar',path_picard,'BuildBamIndex','I='+bam,'O='+bai,'VALIDATION_STRINGENCY=LENIENT']
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		pass
	else:
		prRed('Error in Indexing. Check log file.')
		exit(1)

def SamFormatConverter(path_picard,opts_picard,sam,log):

	print 'From Sam to Bam'
	bam='.'.join(sam.split('.')[:-1])+'.bam'
	args = ['java',opts_picard[0],'-jar',path_picard,'SamFormatConverter','I='+sam,'O='+bam]
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		return bam
	else:
		prRed('Error in Conversion. Check log file.')
	

def SortSam(path_picard,opts_picard,bam,log):

	print 'Sorting bam'
	sort='.'.join(bam.split('.')[:-1])+'.sort.bam'
	args = ['java',opts_picard[0],'-jar',path_picard,'SortSam','I='+bam,'O='+sort,'SORT_ORDER=coordinate']
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		return sort
	else:
		prRed('Error in Sorting. Check log file.')
		exit(1)


def getfasta(infasta,bed,log):
	print 'Getting Fasta from bed'

	outfasta=opts.out_path+ '/'+'.'.join((infasta.split('/')[-1]).split('.')[:-1])+'.bed.fasta'
	args = ['bedtools', 'getfasta', '-fo', outfasta, '-fi', infasta, '-bed' ,bed]
	success = subprocess.call(args,stdout=log,stderr=log)
	
	if not success:
		print 'New Fasta generated from reference and bed: '+ outfasta+'\n'
		return outfasta
	else:
		prRed('Error generating new fasta. Check log file.')
		exit(1)

def Index_fasta(fasta,st_path):
	args = [st_path, 'faidx',fasta]
	success = subprocess.call(args,stdout=log,stderr=log)
	
if __name__ == '__main__':

	parser = argparse.ArgumentParser('A Pipeline to simulate Germline and Somatic variants')
	parser.add_argument('--bam', help="Bam to simulate the variants",default=None)
	parser.add_argument('-fq1', '--fastq1', help="Fastq 1 to simulate the variants",default=None)
	parser.add_argument('-fq2', '--fastq2', help="Fastq 2 to simulate the variants",default=None)
	parser.add_argument('-fa', '--fasta', help="Fasta to generate fastq using pirs",default=None)
	parser.add_argument('-r', '--ref', help="Reference.fasta",default=None)
	parser.add_argument('-b', '--bed', help="Bed file to generate fastq files",default=None)
	parser.add_argument('--vars', help="Germline variant list to simulate in vcf format",default=None)
	parser.add_argument('--vars_som', help="Somatic variant list to simulate in vcf format",default=None)
	parser.add_argument('--dbsnp', help="dbsnp path for variant random extraction",default=None)
	parser.add_argument('--num_snv_dbsnp', help="number of snv to extract from db_snp",default=None)
	parser.add_argument('--num_indel_dbsnp', help="number of indel to extract from db_snp",default=None)
	parser.add_argument('--cosmic', help="cosmic path for variant random extraction",default=None)
	parser.add_argument('--num_snv_cosmic', help="number of snv to extract from cosmic",default=None)
	parser.add_argument('--num_indel_cosmic', help="number of indel to extract from cosmic",default=None)
	parser.add_argument('--min', help="Min threshold of simulating frequence for Germline variants",default='0.5')
	parser.add_argument('--max', help="Max threshold of simulating frequence for Germline variants",default='0.5')
	parser.add_argument('--min_som', help="Min threshold of simulating frequence for Somatic variants",default='0.05')
	parser.add_argument('--max_som', help="Max threshold of simulating frequence for Somatic variants",default='0.25')
	parser.add_argument('--err', help="Simulating frequence error. Variants will be simulated with freq between freq +- err*freq.",default=None)
	parser.add_argument('-a', '--analysis',choices=['Somatic','Germline'],help="Somatic or Germline",default='Germline')
	parser.add_argument('-o', '--out_path',help="output path")
	parser.add_argument('-c', '--cfg',help="configuration file")
	parser.add_argument('--amplicon',help="Amplicon design for fastq simulation",action='store_true')

	global opts
	opts = parser.parse_args()

	if opts.cfg != None:
		conf = json.loads((open(opts.cfg).read()).encode('utf8'))
	else:
		conf= json.loads((open(os.path.dirname(os.path.abspath(__file__))+'/config.json').read()).encode('utf8'))
 

	fastq_path=opts.out_path+'/FASTQ/'
	bam_path=opts.out_path+'/BAM/'
	log_bs_dir=opts.out_path+'/BS_log/'
	sim_log_dir=opts.out_path+'/Simulation_LOGS/'
	log_path=opts.out_path+'/'+ str(datetime.datetime.now())
	log = open(log_path+'.log','w')
	print '\nLog file will be generated: '+ log_path+'\n'

	makedirs([fastq_path,log_bs_dir,bam_path,sim_log_dir])
	
	## Setting path e arguments for each software #####
	path_bwa,bwa_args=Set_bwa(conf)
	path_picard,picard_args=Set_picard(conf)
	path_art,art_args=Set_art(conf)
	addsnv,addindel,bs_args=Set_bamsurgeon(conf)
	samtools_path=Set_samtools(conf)
	fasta=opts.fasta
	if opts.bed != None:
		# Cutting reference.fasta with bed to generate intervals.fasta
		fasta=getfasta(opts.ref,opts.bed,log)
		Index_fasta(fasta,samtools_path)


	if fasta != None:		
		# Generating paired ends fastq from given sequence.fasta
		path_pirs,opts_pirs=Set_pirs_args(conf)
		if opts.amplicon:
			Simulate_fastq_art(path_art,art_args,fasta,log)
		else:
			Simulate_fastq_pirs(path_pirs,opts_pirs,fasta,log)
		
		for file in os.listdir(opts.out_path+'/FASTQ'):
			if '1.fq' in file:
				fastq1=fastq_path+'/'+file
			elif '2.fq' in file:
				fastq2=fastq_path+'/'+file
		# Alignment of generated fastq in bam file, sorting and indexing
		out_bwa=Alignment_bwa(path_bwa,bwa_args,fastq1,fastq2,log)
		outbam=SamFormatConverter(path_picard,picard_args,out_bwa,log)
		status = subprocess.call("rm "+out_bwa, shell=True)
		bam=SortSam(path_picard,picard_args,outbam,log)
		status = subprocess.call("rm "+outbam , shell=True)
		Index_bam(path_picard,picard_args,bam,log)
		
	elif opts.fastq1 != None or opts.fastq2 != None :
		# if fastqs are given
		fastq1=opts.fastq1
		fastq2=opts.fastq2
		# Alignment of given fastq in bam file, sorting and indexing
		out_bwa=Alignment_bwa(path_bwa,bwa_args,fastq1,fastq2,log)
		outbam=SamFormatConverter(path_picard,picard_args,out_bwa,log)
		status = subprocess.call("rm "+out_bwa, shell=True)
		bam=SortSam(path_picard,picard_args,outbam,log)
		status = subprocess.call("rm "+outbam , shell=True)
		Index_bam(path_picard,picard_args,bam,log)

	elif opts.bam != None:
		# if bam file is given
		bam=opts.bam

	vars=opts.vars
	if opts.dbsnp !=None:
		vars=vars_from_db(opts.dbsnp,opts.num_snv_dbsnp,opts.num_indel_dbsnp,opts.out_path+'/Germline_dbsnp.vcf')
	if vars != None:
		print "\nStarting simulation of Germline variants:"
		# if a file containing variant to simulate is given
		#write variant from vcf format to bam surgeon format (bed format)
		snp,indel=Vcf_to_bamsurgeon(vars,opts.min,opts.max,opts.err)
		#if list contains snps
		if os.stat(snp).st_size != 0:
			outbam=bam_path + bam.split('/')[-1].split('.')[0]+'.snp.bam'
			os.chdir(log_bs_dir)
			Bam_surgeon(addsnv,bs_args,bam,outbam,snp,log)
			check_simul_vars(log_bs_dir,outbam,snp,sim_log_dir)
			print_vcf('Germline',log_bs_dir,sim_log_dir,opts.ref)
			#status = subprocess.call("rm "+bam , shell=True)
			#status = subprocess.call("rm "+bam + '.bai' , shell=True)
			bam=SortSam(path_picard,picard_args,outbam,log)
			status = subprocess.call("rm "+outbam, shell=True)
			Index_bam(path_picard,picard_args,bam,log)
		#if list contains indels
		if os.stat(indel).st_size != 0:

			outbam=bam_path + bam.split('/')[-1].split('.')[0]+'.indel.bam'
			os.chdir(log_bs_dir)
			Bam_surgeon(addindel,bs_args,bam,outbam,indel,log)
			check_simul_vars(log_bs_dir,outbam,indel,sim_log_dir)
			print_vcf('Germline',log_bs_dir,sim_log_dir,opts.ref)
			status = subprocess.call("rm "+bam , shell=True)
			status = subprocess.call("rm "+bam + '.bai' , shell=True)
			bam=SortSam(path_picard,picard_args,outbam,log)
			status = subprocess.call("rm "+outbam, shell=True)
			Index_bam(path_picard,picard_args,bam,log)

	germlinebam=bam_path + bam.split('/')[-1].split('.')[0]+'.Germline.bam'
	status = subprocess.call("mv "+bam+ ' '+ germlinebam , shell=True)
	status = subprocess.call("mv "+bam+ '.bai '+ germlinebam +'.bai', shell=True)
	bam=germlinebam
	print "Control sample: "+ germlinebam + '\n'

	#if analysis is somatic
	if opts.analysis == 'Somatic':
		vars_som=opts.vars_som
		if opts.cosmic != None:
			vars_som=vars_from_db(opts.cosmic,opts.num_snv_cosmic,opts.num_indel_cosmic,opts.out_path+'/Somatic_cosmic.vcf')
		print "Starting simulation of Somatic variants:"
		#write variant from vcf format to bam surgeon format (bed format)
		snp,indel=Vcf_to_bamsurgeon(vars_som,opts.min_som,opts.max_som,opts.err)
		somaticbam=bam.split('.')[0]+'.Somatic.bam'
		if os.stat(snp).st_size != 0:
			outbam=bam_path + bam.split('/')[-1].split('.')[0]+'.snp.Somatic.bam'
			os.chdir(log_bs_dir)
			Bam_surgeon(addsnv,bs_args,bam,outbam,snp,log)
			check_simul_vars(log_bs_dir,outbam,snp,sim_log_dir)
			print_vcf('Somatic',log_bs_dir,sim_log_dir,opts.ref)
			bam=SortSam(path_picard,picard_args,outbam,log)
			status = subprocess.call("rm "+outbam , shell=True)
			Index_bam(path_picard,picard_args,bam,log)
		
		if os.stat(indel).st_size != 0:
			outbam=bam_path + bam.split('/')[-1].split('.')[0]+'.indel.Somatic.bam'
			os.chdir(log_bs_dir)
			Bam_surgeon(addindel,bs_args,bam,outbam,indel,log)
			check_simul_vars(log_bs_dir,outbam,indel,sim_log_dir)
			print_vcf('Somatic',log_bs_dir,sim_log_dir,opts.ref)
			if '.sort' in bam:
				status = subprocess.call("rm "+bam , shell=True)
				status = subprocess.call("rm "+bam + '.bai' , shell=True)
			bam=SortSam(path_picard,picard_args,outbam,log)
			status = subprocess.call("rm "+outbam, shell=True)
			Index_bam(path_picard,picard_args,bam,log)

		status = subprocess.call("mv "+bam+ ' '+ somaticbam , shell=True)
		status = subprocess.call("mv "+bam+ '.bai '+ somaticbam +'.bai', shell=True)
		print "Somatic sample: "+ somaticbam
		print "Control sample: "+ germlinebam

	#subprocess.call("rm -r " + log_bs_dir, shell=True)





