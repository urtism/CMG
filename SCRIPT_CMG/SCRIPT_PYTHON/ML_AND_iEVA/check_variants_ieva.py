import argparse
import os
import vcf
import re




def gatk_extr(id,vcf_dir,chr,pos,ref,alt):
	try:
		GT,DP,AF,BQ = ['0','.','.','.']
		try:
			vcf = open(vcf_dir)
		except:
			print vcf_dir +' non presente'
			return GT,DP,AF,BQ
			

		for line in vcf:
			line=line.rstrip()
			if line.startswith('##'):
				continue
			elif line.startswith('#'):
				header = line.split('\t')
				if id in header:
					continue
				else:
					print id+' non trovato in '+ vcf_dir
					return GT,DP,AF,BQ
			else:
				var = line.split('\t')

			if [var[0],var[1],var[3],var[4]] == [chr,pos,ref,alt]:
				for field in var[header.index('INFO')].split(';'):
					if field.startswith('QD='):
						BQ = field.split('=')[1]
				
				GT,AD,DP,GQ,PL = var[header.index(id)].split(':')

				if GT != './.' and GT != '.':
					try:
						AF = str(float(AD.split(',')[1])/(float(AD.split(',')[1])+float(AD.split(',')[0])))
					except:
						AF = '.'
					vcf.close()
				return GT,DP,AF,BQ
			else:
				continue
		vcf.close()
		return GT,DP,AF,BQ
	except Exception as e:
		print id,chr,pos,ref,alt
		raise


def freebayes_extr(id,vcf_dir,chr,pos,ref,alt):
	try:
		GT,DP,AF,BQ = ['0','.','.','.']
		try:
			vcf = open(vcf_dir)
		except:
			print vcf_dir +' non presente'
			return GT,DP,AF,BQ

		for line in vcf:
			line=line.rstrip()
			if line.startswith('##'):
				continue
			elif line.startswith('#'):
				header = line.split('\t')
				if id in header:
					continue
				else:
					print id+' non trovato in '+ vcf_dir
					return GT,DP,AF,BQ
			else:
				var = line.split('\t')

			if [var[0],var[1],var[3],var[4]] == [chr,pos,ref,alt]:
				format =var[8].split(':')

				GT = var[header.index(id)].split(':')[0]
				if GT != './.' and GT != '.':
					AO = var[header.index(id)].split(':')[format.index('AO')]
					RO = var[header.index(id)].split(':')[format.index('RO')]
					QA = var[header.index(id)].split(':')[format.index('QA')]
					QR = var[header.index(id)].split(':')[format.index('QR')]
					DP = var[header.index(id)].split(':')[format.index('DP')]
					try:
						AF = str(float(AO)/(float(AO)+float(RO)))
					except:
						AF = '.'
					try:
						BQ = str(((float(QR)*float(RO)+float(QA)*float(AO))/((float(RO)+float(AO))*(float(RO)+float(AO)))))
					except:
						BQ = '.'
					vcf.close()
				return GT,DP,AF,BQ
			else:
				continue
		vcf.close()
		return GT,DP,AF,BQ
	except Exception as e:
		print id,chr,pos,ref,alt
		raise

def varscan_extr(id,vcf_dir,chr,pos,ref,alt):
	try:
		GT,DP,AF,BQ = ['0','.','.','.']
		try:
			vcf = open(vcf_dir)
		except:
			print vcf_dir +' non presente'
			return GT,DP,AF,BQ

		for line in vcf:
			line=line.rstrip()
			if line.startswith('##'):
				continue
			elif line.startswith('#'):
				header = line.split('\t')
				if id in header:
					continue
				else:
					print id+' non trovato in '+ vcf_dir
					return GT,DP,AF,BQ
			else:
				var = line.split('\t')


			if [var[0],var[1],var[3],var[4]] == [chr,pos,ref,alt]:

				GT,GQ,SDP,DP,RD,AD,FREQ,PVAL,RBQ,ABQ,RDF,RDR,ADF,ADR = var[header.index(id)].split(':')
				if GT != './.' and GT != '.':
					print GT
					FREQ = re.sub ('%','',FREQ)
					try:
						AF = str(float(re.sub(',','.',FREQ))/100.0)
					except:
						AF = '.'
					try:
						BQ = str((float(RBQ)*float(RD)+float(ABQ)*float(AD))/(float(RD)+float(AD)))
					except:
						BQ = '.'
					vcf.close()
				return GT,DP,AF,BQ
			else:
				continue
		vcf.close()
		return GT,DP,AF,BQ
	except Exception as e:
		print id,chr,pos,ref,alt
		raise

def samtools_extr(id,vcf_dir,chr,pos,ref,alt):
	try:
		GT,DP,AF,BQ = ['0','.','.','.']
		try:
			vcf = open(vcf_dir)
		except:
			print vcf_dir +' non presente'
			return GT,DP,AF,BQ

		for line in vcf:
			line=line.rstrip()
			if line.startswith('##'):
				continue
			elif line.startswith('#'):
				header = line.split('\t')
				if id in header:
					continue
				else:
					print id+' non trovato in '+ vcf_dir
					return GT,DP,AF,BQ
			else:
				var = line.split('\t')


			if [var[0],var[1],var[3],var[4]] == [chr,pos,ref,alt]:

				GT,PL,GP,GQ = var[header.index(id)].split(':')

				vcf.close()
				return GT,'','',''
			else:
				continue
		vcf.close()
		return GT,DP,AF,BQ
	except Exception as e:
		print id,chr,pos,ref,alt
		raise

def platypus_extr(id,vcf_dir,chr,pos,ref,alt):
	try:
		GT,DP,AF,BQ = ['0','.','.','.']
		try:
			vcf = open(vcf_dir)
		except:
			print vcf_dir +' non presente'
			return GT,DP,AF,BQ

		for line in vcf:
			line=line.rstrip()
			if line.startswith('##'):
				continue
			elif line.startswith('#'):
				header = line.split('\t')
				if id in header:
					continue
				else:
					print id+' non trovato in '+ vcf_dir
					return GT,DP,AF,BQ
			else:
				var = line.split('\t')


			if [var[0],var[1],var[3],var[4]] == [chr,pos,ref,alt]:

				GT,GL,GOF,GQ,NR,NV = var[header.index(id)].split(':')
				vcf.close()
				return GT,'','',''
			else:
				continue
		vcf.close()
		return GT,DP,AF,BQ
	except Exception as e:
		print id,chr,pos,ref,alt
		raise

def scalpel_extr(id,vcf_dir,chr,pos,ref,alt):
	try:
		GT,DP,AF,BQ = ['0','.','.','.']
		try:
			vcf = open(vcf_dir)
		except:
			print vcf_dir +' non presente'
			return GT,DP,AF,BQ

		for line in vcf:
			line=line.rstrip()
			if line.startswith('##'):
				continue
			elif line.startswith('#'):
				header = line.split('\t')
				if id in header:
					continue
				else:
					print id+' non trovato in '+ vcf_dir
					return GT,DP,AF,BQ
			else:
				var = line.split('\t')


			if [var[0],var[1],var[3],var[4]] == [chr,pos,ref,alt]:

				GT,AD,DP = var[header.index(id)].split(':')
				
				vcf.close()
				return GT,'','',''
			else:
				continue
		vcf.close()
		return GT,DP,AF,BQ
	except Exception as e:
		print id,chr,pos,ref,alt
		raise

def snver_extr(id,vcf_dir,chr,pos,ref,alt):
	try:
		GT,DP,AF,BQ = ['0','.','.','.']
		try:
			vcf = open(vcf_dir)
		except:
			print vcf_dir +' non presente'
			return GT,DP,AF,BQ

		for line in vcf:
			line=line.rstrip()
			if line.startswith('##'):
				continue
			elif line.startswith('#'):
				header = line.split('\t')
				if id in header:
					continue
				else:
					print id+' non trovato in '+ vcf_dir
					return GT,DP,AF,BQ
			else:
				var = line.split('\t')


			if [var[0],var[1],var[3],var[4]] == [chr,pos,ref,alt]:

				GT,PL = var[header.index(id)].split(':')[:2]
				
				return GT,'','',''
			else:
				continue

		vcf.close()
		return GT,DP,AF,BQ
	except Exception as e:
		print id,chr,pos,ref,alt
		raise


def main():

	parser = argparse.ArgumentParser('Controlla nei vcf le varianti del dataset e estrae DP, AF, BQ dal primo vcf che contiene la variante e quali caller la hanno trovata')
	parser.add_argument('--varlist',help="lista delle varianti da controllare")
	parser.add_argument('--cardiodir',help="dir che contiene le cartelle delle run con dentro i vcf Cardio")
	parser.add_argument('--cancerdir',help="dir che contiene le cartelle delle run con dentro i vcf Cancer")

	parser.add_argument('-o','--out',help="file di output")

	global opts
	opts = parser.parse_args()


	varlist = open(opts.varlist,'r')
	newvarlist = open(opts.out,'w')

	for var in varlist:
		var=var.rstrip()
		if var.startswith('Paziente') or var.startswith('Nome'):
			header = var.rstrip()
			newvarlist.write('\t'.join([header,'DP','AF','BQ','GATK','FREEB','VARSCAN','SAMTOOLS','PLATYPUS','SCALPEL']) + '\n')
		else:
			nome,id,chr,pos,ref,alt,HGVSc,HGVSp = var.split('\t')[:8]

			if 'Cardio' in id or 'Conn' in id:
				dir = opts.cardiodir
				pannello = 'Cardio'
			elif 'Cancer' in id:
				dir = opts.cancerdir
				pannello = 'Cancer'

			data = id.split('_')[0]

			for root,dirs,files in os.walk(dir):
				for rundir in dirs:

					if data in rundir:
						Scalpel_vcf = '/'.join([dir,rundir,id+'_Scalpel.norm.vcf'])
						Snver_vcf = '/'.join([dir,rundir,id+'_SNVer.norm.vcf'])
						Freebayes_vcf = '/'.join([dir,rundir,data+'_'+pannello +'_FreeBayes.norm.vcf'])
						Gatk_vcf = '/'.join([dir,rundir,data+'_'+pannello +'_GATK.norm.vcf'])
						Platypus_vcf = '/'.join([dir,rundir,data+'_'+pannello +'_Platypus.norm.vcf']) 
						Samtools_vcf = '/'.join([dir,rundir,data+'_'+pannello +'_Samtools.norm.vcf']) 
						Varscan_vcf = '/'.join([dir,rundir,data+'_'+pannello +'_VarScan.norm.vcf'])

						[Gatk_gt,Gatk_dp,Gatk_af,Gatk_qual] = gatk_extr(id,Gatk_vcf,chr,pos,ref,alt)
						[Freebayes_gt,Freebayes_dp,Freebayes_af,Freebayes_qual] = freebayes_extr(id,Freebayes_vcf,chr,pos,ref,alt)
						[Varscan_gt,Varscan_dp,Varscan_af,Varscan_qual] = varscan_extr(id,Varscan_vcf,chr,pos,ref,alt)
						[Samtools_gt,Samtools_dp,Samtools_af,Samtools_qual] = samtools_extr(id,Samtools_vcf,chr,pos,ref,alt)
						[Platypus_gt,Platypus_dp,Platypus_af,Platypus_qual] = platypus_extr(id,Platypus_vcf,chr,pos,ref,alt)
						[Snver_gt,Snver_dp,Snver_af,Snver_qual] = snver_extr(id,Snver_vcf,chr,pos,ref,alt)
						[Scalpel_gt,Scalpel_dp,Scalpel_af,Scalpel_qual] = scalpel_extr(id,Scalpel_vcf,chr,pos,ref,alt)

						if Gatk_gt != '0':
							print 'gatk'
							DP,AF,BQ = Gatk_dp,Gatk_af,Gatk_qual
						elif Freebayes_gt != '0':
							print 'freebayes'
							DP,AF,BQ = Freebayes_dp,Freebayes_af,Freebayes_qual
						elif Varscan_gt != '0':
							print 'varscan'
							DP,AF,BQ = Varscan_dp,Varscan_af,Varscan_qual
						else:
							DP,AF,BQ = ['.','.','.']

						newvarlist.write('\t'.join([var,DP,AF,BQ,Gatk_gt,Freebayes_gt,Varscan_gt,Samtools_gt,Platypus_gt,Scalpel_gt]) + '\n')
					else:
						continue


main()