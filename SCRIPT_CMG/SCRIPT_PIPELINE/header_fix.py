import argparse
import sys
#from tr import tr
import re


parser = argparse.ArgumentParser('Parse VCF HEADER from FreeBayes,VarScan,GATK to fix it.  Output is to stdout.')
parser.add_argument('-f', '--file', help="freebayes vcf output file name")
parser.add_argument('-v', '--variantcaller', help="variant caller: F = freebayes, G = GATK, V = Varscan, P = Platypus, S = Samtools, L = Scalpel, N = SNVer" )

opts = parser.parse_args()

read=open(opts.file)


for line in read:
	if opts.variantcaller == 'F':
		
		if line.startswith('##FORMAT=<ID=DPR'):
			riga=(line.split(','))
			riga[riga.index('Number=A')]='Number=G'
			line=','.join(riga)
			
		if line.startswith('##FORMAT=<ID=GQ'):
			riga=(line.split(','))
			try:
		 		riga[riga.index('Type=Integer')]='Type=Float'
		 	except:
		 		pass
			line=','.join(riga)
	# 	if line.startswith('chr'):
	# 		riga=line.split('\t')
	# 		replaced = re.sub('[(chr)]', '', riga[0])
	# 		riga[0]=replaced
	# 		line='\t'.join(riga)
	
	elif opts.variantcaller == 'G':
	
		if line.startswith('##FORMAT=<ID=AD'):
			riga=line.split(',')
			try:
				riga[riga.index('Number=.')]='Number=G'
				line=','.join(riga)
			except:
				pass

	elif opts.variantcaller == 'V':
	
		if line.startswith('chr'):
			riga=line.split('\t')
			ref=riga[3]
			alt=riga[4]
			if '/+' in ref:
				riga[3]=alt
				riga[4]=alt+ref.split('/+')[1]
				riga[7]+=";FIX"
			elif '/' in ref:
				riga[3]=ref.split('/')[0]
				riga[7]+=";FIX"
				
			if '/-' in alt:
				riga[4]=alt.split('/-')[0]
				riga[7]+=";FIX"
			elif '/' in alt:
				riga[4]=alt.split('/')[0]
				riga[7]+=";FIX"

			line='\t'.join(riga)

	elif opts.variantcaller == 'P':
	
		# if line.startswith('##FORMAT=<ID=GL'):
		# 	riga=(line.split(','))
		# 	try:
		# 		riga[riga.index('Number=.')]='Number=G'
		# 		line=','.join(riga)
		# 	except:
		# 		pass
		
		if line.startswith('##FORMAT=<ID=NR'):
			try:
				riga=line.split(',')
				riga[riga.index('Number=.')]='Number=A'
				line=','.join(riga)
			except:
				pass

		if line.startswith('##FORMAT=<ID=NV'):
			try:
				riga=line.split(',')
				riga[riga.index('Number=.')]='Number=A'
				line=','.join(riga)
			except:
				pass

		if line.startswith('##INFO=<ID=FR'):
			try:
				riga=line.split(',')
				riga[riga.index('Number=.')]='Number=A'
				line=','.join(riga)
			except:
				pass

		if line.startswith('##INFO=<ID=PP'):
			try:
				riga=line.split(',')
				riga[riga.index('Number=.')]='Number=A'
				line=','.join(riga)
			except:
				pass

		if line.startswith('##INFO=<ID=TR'):
			try:
				riga=line.split(',')
				riga[riga.index('Number=.')]='Number=A'
				line=','.join(riga)
			except:
				pass

		if line.startswith('##INFO=<ID=NF'):
			try:
				riga=line.split(',')
				riga[riga.index('Number=.')]='Number=A'
				line=','.join(riga)
			except:
				pass

		if line.startswith('##INFO=<ID=NR'):
			try:
				riga=line.split(',')
				riga[riga.index('Number=.')]='Number=A'
				line=','.join(riga)
			except:
				pass
		if line.startswith('chr'):
			riga=line.rstrip().split('\t')
			for s in riga[9:]:
				sample = s.split(':')
				gt = sample[0]

				if gt == '1/0':
					gt = '0/1'

				sample[0] = gt
				riga[riga.index(s)] = ':'.join(sample)

			line='\t'.join(riga)

	elif opts.variantcaller == 'S':
	
		if line.startswith('chr'):
			riga=line.split('\t')

			if riga[3].isupper() and riga[4].isupper():
				line='\t'.join(riga)
				pass
			else:
				try:
					riga[3] = riga[3].upper()
					riga[4] = riga[4].upper()
					line='\t'.join(riga)
				except:
					pass
			for s in riga[9:]:
				sample = s.split(':')
				gt = sample[0]
				if gt == '1/0':
					gt = '0/1'
				sample[0] = gt
				riga[riga.index(s)] = ':'.join(sample)
			line='\t'.join(riga)

	elif opts.variantcaller == 'L':
	
		if line.startswith('##FORMAT=<ID=AD'):
			riga=(line.split(','))
			try:
				riga[riga.index('Number=.')]='Number=A'
				line=','.join(riga)
			except:
				pass

	elif opts.variantcaller == 'N':
	
		if line.startswith('##FORMAT=<ID=PL'):
			riga=(line.split(','))
			try:
				riga[riga.index('Number=1')]='Number=G'
				line=','.join(riga)
			except:
				pass
		elif line.startswith('##INFO=<ID=FS') or line.startswith('##INFO=<ID=SP') or line.startswith('##INFO=<ID=PV'):
			riga=(line.split(','))
			try:
				riga[riga.index('Number=.')]='Number=1'
				line=','.join(riga)
			except:
				pass
		elif line.startswith('#CHROM'):
			riga = (line.split('\t'))
			sample = riga[-1].split('/')[-1].split('.')[0]
			riga[-1] = sample
			line='\t'.join(riga)
		elif line.startswith('##'):
			pass
		else:
			riga = (line.rstrip().split('\t'))

			INFO = riga[7].split(';')
			fixedINFO = []

			for field in INFO:
				if field.startswith('FS=') or field.startswith('SP='):
					fixedfield = '.'.join(field.split(','))
				else:
					fixedfield = field
				fixedINFO += [fixedfield] 
			riga[7] = ';'.join(fixedINFO)

			for s in riga[9:]:
				sample = s.split(':')
				gt = sample[0]
				if gt == '1/0':
					gt = '0/1'
				sample[0] = gt
				riga[riga.index(s)] = ':'.join(sample)
			line='\t'.join(riga)

	print line.rstrip()

read.close()