file1=open('/home/jarvis/Scrivania/snp_trap.list','r')
out=open('/home/jarvis/Scrivania/snp_trap.adjust.list','w')


for line in file1:
		line=line.rstrip()
		if line.startswith('Locus'):
			out.write('\t'.join(["Locus","Chromosome","pos","ref","alt","MAF 1000g","MAF","pAA","pAB","pBB"])+'\n')
		else:
			rs=line.split('\t')[0]
			chr=line.split('\t')[1]
			pos=line.split('\t')[2]
			ref=line.split('\t')[3]
			alt=line.split('\t')[4]
			maf=line.split('\t')[5]

			base_maf=maf.split('=')[0]
			if base_maf == ref:
				freq=1-float((maf.split('=')[1]).split('/')[0])
			else:
				freq=float((maf.split('=')[1]).split('/')[0])
			
			pAA=(1-freq)**2
			pAB=2*(1-freq)*freq
			pBB=freq**2
			out.write('\t'.join([rs,chr,pos,ref,alt,maf,str(freq),str(pAA),str(pAB),str(pBB)])+'\n')