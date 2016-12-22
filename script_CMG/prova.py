import statistics

vett_DP=[1,2,3,'',2]
nDP=0.0
v=[]
i=0
for dp in vett_DP:
	if dp and dp is not '':
		nDP= float(nDP)+float(dp)
		v=v+[float(dp)]
		i=i+1
try:
	media1=nDP/i
	media2= statistics.mean(v)
except:
	media1='.'
	media2='.'

print media1,media2
	

