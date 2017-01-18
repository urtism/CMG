import argparse

def add_info_paz(in_file,paz):
	out = open(opts.outfile,'w')
	pazienti= {}
	for paziente in paz:
		paziente = paziente.rstrip()
		if paziente.startswith('NOME_PAZ') or paziente == '':
			continue
		else:
			nome = paziente.split('\t')[0]
			num = paziente.split('\t')[1]
			run = paziente.split('\t')[2]
			id = paziente.split('\t')[4]
			pazienti[nome] = [num,run,id]
	index_sanger = 0
	for line in in_file:
		line = line.rstrip()
		var = line.split('\t')
		HGVSc = '-'
		HGVSp = '-'
		if line.startswith('NOME_PAZ'):
			index_sanger = var.index('HGVS')
			var[index_sanger:index_sanger] = ['HGVSc','HGVSp']
			var.remove(var[index_sanger +2])
			var[int(opts.nome)+1:int(opts.nome)+1] = ['NUM_PAZ','NUM_RUN','ID_PAZ']
			out.write('\t'.join(var) + '\n')
		else:
			nome = var[int(opts.nome)]
			sanger_split = var[index_sanger].split(' ')
			for el in sanger_split:
				if el.startswith('c.'):
					HGVSc = el
				elif el.startswith('p.'):
					HGVSp = el
			var[index_sanger:index_sanger] = [HGVSc,HGVSp]
			del var[index_sanger+2]
			try:
				var[int(opts.nome)+1:int(opts.nome)+1] = pazienti[nome]
			except:
				var[int(opts.nome)+1:int(opts.nome)+1] = ['PAZIENTE NON TROVATO','-','-']
			out.write('\t'.join(var) + '\n')

def main():
	parser = argparse.ArgumentParser('aggiunge il num_paziente, il num_run e il paziente id al file in ingresso -v')
	parser.add_argument('-f','--infile',help="lista di varianti a cui aggiungere id_paziente")
	parser.add_argument('-p','--paz',help="lista dei pazienti illumina")
	parser.add_argument('-o','--outfile',help="file di uscita")
	parser.add_argument('-n','--nome',help="indice della colonna dove e' presente il nome paziente partendo da 0")
	global opts
	opts = parser.parse_args()
	pazienti = open(opts.paz,'r')
	in_file = open(opts.infile,'r')
	add_info_paz(in_file,pazienti)

main()