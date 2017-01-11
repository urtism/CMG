import argparse
import re

def add_ann():
    out = open(opts.outfile + '.tsv','w')
    out_other = open(opts.outfile +'.Other_transcripts.tsv','w')
    vcf = open(opts.input,'r')
    file_list = open(opts.list,'r')
    
    ann_list = file_list.readlines()
    header_ann=[]
    ind_ann = []
    header=[]
    for line in vcf:
        line = line.rstrip()
        file_coor = open(opts.file,'r')
        vet_princ = []
        vet_other = []
        ann_list_other = []
        if line.startswith('##INFO=<ID=ANN') or line.startswith('##INFO=<ID=CSQ'):
            
            start = line.find('Allele')
            end = line.find('">')
            header_ann = line[start:end]
            header_ann = header_ann.split('|')
            continue
        elif line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            try:
                header_ann = line.split('\t')[0:5] + header_ann
            except:
                print line  
            for tag in ann_list:
                tag = tag.rstrip()
                if tag.startswith('#'):
                     continue
                ind_ann = ind_ann + [header_ann.index(tag)]
            var_princ = header_ann
        else:
            info=line.split('\t')[7]
            ann_split=info.split('ANN=')[1].split(',')
            ann_princ,ann_other = search_transcr(ann_split)
            var_princ = line.split('\t')[0:5] + ann_princ.split('|')
            
            for ann in ann_other:
                ann_list_other = ann_list_other + [line.split('\t')[0:5] + ann.split('|')]
            
            try:
                for i in var_princ:
                    if i == '':
                        var_princ[var_princ.index(i)] = '-'
            except:
                print 'gne'
            try:
                for ann in ann_list_other:
                    for i in ann:
                        if i == '':
                            ann[ann.index(i)] = '-'
            except:
                print 'gne'
        var=[]
        for variante in file_coor:
            var_split = variante.rstrip().split('\t')
            if var_split[0] == 'CHROM' and var_princ[0] == '#CHROM':
                for id in ind_ann:
                    try:
                        var = var + [var_princ[id]]
                    except:
                        continue
                out.write('\t'.join(var_split + var) + '\n')
                out_other.write('\t'.join(var_split + var) + '\n')
                
            elif var_princ[0].lstrip('#') == var_split[0] and var_princ[1] == var_split[1] and var_princ[3] == var_split[3] and var_princ[4] == var_split[4]:
                for id in ind_ann:
                    try:
                        var = var + [var_princ[id]]
                    except:
                        continue
                if  var != ['-']:
                    out.write('\t'.join(var_split + var) + '\n')
                else:    
                    for ann in ann_list_other:
                        if ann[0].lstrip('#') == var_split[0] and ann[1] == var_split[1] and ann[3] == var_split[3] and ann[4] == var_split[4]:
                            var=[]   
                            for id in ind_ann:
                                try:
                                    var = var + [ann[id]]
                                except:
                                    continue
                            out_other.write('\t'.join(var_split + var) + '\n')
        file_coor.close()                         
           
def search_transcr(ann_list):
    trsc=open(opts.trs_list,'r')
    trascritti=trsc.readlines()
    ann_other=[]
    ann_princ=''
    for ann in ann_list:
        if trascritti[0].rstrip() in ann.split('|') or trascritti[1].rstrip() in ann.split('|'):
            ann_princ = ann
        else:
            ann_other += [ann]
    return ann_princ,ann_other

def main():
    parser = argparse.ArgumentParser('aggiunge le annotazioni fornite nel file -f (una per riga) al file di input -i in un formato tab delimited in output -o')
    
    parser.add_argument('-i','--input',help="file di input in formato vcf")
    parser.add_argument('-l','--list',help="lista di annotazioni: una per riga")
    parser.add_argument('-f','--file',help="file tab delimited da cui pescare le coordinate del cromosoma")
    parser.add_argument('-t','--trs_list',help="lista di trascritti")
    parser.add_argument('-o','--outfile',help="file di output tab delimited")
    
    global opts
    
    opts = parser.parse_args()
    add_ann()

main()