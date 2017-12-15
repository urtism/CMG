import argparse
import re

def estrai_annotazione():
    vcf = open(opts.vcf,'r')
    out = open(opts.out,'w')
    tsv = open(opts.tsv,'r')
    tags = open(opts.tag_list,'r')
    var_list = dict()
    tag_list=[]
    
    for var in tsv:
        var = var.rstrip()

        if var.startswith('CHROM'):
            header=var

        else:
            chrom = var.split('\t')[0]
            pos = var.split('\t')[1]
            id = var.split('\t')[2]
            ref = var.split('\t')[3]
            alt =  var.split('\t')[4]
            
            var_id = '\t'.join([chrom,pos,ref,alt])
            
            var_list[var_id] = var.split('\t')

    for tag in tags:
        tag=tag.rstrip()
        if tag.startswith('#'):
            continue
        else:
            tag_list+=[tag]

    tags.close()
    header_tot = header.split('\t') + tag_list
    out.write('\t'.join(header_tot)+'\n')

    for line in vcf:
        line = line.rstrip()

        if line.startswith('##INFO=<ID=ANN') or line.startswith('##INFO=<ID=CSQ'):
            start = line.find('Allele')
            end = line.find('">')
            header_ann = (line[start:end]).split('|')
            continue
        elif line.startswith('#'):
            continue
            
        else:
            chrom = line.split('\t')[0]
            pos = line.split('\t')[1]
            id = line.split('\t')[2]
            ref = line.split('\t')[3]
            alt =  line.split('\t')[4]
            var_id = '\t'.join([chrom,pos,ref,alt])
            
            if var_id in var_list.keys():
                info = (line.split('\t')[7]).split('=')[1]
                info_split = info.split('|')

                tags = var_list.get(var_id)
                
                for tag in tag_list:              
                    new_tag = info_split[header_ann.index(tag)]               
                    if new_tag == '':
                        new_tag = '.'
                    tags = tags + [new_tag]
                
                var_list[var_id] = tags
                out.write('\t'.join(var_list.get(var_id))+'\n')

    vcf.close()
    out.close()
    tsv.close()



def split_annotazione():
    vcf = open(opts.vcf,'r')
    trsc = open(opts.trs_list,'r')
    out = open(opts.out+'.vcf','w')
    out_other = open(opts.out+'.other.vcf','w')

    trascritti=trsc.readlines()

    for line in vcf:
        line = line.rstrip()
        line_split = line.split('\t')
        trovato = 0
        if line.startswith('#'):
            out.write(line+'\n')
            if line.startswith('##INFO=<ID=ANN') or line.startswith('##INFO=<ID=CSQ'):
                continue
            else:
                out_other.write(line+'\n')
        else:
            info=line_split[7]
            for ann in info.split('ANN=')[1].split(','):
                for trascritto in trascritti:
                    if trascritto.rstrip() in ann.split('|'):
                        trovato=1
                        line_split[7] = 'ANN='+ ann
                        out.write('\t'.join(line_split)+'\n')
                        break
            if trovato == 0:
                line_split[7]='.'
                out_other.write('\t'.join(line_split)+'\n')

    vcf.close()
    trsc.close()
    out.close()
    out_other.close()



def main():
    parser = argparse.ArgumentParser('aggiunge le annotazioni fornite nel file -f (una per riga) al file di input -i in un formato tab delimited in output -o')
    
    parser.add_argument('-i','--vcf',help="file delle varianti annotate in formato vcf")
    parser.add_argument('-l','--tag_list',help="lista di annotation features da aggiungere: una per riga",default=None)
    parser.add_argument('-f','--tsv',help="file tab delimited a cui aggiungere le annotation features",default=None)
    parser.add_argument('-t','--trs_list',help="lista di trascritti su cui splittare le annotazioni",default=None)
    parser.add_argument('-o','--out',help="file di output")
    parser.add_argument('-S', '--split', help="abilita lo split delle varianti per i trascritti contenuti in list", action='store_true')
    
    global opts
    
    opts = parser.parse_args()
    
    if opts.split:
        split_annotazione()
    else:
        estrai_annotazione()

main()