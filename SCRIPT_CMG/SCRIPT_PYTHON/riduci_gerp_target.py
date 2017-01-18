import os
import glob2
import argparse






def main():
    parser = argparse.ArgumentParser('Adds class to the opts.dataset of variants.  Output is to stdout.')
    
    parser.add_argument('-t','--target',help="path del file dei target")
    parser.add_argument('-g','--gerp',help="path della cartella che contiene i gerpdb divisi per chrom")
    parser.add_argument('--tipo',help="Cardio, Cancer, Exome")
    parser.add_argument('-o','--out',help="path della cartella che contiene i gerpdb filtrati per i target divisi per chrom")
    

    
    global opts
    opts = parser.parse_args()
    
    target=open(opts.target,'r')
    target_list=[]
    for tar in target:
        tar=tar.rstrip()
        if tar.startswith('@'):
            continue
        else:
            target_list=target_list + [tar]
    out_total=open(opts.out + '/gerp.total.'+ opts.tipo +'.rates','w')
    for filename in glob2.glob(os.path.join(opts.gerp,'*.rates')):
        chrom=(filename.split('/')[-1]).split('.')[0]
        print chrom
        #chrom='chrY'
        outdb=open(opts.out + '/gerp.' + chrom +'.'+ opts.tipo +'.rates','w')
        
        gerpchr=open(filename,'r')
        gerp_list=gerpchr.readlines()
#         for elem in gerp_list:
#             i=i+1
#             #print str(i),elem
        
        for targ in target_list:
            print targ
            targ_split= targ.split('\t')
            
            if chrom == targ_split[0]:
                i=int(targ_split[1])
                for line in gerp_list[int(targ_split[1])-1:int(targ_split[2])]:
                    i=i+1
                    outdb.write('\t'.join([chrom,str(i),line]))
                    out_total.write('\t'.join([chrom,str(i),line]))
        
        gerpchr.close()
            
        
main()