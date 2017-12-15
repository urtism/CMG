import re
import argparse
import sys
import statistics
import os
import scipy.stats as stats


def extract_gaps(gaps_file,sample,run,gaps_list,reg,zone,gaps_list_sample):
    gaps = open(gaps_file,'r')
    
    for gap in gaps:
        
        if gap.startswith('#'):
            continue
        else:
            gap_split=gap.split(',')
            chrom=gap_split[0]
            start=gap_split[1]
            stop=gap_split[2]
            reg_id=gap_split[3]
            mean_gap_dp=gap_split[4]
            region_interval=gap_split[5]
            gap_interval=gap_split[6].rstrip()
            gap_lenght=int(stop)-int(start)
            gap_sample='\t'.join([run,sample,reg_id])
            zona='\t'.join([chrom,start,stop])
            
            try:
                if int(start)>=int(reg_id.split('.')[2]) and int(start) < int(reg_id.split('.')[3]) or int(stop)>int(reg_id.split('.')[2]) and int(stop)<=int(reg_id.split('.')[3]):
                    gaps_list['\t'.join([run,sample,chrom,start,stop])]=[reg_id,mean_gap_dp,str(gap_lenght+1)]
                    
                    if gap_sample in gaps_list_sample.keys():
                        num = gaps_list_sample.get(gap_sample)
                        gaps_list_sample[gap_sample] = num + 1
                    else:
                        gaps_list_sample[gap_sample] = 1
                    
                    if reg_id in reg.keys(): 
                        if gaps_list_sample[gap_sample] == 1:
                            num = reg.get(reg_id)
                            reg[reg_id]=num+1
        
                    else:
                        reg[reg_id]=1
                    
                    if zona in zone.keys():
                        num = zone.get(zona)
                        zone[zona]=num+1
                    else:
                        zone[zona]=1 
            except:
                gaps_list['\t'.join([run,sample,chrom,start,stop])]=[reg_id,mean_gap_dp,str(gap_lenght+1)]
                    
                if gap_sample in gaps_list_sample.keys():
                    num = gaps_list_sample.get(gap_sample)
                    gaps_list_sample[gap_sample] = num + 1
                else:
                    gaps_list_sample[gap_sample] = 1
                
                if reg_id in reg.keys(): 
                    num = reg.get(reg_id)
                    reg[reg_id]=num+1
    
                else:
                    reg[reg_id]=1
                
                if zona in zone.keys():
                    num = zone.get(zona)
                    zone[zona]=num+1
                else:
                    zone[zona]=1 
    
    gaps.close()       
            
        


def main():
    
    parser = argparse.ArgumentParser('Parse VCF output from Variant callers to output a variant_dataset.txt.  Output is to stdout.')
    parser.add_argument('-i', '--gaps_path', help="path della cartella da cui estrarre i gaps (per esempio Cardio)")
    parser.add_argument('-o','--out',help="path di output")
    
    global opts 
    opts = parser.parse_args()
    gaps_hash = {}
    regions_hash = {}
    gap_zone_hash = {}
    gaps_list_sample ={}
    i=0
    for dir_run in os.listdir(opts.gaps_path):
        #print '\nrun',dir_run
        dir_split = dir_run.split('_')
        nome_run = '_'.join(dir_split[0:4])
        if dir_run == 'Versioni_2.0':
            continue
        for dir_sample in os.listdir(opts.gaps_path+'/'+dir_run):
            if dir_sample.endswith('pdf'):
                continue
            #print 'sample',dir_sample
            gaps_name = '-'.join(dir_sample.split('_'))+'_S1.gaps.csv'
            #print 'gaps_name',gaps_name
            
            for sample_file in os.listdir(opts.gaps_path+'/'+dir_run + '/'+ dir_sample):
                if sample_file.endswith('gaps.csv'):
                    i+=1
                    old_path = opts.gaps_path+'/'+dir_run + '/'+ dir_sample+'/'+ sample_file
                    new_path = opts.gaps_path+'/'+dir_run + '/'+ dir_sample+'/'+ gaps_name
                    os.rename(old_path, new_path)
                    print 'Analizzo Run:'+nome_run,'Gap: '+ gaps_name
                    extract_gaps(new_path,dir_sample,nome_run,gaps_hash,regions_hash,gap_zone_hash,gaps_list_sample)
                    
                    
    out = open(opts.out,'w')
    print 'Printo i dati...'
    out.write('\t'.join(['RUN_ID','SAMPLE_ID','CHROM', 'START', 'STOP','REGION_ID','MEAN_DP','GAP_LENGHT','REGION_FREQ','GAP_FREQ','GAP_IN_REG_ID'])+ '\n')
    for el in gaps_hash.keys():
        reg_id=gaps_hash.get(el)[0]
        num=regions_hash.get(reg_id)
        zona='\t'.join(el.split('\t')[2:])
        num_gap=gap_zone_hash.get(zona)
        num_gap_sample = gaps_list_sample.get('\t'.join([el.split('\t')[0],el.split('\t')[1],reg_id]))
        out.write(el + '\t' + '\t'.join(gaps_hash.get(el)+[str(num)+'/'+str(i)]+[str(num_gap)+'/'+str(i)]+[str(num_gap_sample)]) +  '\n')
        print reg_id + ' DONE'

    
                
main()
