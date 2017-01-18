import re
import string
import argparse
import sys
import statistics
import os
import scipy.stats as stats


''' //////////// CLASSI ////////////'''

class Caller():
    GT=''
    AO=''
    RO=''
    AO_f=''
    AO_r=''
    DP_f=''
    DP_r=''
    DP=''
    QB=''
    Call='0'
    AF=''
    FILTER=''
    GQ=''
    

class Freebayes(Caller):
       
    #PROVENIENTI DAL FORMAT
    A_QB=''
    R_QB=''

class GATK(Caller):    
    pass
    

class Varscan(Caller):
    QB_A=''

class Features():

    GT_Varscan='.'
    GT_Freebayes='.'
    GT_GATK='.'
    
    DP_Varscan='.'
    DP_Freebayes='.'
    DP_GATK='.'
    
    DP_media='.'
    DP_mediana='.'

    QB_GATK='.'
    QB_Varscan='.'
    QB_Freebayes='.'

    AF_GATK='0'
    AF_Varscan='0'
    AF_Freebayes='0'
    AF_media='.'
    AF_mediana='.'
    
    CallGATK='0'
    CallVarscan='0'
    CallFreebayes='0'
    
    FILTER_GATK='.'
    FILTER_Varscan='.'
    FILTER_Freebayes='.'

''' //////////// FUNZIONI ////////////'''

def get_info_Freebayes(chrom,pos,ref,alt,filter,info,format,sample,freebayes):
    '''estrae le informazioni dal vcf di freebayes'''
    
    freebayes.GT=sample[format.index('GT')]
    if freebayes.GT=='.' :
        freebayes.GT='./.'
    else:
        freebayes.GQ=sample[format.index('GQ')]
        freebayes.AO=float(sample[format.index('AO')])
        freebayes.RO=float(sample[format.index('RO')])
        freebayes.DP=float(sample[format.index('DP')])
        freebayes.AF=freebayes.AO/freebayes.DP
        try:
            freebayes.A_QB=float(sample[format.index('QA')])/freebayes.AO
        except:
            freebayes.A_QB='.'
        try:    
            freebayes.R_QB=float(sample[format.index('QR')])/freebayes.RO
        except:
            freebayes.R_QB='.'
        try:    
            freebayes.QB=(float(sample[format.index('QA')]) + float(sample[format.index('QR')]))/freebayes.DP
        except:
            if sample[format.index('QA')] == '0':
                freebayes.QB=freebayes.R_QB
            elif sample[format.index('QR')] == '0':
                freebayes.QB=freebayes.R_QA

        freebayes.FILTER=filter
    
        freebayes.Call=1

def get_info_GATK(chrom,pos,ref,alt,filter,info,format,sample,GATK):
    '''estrae le informazioni dal vcf di GATK'''
    
    GATK.GT=sample[format.index('GT')]
    if GATK.GT=='./.':
        pass
    else:
        GATK.AO=float((sample[format.index('AD')]).split(',')[1])
        GATK.RO=float((sample[format.index('AD')]).split(',')[0])
        GATK.DP=float(sample[format.index('DP')])
        GATK.AF=GATK.AO/GATK.DP        
        GATK.QB=sample[format.index('SQD')]
        GATK.FILTER=filter
        GATK.Call=1
    
def get_info_varscan(chrom,pos,ref,alt,filter,info,format,sample,varscan):
    '''estrae le informazioni dal vcf di varscan'''
    
    varscan.GT=sample[format.index('GT')]
    if varscan.GT== './.':
        #print chrom,pos,ref,alt,varscan.GT
        pass
    else:
        varscan.AO=float(sample[format.index('AD')])
        varscan.RO=float(sample[format.index('RD')])
        varscan.DP=varscan.RO + varscan.AO
        varscan.Call=1
        varscan.QB_R=float(sample[format.index('RBQ')])
        varscan.QB_A=float(sample[format.index('ABQ')])
        Varscan.GQ=float(sample[format.index('GQ')])
        varscan.AF=float(varscan.AO/(varscan.DP))
        Varscan.FILTER=filter
    Varscan.Call=1

def set_features(dictionary):
    '''setta i valori delle features in base alle info estratte dai vcf'''
    for variante in dictionary.keys():
        features=Features()
        varc_array=dictionary.get(variante)
        
        vett_MBQ=[]
        vett_DP=[]
        vett_AO=[]
        vett_RO=[]
        index=0
        #print variante
        #print varc_array
        for varcall in varc_array:
            if varcall is not "":
                vett_MBQ=vett_MBQ+[varcall.QB]
                vett_DP=vett_DP+[varcall.DP]
                vett_AO=vett_AO+[varcall.AO]
                vett_RO=vett_RO+[varcall.RO]
                
                if index == 0:
     
                    features.GT_Freebayes=varc_array[0].GT
                    features.AO_Freebayes=varc_array[0].AO
                    features.RO_Freebayes=varc_array[0].RO
                    features.DP_Freebayes=varc_array[0].DP
                    features.QB_Freebayes=varc_array[0].QB
                    features.GQ_Freebayes=varc_array[0].GQ
                    features.CallFreebayes=varc_array[0].Call
                    features.AF_Freebayes=varc_array[0].AF
                    features.FILTER_Freebayes=varc_array[0].FILTER

                if index == 2:

                    features.GT_GATK=varc_array[2].GT
                    features.AO_GATK=varc_array[2].AO
                    features.RO_GATK=varc_array[2].RO
                    features.DP_GATK=varc_array[2].DP
                    features.QB_GATK=varc_array[2].QB
                    features.GQ_GATK=varc_array[2].GQ
                    features.CallGATK=varc_array[2].Call
                    features.AF_GATK=varc_array[2].AF                    
                    features.FILTER_GATK=varc_array[2].FILTER

                elif index == 1:
                    
                    features.GT_Varscan=varc_array[1].GT
                    features.AO_Varscan=varc_array[1].AO
                    features.RO_Varscan=varc_array[1].RO
                    features.DP_Varscan=varc_array[1].DP
                    features.QB_Varscan=varc_array[1].QB_A
                    features.GQ_Varscan=varc_array[1].GQ
                    features.CallVarscan=varc_array[1].Call
                    features.AF_Varscan=varc_array[1].AF
                    features.FILTER_Varscan=varc_array[1].FILTER


            index = index + 1    
        vett_AF_media=[features.AF_GATK,features.AF_Varscan,features.AF_Freebayes]
        AF_media=0
        AF_median=0
        
        v=[]
        for dp in vett_DP:
            if dp and dp is not '':
                v=v+[float(dp)]
        try:
            features.DP_media=statistics.mean(v)
        except:
            features.DP_media='.'
        try:
            features.DP_mediana=statistics.median(v)
        except:
            features.DP_mediana='.'        
        
        v=[]
        for af in vett_AF_media:
            if af and af is not '0':
                v=v+[float(af)]
        try:
            features.AF_media= statistics.mean(v)
        except:
            features.AF_media='.'
        try:
            features.AF_mediana= statistics.median(v)
        except:
            features.AF_mediana='.'    
        
        v=[]
        for bq in vett_MBQ:
            if bq and bq is not '.':
                v=v+[float(bq)]
        try:
            features.MBQ_media=statistics.mean(v)
        except:
            features.MBQ_media='.'
        try:
            features.MBQ_mediana=statistics.median(v)
        except:
            features.MBQ_mediana='.'    

        dictionary[variante]= varc_array + [features]

def switch(dictionary,ID,index,chrom,pos,ref,alt,filter,info,format,sample):
    '''tramite index richiama la funzione di estrazione delle informazioni del variant caller associato all'indice'''
    if dictionary.has_key(ID):
        vettore=dictionary[ID]
    else:
        vettore=['','','']

    if index==0:
        # print 'Freebayes'
        freebayes=Freebayes()
        get_info_Freebayes(chrom,pos,ref,alt,filter,info,format,sample,freebayes)
        vettore[0]=freebayes

    elif index==2:
        # print 'gatk'
        gatk=GATK()
        get_info_GATK(chrom,pos,ref,alt,filter,info,format,sample,gatk)
        vettore[2]=gatk
    elif index==1 :
        # print 'varscan'
        varscan=Varscan()
        get_info_varscan(chrom,pos,ref,alt,filter,info,format,sample,varscan)
        vettore[1]=varscan
    dictionary[ID]=vettore

def read(iterable,index,dictionary):
    '''legge il vcf e splitta le varie sezioni'''
    for line in iterable:
        line = line.rstrip()
        if line.startswith('#'):
            continue
        else:
            ind=0
            var = line.split("\t")
            #print var[1]
            chrom=var[0]
            pos=var[1]
            ref=var[3]
            alt=var[4]
            filter=var[6]
            ID='\t'.join([chrom,pos,ref,alt])
            INFO=var[7]
            info=INFO.split(";") 
            FORMAT=var[8]
            format=FORMAT.split(":")
            SAMPLE = var[-1]
            sample = SAMPLE.split(':')
            if alt != '*':
                switch(dictionary,ID,index,chrom,pos,ref,alt,filter,info,format,sample)
                
def control(dictionary):
    ''' esegue un controllo sulle varianti, se non hanno variant caller che le chiama vengono eliminate'''
    for variante in dictionary.keys():
        if dictionary[variante][:3] == ['','','']:
            del dictionary[variante]
                
def print_var(dictionary,out,sample_name):

    lista_features=open(opts.listaFeatures,'r')
    dataset_varianti=open(out + '/' + sample_name + '.tsv','w')
    header=[]
    features_variante=[]
    
    for line in lista_features:
        line = line.rstrip()
        if line.startswith('#'):
            continue
        else:
            header=header+[line]
            features_variante=features_variante+['features.'+line]

    #print features_variante
    dataset_varianti.write('CHROM\tPOS\tID\tREF\tALT\t' + '\t'.join(header)+ '\n')
    for variante in dictionary.keys():
        features = dictionary.get(variante)[-1]
        #print features_variante
        features_variante_eval=[]
        for feat in features_variante:
            #print feat
             feat_eval=str(eval(feat))
             features_variante_eval=features_variante_eval + [feat_eval]
        #CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT    20151202_01_Cardio
        var=variante.split('\t')[0]+'\t'+variante.split('\t')[1]+'\t' +sample_name +'\t'+variante.split('\t')[2]+'\t'+variante.split('\t')[3]+ '\t' + '\t'.join(features_variante_eval)
        dataset_varianti.write(var+ '\n')
#         for gt in [features.GT_GATK,features.GT_Varscan,features.GT_Freebayes]:
#             if gt != './.' and gt != '.':
#                 dataset_varianti.write(var+ '\n')
#                 break
    dataset_varianti.close()
    lista_features.close()
    
def print_vcf(varianti,out):
    dataset_varianti_vcf=open(out+ '/TOTAL.vcf','w')
    dataset_varianti_vcf.write('##fileformat=VCFv4.2\n'+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLES\n')
    for variante in varianti.keys():
        var_vcf=variante.split('\t')[0]+'\t'+variante.split('\t')[1]+'\t.\t'+variante.split('\t')[2]+'\t'+variante.split('\t')[3]+'\t.\t.\t.\t.\t.'
        dataset_varianti_vcf.write(var_vcf+ '\n')
    dataset_varianti_vcf.close()


def samples_name_extract(vcf):
    samples = []
    for line in vcf:
        line=line.rstrip()
        if line.startswith('#CHROM'):
            line_split = line.split('\t')
            samples = line_split[9:]
    return samples
        
def split_vcf(vcf_dir,samples):
    vcf_name = os.path.basename(vcf_dir)
    vcf_path = os.path.dirname(vcf_dir)
    vcf = open(vcf_dir,'r')
    header = []
    header_chrom = []
    varianti = []
    
    if 'Free' in vcf_name:
            variant_caller = 'FreeBayes'
    elif 'GATK' in vcf_name:
            variant_caller = 'GATK'
    elif 'VarScan' in vcf_name:
            variant_caller = 'VarScan'

    for line in vcf:
        line=line.rstrip()
        if line.startswith('##'):
             header = header + [line]
        elif line.startswith('#CHROM'):
             header_chrom = line.split('\t')
        else:
             varianti = varianti + [line]
    i=0
    for sample in samples:
        print sample
        try: 
            os.mkdir(opts.out_path +'/' + sample)
        except:
            pass
        
        sample_vcf = open(opts.out_path +'/' + sample +'/' + sample + '_'+variant_caller +'.vcf','w')
        
        sample_vcf.write('\n'.join(header) + '\n')
        sample_vcf.write('\t'.join(header_chrom[0:9] + [sample])  +'\n')
        gvcf = open(opts.gvcf_path +'/' + sample +'.g.vcf','r')
        all_gvcf=gvcf.readlines()
        gvcf.close()

        sample_gvcf = []

        for riga in all_gvcf:
            riga = riga.rstrip()
            if riga.startswith('#'):
                continue
            else:
                riga_split = riga.split('\t')
                if riga_split[4] != '<NON_REF>':
                    sample_gvcf += [riga]

        for variante in varianti:
            
            variante_split = variante.split('\t')
            chrom = variante_split[0]
            pos = variante_split[1]
            id = variante_split[2]
            ref = variante_split[3]
            alt = variante_split[4]
            if variant_caller == 'GATK':
                sSB = '.'
                sQD = '.'
                for line in sample_gvcf:
                    line = line.rstrip()
                    if line.startswith(chrom+'\t'+pos) and line.split('\t')[4]:
                        format = line.split('\t')[-2]
                        sformat = line.split('\t')[-1]
                        format_qual = line.split('\t')[5]
                        info = line.split('\t')[7].split(';')
                        for elem in info:
                            if elem.startswith('DP='):
                                DP = elem.split('=')[1]
                                break
                        try:
                            sSB = sformat.split(':')[format.split(':').index('SB')]
                            sQD = round(float(format_qual)/float(DP) ,2)
                            break
                        except:
                            print 'ci sono problemi',sample,chrom,pos,format,sformat
        
                variante_common = variante_split[0:8] + [variante_split[8]+':SB:SQD']
                format_sample = variante_split[header_chrom.index(sample)]  +':'+ sSB + ':' + str(sQD)
            else:
                variante_common = variante_split[0:9]
                format_sample = variante_split[header_chrom.index(sample)]
                
            sample_vcf.write('\t'.join(variante_common + [format_sample]) + '\n')
            
        sample_vcf.close()

def main():

    parser = argparse.ArgumentParser('Parse VCF output from Variant callers to output a variant_dataset.txt.  Output is to stdout.')
    parser.add_argument('-f', '--freebayes', help="Freebayes vcf output file name")
    parser.add_argument('-g', '--gatk', help="gatk vcf output file name")
    parser.add_argument('-v', '--varscan', help="Varscan vcf output file name")
    parser.add_argument('-l', '--listaFeatures', help="Lista di features da stampare",default=None)
    parser.add_argument('-s', '--split', help="Split vcf per samples", action='store_true')
    parser.add_argument('-F', '--feat_extraction', help="Enable features extraction", action='store_true')
    parser.add_argument('-a', '--amplicon',help="Amplicon design", action='store_true')
    parser.add_argument('-o', '--out_path',help="path di output")
    parser.add_argument('-G', '--gvcf_path',help="gvcf path")

    global opts 
    opts = parser.parse_args()
    callers = [opts.gatk,opts.varscan,opts.freebayes]
    samples = samples_name_extract(open(opts.freebayes,'r'))
    
    try:
        os.mkdir(opts.out_path)
    except:
        pass
    
    if opts.split:
        print 'Splitto le varianti per campione...'
        for vcf_dir in callers:
            split_vcf(vcf_dir,samples)
        print'Done'

    if opts.feat_extraction:
        varianti_total = dict()
        for dir_sample in os.listdir(opts.out_path):
            varianti = dict()
            vcf_path = opts.out_path +'/' + dir_sample
            if os.path.isdir(vcf_path):
                print "Analizzo le varianti da: " + vcf_path
                for vcf_name in os.listdir(vcf_path) :
                    print vcf_name
                    if 'Free' in vcf_name:
                        index = 0
                    elif 'GATK' in vcf_name:
                        index = 2
                    elif 'VarScan' in vcf_name:
                        index = 1
                    
                    in_file = open(vcf_path + '/' + vcf_name)
                    vcfreader = read(in_file,index,varianti)
                
                set_features(varianti)
                print_var(varianti,opts.out_path,dir_sample)
            for var in varianti.keys():
                #varianti_total[var] = var.split('\t')[0]+'\t'+var.split('\t')[1]+'\t.\t'+var.split('\t')[2]+'\t'+var.split('\t')[3]+'\t.\t.\t.\t.\t.'
                varianti_total[var] = ''
        print_vcf(varianti_total,opts.out_path)
main()