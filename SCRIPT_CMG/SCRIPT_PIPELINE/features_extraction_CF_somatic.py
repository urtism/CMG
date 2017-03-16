import re
import string
import argparse
import sys
import statistics
import scipy.stats as stats


''' //////////// CLASSI ////////////'''

class Caller():
    GT_n=''
    GT_t=''
    AO_t=''
    AO_n=''
    RO_t=''        
    RO_n=''
    AO_f=''
    AO_r=''
    DP_n=0
    DP_t=0
    QB_t=''
    QB_n=''
    Somatic=''
    Call=''
    AF_t=''
    AF_n=''
    FILTER=''

class Mutect(Caller):
    t_ad=''
    t_af=''
    t_QSS=''
    
    n_ad=''
    n_af=''
    n_QSS=''
    
class Vardict(Caller):
    STATUS=''
    
class Varscan(Caller):
    pass
    
class Features():

    GT_t_Varscan='.'
    GT_t_Vardict='.'
    GT_t_Mutect='.'
    
    GT_n_Varscan='.'
    GT_n_Vardict='.'
    GT_n_Mutect='.'
    
    DP_t_Varscan='.'
    DP_t_Vardict='.'
    DP_t_Mutect='.'

    DP_n_Varscan='.'
    DP_n_Vardict='.'
    DP_n_Mutect='.'
    
    DP=float(0)
    
    MBQT='.'
    
    QB_n_Mutect='.'
    QB_n_Varscan='.'
    QB_n_Vardict='.'
       
    AF_t_Mutect='0'
    AF_n_Mutect='.'
    AF_t_Varscan='0'
    AF_n_Varscan='.'
    AF_t_Vardict='0'
    AF_n_Vardict='.'
    
    Delta_Mutect='.'
    Delta_Varscan='.'
    Delta_Vardict='.'
    
    Delta_perc_Mutect='.'
    Delta_perc_Varscan='.'
    Delta_perc_Vardict='.'

    MBQT='.'
    MBQT_mediana='.'

    MBQN_media='.'
    MBQN_mediana='.'

    delta_media='.'
    delta_mediana='.'
    Delta_perc_media='.'
    Delta_perc_mediana='.'

    SomaticMutect='.'
    SomaticVarscan='.'
    SomaticVardict='.'
    
    CallMutect='0'
    CallVarscan='0'
    CallVardict='0'
    
    SS_Mutect='.'
    SS_Varscan='.'
    SS_Vardict='.'
    
    INFO_Mutect='.'
    INFO_Varscan='.'
    INFO_Vardict='.'

    AO_tum_media='.'
    RO_tum_media='.'
    AO_tum_mediana='.'
    RO_tum_mediana='.'
    
    DP_tum_media='.'
    DP_tum_mediana='.'

    AF_tum_media='.'
    AF_tum_mediana='.'

    AO_norm_media='.'
    RO_norm_media='.'
    AO_norm_mediana='.'
    RO_norm_mediana='.'

    DP_norm_media='.'
    DP_norm_mediana='.'

    AF_norm_media='.'
    AF_norm_mediana='.'


   
    FILTER_Mutect='.'     
    STATUS_Vardict='.'






''' //////////// FUNZIONI ////////////'''


def get_info_Mutect(chrom,pos,ref,alt,filter,info,format,tumor,normal,Mutect):
    '''estrae le informazioni dal vcf di Mutect'''
    Mutect.GT_t=tumor[format.index('GT')]
    Mutect.GT_n=normal[format.index('GT')]
    Mutect.AO_t=float((tumor[format.index('AD')]).split(',')[1])
    Mutect.RO_t=float((tumor[format.index('AD')]).split(',')[0])
    Mutect.DP_t=Mutect.AO_t+Mutect.RO_t
    Mutect.AF_t=Mutect.AO_t/Mutect.DP_t
    Mutect.DP_n=float((normal[format.index('AD')]).split(',')[0]) + float((normal[format.index('AD')]).split(',')[1])
    

    if normal is not 'null' and Mutect.DP_n is not '0':
        try:
            Mutect.AO_n=float((normal[format.index('AD')]).split(',')[1])
        except:
            Mutect.AO_n=float(0)
        try:
            Mutect.RO_n=float((normal[format.index('AD')]).split(',')[0])
        except:
            Mutect.RO_n=float(0)
        try:
            Mutect.AF_n=Mutect.AO_n/Mutect.DP_n
        except:
            Mutect.AF_n='.'
    else:
        Mutect.DP_n=0
        Mutect.AF_n='.'
        Mutect.AO_n='.'
        Mutect.RO_n='.'

    
  
    try:
        Mutect.QB_t=(float((tumor[format.index('QSS')]).split(',')[0]) + float((tumor[format.index('QSS')]).split(',')[1]))/Mutect.DP_t
    except:
        Mutect.QB_t='.'
    try:
        Mutect.QB_n=(float((normal[format.index('QSS')]).split(',')[0]) + float((normal[format.index('QSS')]).split(',')[1]))/Mutect.DP_n
    except:
        Mutect.QB_n='.'
    
    Mutect.Call=1    
    Mutect.FILTER=filter

    
def get_info_varscan(chrom,pos,ref,alt,filter,info,format,tumor,normal,varscan):
    '''estrae le informazioni dal vcf di varscan'''
    varscan.GT_t=tumor[format.index('GT')]
    varscan.GT_n=normal[format.index('GT')]
    varscan.AO_t=float(tumor[format.index('AD')])
    varscan.RO_t=float(tumor[format.index('RD')])
    varscan.DP=float(tumor[format.index('DP')])
    varscan.DP_t=float(tumor[format.index('DP')])

    if normal is not 'null' and normal[format.index('DP')] is not '0':
        varscan.AO_n=float(normal[format.index('AD')])
        varscan.RO_n=float(normal[format.index('RD')])
        varscan.DP_n=float(normal[format.index('DP')])
        varscan.AF_n=float(varscan.AO_n/varscan.DP_n)
        
    else:
        varscan.AO_n='.'
        varscan.RO_n='.'
        varscan.AF_n='.'
        varscan.DP_n=0

    for ind in info:
        if ind.startswith('SOMATIC'):
            varscan.Somatic=1
            break
        else:
            varscan.Somatic=0
                
   
    varscan.Call=1
    varscan.QB_t='.'
    varscan.QB_n='.'
    
    try:
        varscan.AF_t=float(varscan.AO_t/(varscan.DP_t))
    except:
        varscan.AF_t=float(0)    


def get_info_vardict(chrom,pos,ref,alt,filter,info,format,tumor,normal,vardict):
    '''estrae le informazioni dal vcf di vardict'''
    vardict.GT_t=tumor[format.index('GT')]
    vardict.GT_n=normal[format.index('GT')]
    vardict.AO_t=float(tumor[format.index('VD')])
    vardict.RO_t=float(tumor[format.index('RD')].split(',')[0])+float(tumor[format.index('RD')].split(',')[1])
    vardict.DP_t=float(tumor[format.index('DP')])
    vardict.DP=float(tumor[format.index('DP')])

    if normal is not 'null' and normal[format.index('DP')] is not '0' :
        
        vardict.AO_n=float(normal[format.index('VD')])
        vardict.RO_n=float(normal[format.index('RD')].split(',')[0])+float(normal[format.index('RD')].split(',')[1])
        vardict.DP_n=float(normal[format.index('DP')])
        vardict.AF_n=float(vardict.AO_n/(vardict.DP_n))
    else:
        vardict.AO_n='.'
        vardict.RO_n='.'
        vardict.DP_n=0
        vardict.AF_n='.'
    for ind in info:
        if ind.startswith("STATUS"):    
            vardict.STATUS=ind.split('=')[1]
                
    if "Somatic" in vardict.STATUS:
        vardict.Somatic=1
    else:
        vardict.Somatic=0    

    vardict.QB=float(tumor[format.index('QUAL')])
    vardict.QB_n=float(normal[format.index('QUAL')])

    vardict.Call=1

    try:
        vardict.AF_t=float(vardict.AO_t/(vardict.AO_t + vardict.RO_t))
    except:
        vardict.AF_t=float(0)

def set_features_snp(dictionary):
    '''setta i valori delle features in base alle info estratte dai vcf'''
    for variante in dictionary.keys():
        features=Features()
        varc_array=dictionary.get(variante)
        
        vett_MBQ_t=[]
        vett_MBQ_n=[]
        vett_DP_t=[]
        vett_DP_n=[]
        vett_AO_tum=[]
        vett_AO_norm=[]
        vett_RO_tum=[]
        vett_RO_norm=[]


        index=0

        for varcall in varc_array:
            if varcall is not "":
                vett_MBQ_t=vett_MBQ_t+[varcall.QB_t]
                vett_MBQ_n=vett_MBQ_n+[varcall.QB_n]
                vett_DP_t=vett_DP_t+[varcall.DP_t]
                vett_DP_n=vett_DP_n+[varcall.DP_n]
                vett_AO_tum=vett_AO_tum+[varcall.AO_t]
                vett_RO_tum=vett_RO_tum+[varcall.RO_t]
                vett_AO_norm=vett_AO_norm+[varcall.AO_n]
                vett_RO_norm=vett_RO_norm+[varcall.RO_n]
                
                if index == 0:

                    features.GT_t_Mutect=varc_array[0].GT_t
                    features.GT_n_Mutect=varc_array[0].GT_n
                    features.QB_t_Mutect=varc_array[0].QB_t
                    features.QB_n_Mutect=varc_array[0].QB_n
                    features.AF_t_Mutect=varc_array[0].AF_t
                    features.AF_n_Mutect=varc_array[0].AF_n
                    
                    if features.AF_n_Mutect is not ".":
                        features.Delta_Mutect=float(features.AF_t_Mutect) - float(features.AF_n_Mutect)

                    try:
                        features.Delta_perc_Mutect=float(features.AF_t_Mutect) / float(features.AF_n_Mutect)
                    except:
                        if features.AF_n_Mutect == 0:
                            features.Delta_perc_Mutect='100'
                        elif features.AF_n_Mutect is '.':
                            features.Delta_perc_Mutect='.'
                        
                    features.CallMutect=varc_array[0].Call
                    features.FILTER_Mutect=varc_array[0].FILTER


                elif index == 1:
                    
                    features.GT_t_Varscan=varc_array[1].GT_t
                    features.GT_n_Varscan=varc_array[1].GT_n
                    features.AF_t_Varscan=varc_array[1].AF_t
                    features.AF_n_Varscan=varc_array[1].AF_n

                    if features.AF_n_Varscan is not '.':
                            features.Delta_Varscan=float(features.AF_t_Varscan) - float(features.AF_n_Varscan)
                                    
                    try:
                        features.Delta_perc_Varscan=float(features.AF_t_Varscan) / float(features.AF_n_Varscan)
                    except:
                        if features.AF_n_Varscan == 0.0:
                            features.Delta_perc_Varscan='100'
                        elif features.AF_n_Varscan is '.':
                            features.Delta_perc_Varscan='.'

                    features.SomaticVarscan=varc_array[1].Somatic
                    features.CallVarscan=varc_array[1].Call
                                            
                elif index == 2:

                    features.GT_t_Vardict=varc_array[2].GT_t
                    features.GT_n_Vardict=varc_array[2].GT_n
                    features.QB_t_Vardict=varc_array[2].QB_t
                    features.QB_n_Vardict=varc_array[2].QB_n
                    features.AF_t_Vardict=varc_array[2].AF_t
                    features.AF_n_Vardict=varc_array[2].AF_n

                    if features.AF_n_Vardict is not ".":
                        features.Delta_Vardict=float(features.AF_t_Vardict) - float(features.AF_n_Vardict)

                    try:
                        features.Delta_perc_Vardict=float(features.AF_t_Vardict) / float(features.AF_n_Vardict)
                    except:
                        if features.AF_n_Vardict == 0.0:
                            features.Delta_perc_Vardict='100'
                        elif features.AF_n_Vardict is '.':
                            features.Delta_perc_Vardict='.'

                    features.SomaticVardict=varc_array[2].Somatic
                    features.CallVardict=varc_array[2].Call
                    features.STATUS_Vardict=varc_array[2].STATUS


            index = index + 1    

        vett_delta=[features.Delta_Mutect,features.Delta_Varscan,features.Delta_Vardict]
        vett_delta_perc=[features.Delta_perc_Mutect,features.Delta_perc_Varscan,features.Delta_perc_Vardict]
        vett_AF_tum_media=[features.AF_t_Mutect,features.AF_t_Varscan,features.AF_t_Vardict]
        vett_AF_norm_media=[features.AF_n_Mutect,features.AF_n_Varscan,features.AF_n_Vardict]
        delta_m=0
        delta_m_perc=0
        
        AF_med=0
    
        AO_tum_media='.'
        RO_tum_media='.'
        AO_tum_mediana='.'
        RO_tum_mediana='.'

        AO_norm_media='.'
        RO_norm_media='.'
        AO_norm_mediana='.'
        RO_norm_mediana='.'

        i=0
        nDP=0
        v=[]
        for dp in vett_DP_t:
            if dp and dp is not '':
                nDP= float(nDP)+float(dp)
                v=v+[float(dp)]
                i=i+1
        try:
            features.DP_tum_media=nDP/i
        except:
            features.DP_tum_media='.'
        try:
            features.DP_tum_mediana= int(statistics.median(v))
        except:
            features.DP_tum_mediana='.'

        i=0
        nDP=0
        v=[]
        for dp in vett_DP_n:
            if dp and dp is not '':
                nDP= float(nDP)+float(dp)
                v=v+[float(dp)]
                i=i+1
        try:
            features.DP_norm_media=nDP/i
        except:
            features.DP_norm_media='.'
        try:
            features.DP_norm_mediana= int(statistics.median(v))
        except:
            features.DP_norm_mediana='.'
        
        i=0
        v=[]
        nAO=0
        for ao in vett_AO_tum:
            if ao is not '':
                nAO= int(nAO)+int(ao)
                v=v+[int(ao)]
                i=i+1
        try:
            features.AO_tum_media=nAO/i
        except:
            features.AO_tum_media='.'
        try:
            features.AO_tum_mediana= int(statistics.median(v))
        except:
            features.AO_tum_mediana='.'

        i=0
        v=[]
        nAO=0        
        for ao in vett_AO_norm:
            if ao is not '.':
                nAO= int(nAO)+int(ao)
                v=v+[float(ao)]
                i=i+1

        try:
            features.AO_norm_media=nAO/i
            #print 'media',features.AO_norm_media
        except:
            features.AO_norm_media='.'
        try:
            features.AO_norm_mediana= int(statistics.median(v))
        except:
            features.AO_norm_mediana='.'
        i=0
        v=[]
        nRO=0
        for ro in vett_RO_tum:
            if ro is not '.':
                nRO= int(nRO)+int(ro)
                v=v+[int(ro)]
                i=i+1
        try:
            features.RO_tum_media=nRO/i
        except:
            features.RO_tum_media='.'
        try:
            features.RO_tum_mediana= statistics.median(v)
        except:
            features.RO_tum_mediana='.'

        i=0
        v=[]
        nRO=0
        for ro in vett_RO_norm:
            if ro is not '.':
                nRO= int(nRO)+int(ro)
                v=v+[int(ro)]
                i=i+1
        try:
            features.RO_norm_media=nRO/i
        except:
            features.RO_norm_media='.'
        try:
            features.RO_norm_mediana= statistics.median(v)
        except:
            features.RO_norm_mediana='.'

        i=0
        nMBQT=0
        v=[]
        for bq in vett_MBQ_t:
            #print variante,vett_MBQ,bq
            if bq and bq is not '.':
                nMBQT= float(nMBQT)+float(bq)
                if bq is not '0':
                    v=v+[float(bq)]
                i=i+1
        try:
            features.MBQT=nMBQT/i
        except:
            features.MBQT='.'
        try:
            features.MBQT_mediana=round(statistics.median(v))
        except:
            features.MBQT_mediana='.'    

        i=0
        nMBQN=0
        v=[]
        for bq in vett_MBQ_n:
            #print variante,vett_MBQ,bq
            if bq is not '.':
                try:
                    nMBQN= float(nMBQN)+float(bq)
                except:
                    print variante,vett_MBQ_n,bq
                if bq is not '0':
                    v=v+[float(bq)]
                i=i+1
        try:
            features.MBQN_media=nMBQN/i
        except:
            features.MBQN_media='.'
        try:
            features.MBQN_mediana=round(statistics.median(v))
        except:
            features.MBQN_mediana='.'

        i=0
        v=[]
        for delta in vett_delta:
            if delta is not '.':
                delta_m = float(delta_m)+ float(delta)
                v=v+[float(delta)]
                i=i+1
        try:
            features.delta_media = delta_m/i
        except:
            features.delta_media='.'
        try:
            features.delta_mediana=statistics.median(v)
        except:
            features.delta_mediana='.'
        
        i=0
        v=[]
        for delta_perc in vett_delta_perc:
            if delta_perc and delta_perc is not '.':
                delta_m_perc = float(delta_m_perc) + float(delta_perc)
                if delta_perc is not '0':
                    v=v+[float(delta_perc)]
                i=i+1
        try:
            features.Delta_perc_media= delta_m_perc/i
        except:
            features.Delta_perc_media='.'
        try:
            features.Delta_perc_mediana= statistics.median(v)
        except:
            features.Delta_perc_mediana='.'

        i=0
        v=[]
        for af in vett_AF_tum_media:
            if af and af is not '0':
                AF_med=float(AF_med) + float(af)
                v=v+[float(af)]
                i=i+1
        try:
            features.AF_tum_media= AF_med/i
        except:
            features.AF_tum_media='.'
        try:
            features.AF_tum_mediana= statistics.median(v)
        except:
            features.AF_tum_mediana='.'

        i=0
        AF_med=0
        v=[]
        for af in vett_AF_norm_media:
            if af is not '.':
                AF_med=float(AF_med) + float(af)
                v=v+[float(af)]
                i=i+1
        try:
            features.AF_norm_media= AF_med/i
        except:
            features.AF_norm_media='.'
        try:
            features.AF_norm_mediana= statistics.median(v)
        except:
            features.AF_norm_mediana='.'


        #print vett_delta,"\tdeltamedia",features.delta_media,"\n",vett_delta_perc,"\tdelta_perc_media",features.Delta_perc_media,"\n",vett_AF_media,"\taf_media",features.AF_media,"\n",vett_STRB_media,"\tstrb",features.STRBIAS_media
        dictionary[variante]= varc_array + [features]

def switch_snp(dictionary,ID,index,chrom,pos,ref,alt,filter,info,format,tumor,normal):
    '''tramite index richiama la funzione di estrazione delle informazioni del variant caller associato all'indice'''
    if dictionary.has_key(ID):
        vettore=dictionary[ID]
    else:
        vettore=['','','']

    if index==0:
        # print 'Mutect'
        mutect=Mutect()
        get_info_Mutect(chrom,pos,ref,alt,filter,info,format,tumor,normal,mutect)
        if Mutect.AF_t != 0.0 or Mutect.AF_t != '.': 
            vettore[0]=mutect
    elif index==3:
        # print 'vardict'
        vardict=Vardict()
        get_info_vardict(chrom,pos,ref,alt,filter,info,format,tumor,normal,vardict)
        if vardict.AF_t != 0.0 or vardict.AF_t != '.': 
            vettore[2]=vardict
    elif index==1 or index==2:
        # print 'varscan'
        varscan=Varscan()
        get_info_varscan(chrom,pos,ref,alt,filter,info,format,tumor,normal,varscan)
        if varscan.AF_t != 0.0 or varscan.AF_t != '.': 
            vettore[1]=varscan
    dictionary[ID]=vettore

def read(iterable,index,dictionary):
    '''legge il vcf e splitta le varie sezioni'''
    normal_tumor=0
    chrom=''
    alt=''
    pos=''
    ref=''
    filter=''
    format=''
    info=''
    for line in iterable:
        line.rstrip()
        #print line
        if line.startswith('#CHROM'):
            parts = line.split("\t")
            if 'NORMAL' in parts[9] or parts[9] == opts.normal:
                #print "prima normal",parts[9]
                normal_tumor=1
            else: 
                #print "prima tumor"
                normal_tumor=0
        elif line.startswith('#'):
            continue
        else:
            #print line
            ind=0
            parts = line.split("\t")
            #print parts[1]
            chrom=parts[0]
            pos=parts[1]
            ref=parts[3]
            alt=parts[4]
            filter=parts[6]
            ID='\t'.join([chrom,pos,ref,alt])
            INFO=parts[7]
            info=INFO.split(";") 
            FORMAT=parts[8]
            format=FORMAT.split(":")
            tumor=''
            normal=''
            if normal_tumor==1:
                NORMAL=parts[9]
                TUMOR=parts[10].rstrip()
                tumor=TUMOR.split(":")
                normal=NORMAL.split(":")
            else:
                NORMAL=parts[10].rstrip()
                TUMOR=parts[9]
                tumor=TUMOR.split(":")
                normal=NORMAL.split(":")
            
            if NORMAL.startswith('.:.:.:') or len(normal)==1:
                normal='null'
            if TUMOR.startswith('.:.:.:') or len(tumor)==1 :
                continue
            else:
                if "VD" in format:
                    if tumor[format.index('VD')] is '0':
                        continue
                    else:
                        switch_snp(dictionary,ID,index,chrom,pos,ref,alt,filter,info,format,tumor,normal)
                else:    
                    switch_snp(dictionary,ID,index,chrom,pos,ref,alt,filter,info,format,tumor,normal)

def max_delta_perc(dictionary):
    ''' calcolo il max delta percentuale'''
    buffer_delta_perc='0'
    max_delta_perc='0'
    
    for variante in dictionary.keys():
        features = dictionary.get(variante)[5]
        if features.Delta_perc_Mutect == '.':
            mutect=0
        else:
            mutect=features.Delta_perc_Mutect
        if features.Delta_perc_Vardict == '.':
            vard=0
        else:
            vard=features.Delta_perc_Vardict
        if features.Delta_perc_Varscan == '.':
            varsc=0
        else:
            varsc=features.Delta_perc_Varscan
        if features.Delta_perc_Radia == '.':
            rad=0
        else:
            rad=features.Delta_perc_Radia
        if features.Delta_perc_Somaticsniper == '.':
            som=0
        else:
            som=features.Delta_perc_Somaticsniper
        buffer_delta_perc= max([float(mutect),float(varsc),float(vard),float(rad),float(som)])
        if float(max_delta_perc) < float(buffer_delta_perc):
            max_delta_perc=buffer_delta_perc
    #print max_delta_perc        

def control(dictionary):
    ''' esegue un controllo sulle varianti, se non hanno variant caller che le chiama vengono eliminate'''
    for variante in dictionary.keys():
        if dictionary[variante][:3] == ['','','']:
            #print "sto cancellando:",variante
            del dictionary[variante]
                
def print_var_snp_complete(dictionary):
    varianti_tsv=open(opts.out+ '.features.tsv','w')
    varianti_tsv.write('\t'.join(["SAMPLE_NORMAL_ID","SAMPLE_TUMOR_ID","CHROM","POS","REF","ALT","CallMutect","CallVarscan","CallVardict",
            "SomaticMutect","SomaticVarscan","SomaticVardict",
            "GT_t_Mutect","GT_t_Varscan","GT_t_Vardict",
            "GT_n_Mutect","GT_n_Varscan","GT_n_Vardict","DP_median","BQ_t_Mutect","BQ_t_Vardict","MBQT_media","MBQT_mediana",
            "AF_Mutect","AF_Varscan","AF_Vardict","AF_media","AF_median",
            "Delta_Mutect","Delta_Varscan","Delta_Vardict","Delta_media","Delta_median",
            "Delta_perc_Mutect","Delta_perc_Varscan","Delta_perc_Vardict","Delta_perc_media","Delta_perc_median",
            "DP_n/t_Mutect","DP_n/t_varscan","DP_n/t_vardict","DP_n/t_media","DP_n/t_median",
            "T_lod_Mutect","N_lod_Mutect","Delta_lod_Mutect","T/N_lod_Mutect",
            "FOXOG_Mutect","GQ_Mutect","PGT_Mutect","PID_Mutect","PL_Mutect","HCNT_Mutect",
            "MAX_ED_Mutect","MIN_ED_Mutect","PON_Mutect","RPA_Mutect","RU_Mutect","STR_Mutect",
            "STRBIAS_Mutect","STRBIAS_Varscan","STRBIAS_Vardict","STRBIAS_medio","STRBIAS_median",
            "ODDRATIO_Vardict","SBF_Vardict","SHIFT3_Vardict","MSI_Vardict","MSILEN_Vardict","SOR_Vardict",
            "LSEQ_Vardict","RSEQ_Vardict","STATUS_Vardict","PMEAN_Vardict","PSTD_Vardict","QSTD_Vardict",
            "MQ_Vardict","SN_Vardict","HIAF_Vardict","NM_Vardict",
            "LOH_Varscan","LOH_Vardict"]) +'\n')


    for variante in dictionary.keys():
        features = dictionary.get(variante)[-1]

        varianti_tsv.write('\t'.join([opts.normal,opts.tumor,variante,str(features.CallMutect),str(features.CallVarscan),str(features.CallVardict),
            str(features.SomaticMutect),str(features.SomaticVarscan),str(features.SomaticVardict),
            str(features.GT_t_Mutect),str(features.GT_t_Varscan),str(features.GT_t_Vardict),
            str(features.GT_n_Mutect),str(features.GT_n_Varscan),str(features.GT_n_Vardict),
            str(features.DP_median),str(features.QB_t_Mutect),str(features.QB_t_Vardict),str(features.MBQT),str(features.MBQT_median),
            str(features.AF_t_Mutect),str(features.AF_t_Varscan),str(features.AF_t_Vardict),str(features.AF_media),str(features.AF_median),
            str(features.Delta_Mutect),str(features.Delta_Varscan),str(features.Delta_Vardict),str(features.delta_media),str(features.delta_median),
            str(features.Delta_perc_Mutect),str(features.Delta_perc_Varscan),str(features.Delta_perc_Vardict),str(features.Delta_perc_media)]).rstrip() +'\n')

def print_var_snp_reduced(dictionary):
    
    varianti_tsv=open(opts.out+ '.features.tsv','w')
    
    varianti_tsv.write('\t'.join(["CHROM","POS","ID","REF","ALT","CallMutect","CallVarscan","CallVardict",
            "SomaticVarscan","SomaticVardict",
            "FILTER_Mutect","STATUS_Vardict",
            "GT_t_Mutect","GT_t_Varscan","GT_t_Vardict",
            "GT_n_Mutect","GT_n_Varscan","GT_n_Vardict",
            "DP_TUM",
            "AF_TUM_Mutect","AF_TUM_Varscan","AF_TUM_Vardict","AF_TUM_media",
            "AO_TUM","RO_TUM","BQ_TUM_media",
            "DP_NORM",
            "AF_NORM_Mutect","AF_NORM_Varscan","AF_NORM_Vardict","AF_NORM_media",
            "AO_NORM","RO_NORM","BQ_NORM_media",
            "Delta_mediana"]) + '\n')

    for variante in dictionary.keys():
        features = dictionary.get(variante)[-1]

        varianti_tsv.write('\t'.join([variante.split('\t')[0],variante.split('\t')[1],opts.tumor,variante.split('\t')[2],variante.split('\t')[3],
            str(features.CallMutect),str(features.CallVarscan),str(features.CallVardict),
            str(features.SomaticVarscan),str(features.SomaticVardict),
            str(features.FILTER_Mutect),str(features.STATUS_Vardict),
            str(features.GT_t_Mutect),str(features.GT_t_Varscan),str(features.GT_t_Vardict),
            str(features.GT_n_Mutect),str(features.GT_n_Varscan),str(features.GT_n_Vardict),
            str(features.DP_tum_median),
            str(features.AF_t_Mutect),str(features.AF_t_Varscan),str(features.AF_t_Vardict),str(features.AF_tum_media),
            str(features.AO_tum_media),str(features.RO_tum_media),str(features.MBQT_median),
            str(features.DP_norm_median),
            str(features.AF_n_Mutect),str(features.AF_n_Varscan),str(features.AF_n_Vardict),str(features.AF_norm_media),
            str(features.AO_norm_media),str(features.RO_norm_media),str(features.MBQN_median),
            str(features.delta_median)]).rstrip() + '\n')
    
    varianti_tsv.close()

def print_vcf(varianti):
    varianti_vcf=open(opts.out+ '.vcf','w')
    varianti_vcf.write('##fileformat=VCFv4.2\n'+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLES\n')
    for variante in varianti.keys():
        var_vcf=variante.split('\t')[0]+'\t'+variante.split('\t')[1]+'\t.\t'+variante.split('\t')[2]+'\t'+variante.split('\t')[3]+'\t.\t.\t.\t.\t.'
        varianti_vcf.write(var_vcf+ '\n')
    varianti_vcf.close()

def print_var(dictionary):

    lista_features=open(opts.listaFeatures,'r')
    varianti_tsv=open(opts.out+ '.features.tsv','w')
    
    header=[]
    features_variante=[]
    
    for line in lista_features:
        line = line.rstrip()
        if line.startswith('#'):
            continue
        else:
            header=header+[line]
            features_variante=features_variante+['features.'+line]

    dataset_varianti.write('CHROM\tPOS\tID\tREF\tALT\t' + '\t'.join(header)+ '\n')
    for variante in dictionary.keys():
        features = dictionary.get(variante)[-1]
        features_variante_eval=[]
        for feat in features_variante:
             feat_eval=str(eval(feat))
             features_variante_eval=features_variante_eval + [feat_eval]
        var=variante.split('\t')[0]+'\t'+variante.split('\t')[1]+'\t' +sample_name +'\t'+variante.split('\t')[2]+'\t'+variante.split('\t')[3]+ '\t' + '\t'.join(features_variante_eval)
        if features.AF_media != '.':
                dataset_varianti.write(var+ '\n')
                break
    dataset_varianti.close()
    lista_features.close()
        
def main():


    
    parser = argparse.ArgumentParser('Parse VCF output from Variant callers to output a variant_dataset.txt.  Output is to stdout.')
    parser.add_argument('-m', '--mutect', help="Mutect vcf output file name")
    parser.add_argument('-d', '--vardict', help="Vardict vcf output file name")
    parser.add_argument('-v', '--varscan_snp', help="Varscan vcf output file name")
    parser.add_argument('-i', '--varscan_indel', help="Varscan vcf output file name")
    parser.add_argument('-n','--normal',help="Name of normal sample")
    parser.add_argument('-t','--tumor',help="Name of tumor sample")
    parser.add_argument('-a','--amplicon',help="Amplicon design", action='store_true')
    parser.add_argument('-c','--complete',help="Print complete info", action='store_true')
    parser.add_argument('-o', '--out', help="file name in output. It returns file_name.features.tsv and file_name.vcf ")
    
    global opts 
    opts = parser.parse_args()
    
    callers = [opts.mutect,opts.varscan_snp,opts.varscan_indel,opts.vardict]
    varianti = dict() 
    index=0;
    for vcf in callers:
        #print vcf
        in_file = open(vcf,'r')
        vcfreader = read(in_file,index,varianti)
        index = index + 1
    set_features_snp(varianti)
    control(varianti)
    if opts.complete:
        print_var_snp_complete(varianti)
    else:
        print_var_snp_reduced(varianti)
    print_vcf(varianti)
main()
