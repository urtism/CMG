import os
import tarfile
import argparse
import pandas as pd
import hashlib
import seaborn
import re
import numpy as np
import matplotlib.pyplot as plt
from six.moves import urllib
from joblib import dump, load

from pandas.api.types import is_object_dtype
from pandas.plotting import scatter_matrix

from ggplot import *

from functools import partial

from sklearn.model_selection import StratifiedShuffleSplit, train_test_split, StratifiedKFold, cross_val_predict, cross_val_score, cross_validate, RandomizedSearchCV, GridSearchCV
from sklearn.feature_selection import mutual_info_classif, VarianceThreshold, RFE, RFECV
from sklearn.preprocessing import LabelBinarizer, OneHotEncoder, Imputer, StandardScaler, MaxAbsScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.base import clone, BaseEstimator, TransformerMixin
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, auc, precision_recall_curve, roc_curve, make_scorer, recall_score, precision_score, balanced_accuracy_score, f1_score

class CombinedAttributesAdder(BaseEstimator, TransformerMixin):
    def __init__(self, iQual = False, iAMQ = False, iAAS = False, iAXS = False, iAXS0 = False, iAMQ0 = False, iACR = False, iGC = False, ID = True):  #   no  *args   or  **kargs
        self.iQual = iQual
        self.iAMQ = iAMQ
        self.iAAS = iAAS
        self.iAXS = iAXS
        self.iAXS0 = iAXS0
        self.iAMQ0 = iAMQ0
        self.iACR = iACR
        self.iGC = iGC
        self.ID = ID

    def fit(self, DS, y=None):
        return self        #   nothing else    to  do
    def transform(self, DS, y=None):

        X = DS.values
        L = DS.columns.tolist()            

        if self.ID:
            try:
                ID = (X[:, L.index('SAMPLE_ID') ]+'-'+X[:,  L.index('VAR_ID')] )
                X = np.c_[X, ID]
                L += ['ID']
            except:
                ID = DS.index
                X = np.c_[X, ID]
                L += ['ID']

        if self.iGC:

            iGC = (X[:, L.index('IEVA-iGC') ]/100.0)
            X = np.c_[X, iGC]
            L += ['IEVA-iGCfraz']

        if self.iQual:
            try:
                iQual = (X[:, L.index('IEVA-iQual REF') ]-X[:,  L.index('IEVA-iQual ALT')] )/ (X[:, L.index('IEVA-iQual REF')] + X[:, L.index('IEVA-iQual ALT')])
            except:
                iQual = np.empty(len(X))
            X = np.c_[X, iQual]
            L += ['IEVA-iQual']

        if self.iAMQ:
            try:
                iAMQ = (X[:, L.index('IEVA-iAMMQ REF') ]-X[:,  L.index('IEVA-iAMMQ ALT')] )/ (X[:, L.index('IEVA-iAMMQ REF')] + X[:, L.index('IEVA-iAMMQ ALT')])
            except:
                iAMQ = np.empty(len(X))
            X = np.c_[X, iAMQ]
            L += ['IEVA-iAMQ']
        if  self.iAAS:
            try:
                iAAS = (X[:, L.index('IEVA-iAAS REF') ]-X[:,  L.index('IEVA-iAAS ALT')] )/ (X[:, L.index('IEVA-iAAS REF')] + X[:, L.index('IEVA-iAAS ALT')])
            except:
                iAAS = np.empty(len(X))
            X = np.c_[X, iAAS]
            L += ['IEVA-iAAS']

        if  self.iAXS:
            try:
                iAXS = (X[:, L.index('IEVA-iAXS REF') ]-X[:,  L.index('IEVA-iAXS ALT')] )/ (X[:, L.index('IEVA-iAXS REF')] + X[:, L.index('IEVA-iAXS ALT')])
            except:
                iAXS = np.empty(len(X))
            X = np.c_[X, iAXS]
            L += ['IEVA-iAXS']

        if  self.iAXS0:
            try:
                iAXS0 = (X[:, L.index('IEVA-iAXS0 REF') ]-X[:,  L.index('IEVA-iAXS0 ALT')] )/ (X[:, L.index('IEVA-iAXS0 REF')] + X[:, L.index('IEVA-iAXS0 ALT')])
                  
            except:
                iAXS0 = np.empty(len(X))
            X = np.c_[X, iAXS0]
            L += ['IEVA-iAXS0']

        if  self.iAMQ0:
            try:
                iAMQ0 = (X[:, L.index('IEVA-iAMQ0 REF') ]-X[:,  L.index('IEVA-iAMQ0 ALT')] )/ (X[:, L.index('IEVA-iAMQ0 REF')] + X[:, L.index('IEVA-iAMQ0 ALT')])
            except:
                iAMQ0 = np.empty(len(X))
            X = np.c_[X, iAMQ0]
            L += ['IEVA-iAMQ0']

        if  self.iACR:
            try:
                iACR = (X[:, L.index('IEVA-iACR REF') ]-X[:,  L.index('IEVA-iACR ALT')] )/ (X[:, L.index('IEVA-iACR REF')] + X[:, L.index('IEVA-iACR ALT')])
            except:
                iACR = np.empty(len(X))
            X = np.c_[X, iACR]
            L += ['IEVA-iACR']
        return X,L

def Split_datasets(ds,size,outpath,stratify):

    if stratify is not None:
        split = StratifiedShuffleSplit(n_splits=1, test_size=size)
        for train_index, test_index  in  split.split(ds, ds[stratify]):
            train_set =   ds.iloc[train_index]
            test_set  =   ds.iloc[test_index]
    else:
        train_set, test_set = train_test_split(ds, test_size=size, random_state=41)

    try:
        test_set.set_index(['SAMPLE_ID','VAR_ID'], inplace =True)
        train_set.set_index(['SAMPLE_ID','VAR_ID'], inplace =True)
    except:
        pass

    test_set.to_csv(path_or_buf=outpath+'.TEST.csv',sep='\t')
    train_set.to_csv(path_or_buf=outpath+'.TRAIN.csv',sep='\t')
    return train_set,test_set

def print_stats(dataset,OUTPUT = False):

    if dataset is None:
        dataset_list= ['/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/INDEL/PANDAS/INTARGET/BK-ENR.SNV.Freebayes+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/INDEL/PANDAS/INTARGET/BK-ENR.SNV.GATK+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/INDEL/PANDAS/INTARGET/BK-ENR.SNV.Platypus.+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/INDEL/PANDAS/INTARGET/BK-ENR.SNV.Samtools+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/INDEL/PANDAS/INTARGET/BK-ENR.SNV.Scalpel+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/INDEL/PANDAS/INTARGET/BK-ENR.SNV.SNVer+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/INDEL/PANDAS/INTARGET/BK-ENR.SNV.VarScan+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/SNV/PANDAS/INTARGET/BK-ENR.SNV.Freebayes+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/SNV/PANDAS/INTARGET/BK-ENR.SNV.GATK.target+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/SNV/PANDAS/INTARGET/BK-ENR.SNV.Platypus+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/SNV/PANDAS/INTARGET/BK-ENR.SNV.Samtools+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/SNV/PANDAS/INTARGET/BK-ENR.SNV.SNVer+-10.tsv',
            '/home/jarvis/Scrivania/TEST/bam-ieva/ENRICHMENT/out/SNV/PANDAS/INTARGET/BK-ENR.SNV.VarScan+-10.tsv']
    else:
        dataset_list =[dataset]



    for DATASET_PATH in dataset_list:
        
        try:
            dataset = pd.read_csv(DATASET_PATH,sep='\t')
            out = open(DATASET_PATH+'.STATS','w')
        except:
            dataset = DATASET_PATH.copy()
            if OUTPUT:
                out = open(OUTPUT+'.STATS','w')
        #corr_matrix = dataset.corr()
        out.write('\t'.join(['PASS','FILTER'])+'\n')
        try:
            
            out.write('\t'.join([str(dataset['CLASS'].value_counts()[0]),str(dataset['CLASS'].value_counts()[1])])+'\n\n')
        except:
            out.write('\t'.join([str(dataset['CLASS'].value_counts()[0]),'0'])+'\n\n')


        out.write('\t'.join(['LABEL','MAX','MIN','MEAN','VALUES','MISSING VALUES','MISSING VALUE %'])+'\n')

        for f in dataset.columns.values.tolist():
            if is_object_dtype(dataset[f]):
                MEAN = '-'
            else:
                MEAN = str(round(dataset[f].mean(),3))

            MAX = str(dataset[f].max())
            MIN = str(dataset[f].min())
            VALUES = str(dataset[f].count())
            MVALUES = str(dataset[f].isnull().sum())
            try:
                MPERC = str(round(float(dataset[f].isnull().sum())/len(dataset[f]),3))
            except:
                MPERC = '0.0'
            out.write('\t'.join([f,MAX,MIN,MEAN,VALUES,MVALUES,MPERC])+'\n')

def print_full(x):
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', -1)
    print(x)
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')
    pd.reset_option('display.width')
    pd.reset_option('display.float_format')
    pd.reset_option('display.max_colwidth')


def drop_pseudonucleotidi(dataset):

    dataset = dataset.drop([
        'IEVA-iPNC AA',
        'IEVA-iPNC AC',
        'IEVA-iPNC AG',
        'IEVA-iPNC AT',
        'IEVA-iPNC CA',
        'IEVA-iPNC CC',
        'IEVA-iPNC CG',
        'IEVA-iPNC CT',
        'IEVA-iPNC GA',
        'IEVA-iPNC GC',
        'IEVA-iPNC GG',
        'IEVA-iPNC GT',
        'IEVA-iPNC TA',
        'IEVA-iPNC TC',
        'IEVA-iPNC TG',
        'IEVA-iPNC TT'],axis=1)

    return dataset

if __name__ == '__main__':

    parser = argparse.ArgumentParser('Questo script contiene le funzioni base per la manipolazione del dataset.\n')
    parser.add_argument('-D', '--dataset', default=None, help="Dataset da analizzare")
    parser.add_argument('-s', '--split', default=None, help="Split del dataset in TESTset e TRAINset, indicare la test_size in decimali (es. 0.3 se si vuole avere un 70/30)")
    parser.add_argument('-st', '--stratsplit', default=None, help="Label utilizzato per stratificare il dataset splittato")
    parser.add_argument('-m', '--model', default=None, help="Classificatore allenato da utilizzare per la classificazione di un nuovo dataset")
    parser.add_argument('-tr', '--train', default=None, help="train set")
    parser.add_argument('-ts', '--test', default=None, help="test set")
    parser.add_argument('-p', '--pipe', default=None, help=" list of action to do in data cleaning and features extraction")
    parser.add_argument('-dl', '--droplist',default=None, help="List of features to be dropped from the dataset")
    parser.add_argument('-O', '--outpath',default=None, help="Path di output")

    global opts
    opts = parser.parse_args()

    if opts.dataset is not None:

        try:
            fulldataset = pd.read_csv(opts.dataset,sep='\t')
            fulldataset.set_index(['ID'], inplace =True)
            dataset = fulldataset.drop('CLASS',axis=1)
            print_stats(opts.dataset)


        except:
            fulldataset = pd.read_csv(opts.dataset,sep='\t')
            attr_adder = CombinedAttributesAdder()
            dataset_extra_attribs, dataset_extra_columns = attr_adder.transform(fulldataset)
            fulldataset = pd.DataFrame(dataset_extra_attribs, columns=dataset_extra_columns).reindex()
            
            if 'VAR_ID' in fulldataset.columns or 'SAMPLE_ID' in fulldataset.columns:
                fulldataset = fulldataset.drop(['VAR_ID','SAMPLE_ID'], axis=1)

            fulldataset.set_index(['ID'], inplace =True)
            fulldataset.to_csv(path_or_buf=opts.outpath+'.ID.csv',sep='\t')
            dataset = fulldataset.drop(['CLASS'],axis=1)
            print 'ATTENZIONE:\nUtilizare come nuovo dataset il dataset: '+opts.outpath+'.ID.csv'
            exit()

    OUTPUT = opts.outpath

    try:
        pipe = opts.pipe.split(',')
    except:
        pipe = []

##1) Divido il dataset in train e test set
    if opts.split is not None:
        train_set,test_set = Split_datasets(fulldataset,float(opts.split),opts.outpath,stratify=opts.stratsplit)

#PROVA: invece di eliminare le features ref, ho generato della features miste con ref e alt con la formula ref-alt/ref+alt [-1,+1]
    if 'SELECT' in pipe:
        toselect = ['GATK_INFO-FS',
            'GATK_INFO-QD',
            'GATK_FORMAT-DP',
            'GATK_FORMAT-AF',
            'GATK_FORMAT-GQ',
            'GATK_FORMAT-GT',
            'GATK_FORMAT-PL-0/0',
            'GATK_FORMAT-PL-0/1',
            'GATK_FORMAT-PL-1/1',
            'FB_FORMAT-GT',
            'FB_FORMAT-AF',
            'FB_FORMAT-DP',
            'FB_FORMAT-GL-0/0',
            'FB_FORMAT-GL-0/1',
            'FB_FORMAT-GL-1/1',
            'FB_FORMAT-GQ',
            'FB_FORMAT-QA',
            'VS_FORMAT-GT',
            'VS_FORMAT-GQ',
            'VS_FORMAT-FREQ',
            'VS_FORMAT-PVAL',
            'VS_FORMAT-ABQ',
            'VS_FORMAT-DP',
            'PT_INFO-SbPval',
            'PT_INFO-MQ',
            'PT_INFO-QD',
            'PT_FORMAT-GT',
            'PT_FORMAT-GL-0/0',
            'PT_FORMAT-GL-0/1',
            'PT_FORMAT-GL-1/1',
            'PT_FORMAT-GQ',
            'PT_FORMAT-NR',
            'PT_FORMAT-AF',
            'ST_INFO-MQ',
            'ST_FORMAT-GQ',
            'ST_FORMAT-PL-0/0',
            'ST_FORMAT-PL-0/1',
            'ST_FORMAT-PL-1/1',
            'ST_FORMAT-GT',
            'SV_INFO-AF',
            'SV_INFO-DP',
            'SV_FORMAT-GT',
            'SV_FORMAT-PL-0/0',
            'SV_FORMAT-PL-0/1',
            'SV_FORMAT-PL-1/1']

        for c in dataset.columns:
            if c.startswith('IEVA'):
                toselect += [c]

        dataset1 = dataset[toselect]

        OUTPUT += '.SELECT'

        dataset_to_print=pd.concat([dataset1,fulldataset['CLASS']],axis=1)
        try:
            dataset_to_print.set_index(['ID'], inplace =True)
        except:
            pass
        dataset_to_print.to_csv(path_or_buf=OUTPUT+'.csv',sep='\t')
        print_stats(OUTPUT+'.csv')

        for f in toselect:
            try:
                if f.startswith('IEVA'):
                    continue
                else:
                    dataset = dataset.drop(f,axis=1)
            except:
                continue

        dataset_to_print=pd.concat([dataset,fulldataset['CLASS']],axis=1)
        try:
            dataset.set_index(['ID'], inplace =True)
            dataset_to_print.set_index(['ID'], inplace =True)
        except:
            pass
        
        dataset_to_print.to_csv(path_or_buf=OUTPUT+'.NOSELECT.csv',sep='\t')
        print_stats(OUTPUT+'.NOSELECT.csv')


    if 'NF' in pipe:

        attr_adder = CombinedAttributesAdder(iQual = True, iAMQ = True, iAAS = True, iAXS = True, iAXS0 = True, iAMQ0 = True, iACR = True, iGC = True )
        dataset_extra_attribs, dataset_extra_columns = attr_adder.transform(dataset)
        dataset = pd.DataFrame(dataset_extra_attribs, columns=dataset_extra_columns)

        gt_map = {'0/0': '0', '0/1': '1','1/0': '1','1/1': '2', '':'','Undef':''}
        #dataset['FORMAT-GT'] = dataset['FORMAT-GT'].map(gt_map)

        dataset['GATK_FORMAT-GT'] = dataset['GATK_FORMAT-GT'].map(gt_map)
        dataset['FB_FORMAT-GT'] = dataset['FB_FORMAT-GT'].map(gt_map)
        dataset['VS_FORMAT-GT'] = dataset['VS_FORMAT-GT'].map(gt_map)
        dataset['PT_FORMAT-GT'] = dataset['PT_FORMAT-GT'].map(gt_map)
        dataset['ST_FORMAT-GT'] = dataset['ST_FORMAT-GT'].map(gt_map)
        dataset['SV_FORMAT-GT'] = dataset['SV_FORMAT-GT'].map(gt_map)


        try:
            dataset.set_index(['ID'], inplace =True)
        except:
            pass
        try:
            dataset = dataset.drop(['IEVA-iGC'], axis=1)
        except:
            pass
        
        OUTPUT += '.NF'

        if 'VAR_ID' in dataset.columns or 'SAMPLE_ID' in dataset.columns:
            dataset = dataset.drop(['VAR_ID','SAMPLE_ID'], axis=1)

        dataset_to_print=pd.concat([dataset,fulldataset['CLASS']],axis=1)
        try:
            dataset_to_print.set_index(['ID'], inplace =True)
        except:
            pass
        dataset_to_print.to_csv(path_or_buf=OUTPUT+'.csv',sep='\t')
        print_stats(OUTPUT+'.csv')


# #2) Gestire i missing values:2.0)  metto da parte le features di tipo object; 2.1) elimino le features con piu del 30% di missing values; 2.2) elimino le features riferite alla REF e all'ALT
    

# # #2.3)sostituisco con 0 i missing values nella features IEVA-iSRL ( lunghezza degli short repeats)

    

    #2.4)utilizzo imputer di sklearn con la mediana
    if 'MV' in pipe:

        dataset['IEVA-iSRL'] = dataset['IEVA-iSRL'].fillna(0.0)

        dataset['GATK_FORMAT-GT'] = dataset['GATK_FORMAT-GT'].fillna(-1)
        dataset['FB_FORMAT-GT'] = dataset['FB_FORMAT-GT'].fillna(-1)
        dataset['VS_FORMAT-GT'] = dataset['VS_FORMAT-GT'].fillna(-1)
        dataset['PT_FORMAT-GT'] = dataset['PT_FORMAT-GT'].fillna(-1)
        dataset['ST_FORMAT-GT'] = dataset['ST_FORMAT-GT'].fillna(-1)
        dataset['SV_FORMAT-GT'] = dataset['SV_FORMAT-GT'].fillna(-1)

        dataset['GATK_FORMAT-DP'] = dataset['GATK_FORMAT-DP'].fillna(0.0)
        dataset['FB_FORMAT-DP'] = dataset['FB_FORMAT-DP'].fillna(0.0)
        dataset['VS_FORMAT-DP'] = dataset['VS_FORMAT-DP'].fillna(0.0)
        dataset['SV_INFO-DP'] = dataset['SV_INFO-DP'].fillna(0.0)

        dataset['GATK_FORMAT-AF'] = dataset['GATK_FORMAT-AF'].fillna(0.0)
        dataset['FB_FORMAT-AF'] = dataset['FB_FORMAT-AF'].fillna(0.0)
        dataset['VS_FORMAT-FREQ'] = dataset['VS_FORMAT-FREQ'].fillna(0.0)
        dataset['PT_FORMAT-AF'] = dataset['FB_FORMAT-AF'].fillna(0.0)
        dataset['SV_INFO-AF'] = dataset['SV_INFO-AF'].fillna(0.0)


        num_dataset = dataset.drop(dataset.select_dtypes(['object']),axis=1)

        imputer = Imputer(strategy="median")
        imputer.fit(num_dataset)
        I = imputer.transform(num_dataset)
        X = pd.DataFrame(I, columns=num_dataset.columns)
        X.set_index(dataset.index, inplace =True)

        #2.5)riaggiungo le features di tipo object che possono servire

        dataset = pd.concat([X,dataset.select_dtypes(['object'])], axis=1)


        OUTPUT += '.MV'
        dataset_to_print=pd.concat([dataset,fulldataset['CLASS']],axis=1)
        try:
            dataset.set_index(['ID'], inplace =True)
            dataset_to_print.set_index(['ID'], inplace =True)
        except:
            pass
        dataset_to_print.to_csv(path_or_buf=OUTPUT+'.csv',sep='\t')
        print_stats(OUTPUT+'.csv')


    if 'NORM' in pipe:
        refill=[]

        donotNORM = [
            'IEVA-iGCfraz',
            'IEVA-iMIU',
            'IEVA-iDUP',
            'IEVA-iSB',
            'IEVA-iCRT',
            'IEVA-iQRT',
            'IEVA-iMRT',
            'IEVA-iPRT',
            'IEVA-iQual',
            'IEVA-iAMQ',
            'IEVA-iAAS',
            'IEVA-iAXS',
            'IEVA-iAXS0',
            'IEVA-iAMQ0',
            'IEVA-iACR',
            'IEVA-iRM',
            'IEVA-SR']

        donotNORM += ['GATK_FORMAT-GT','FB_FORMAT-GT','VS_FORMAT-GT','PT_FORMAT-GT','ST_FORMAT-GT','SV_FORMAT-GT']
        donotNORM += ['GATK_FORMAT-AF','FB_FORMAT-AF','VS_FORMAT-FREQ','PT_FORMAT-AF','SV_INFO-AF']
        donotNORM += ['GATK_INFO-InbreedingCoeff','GATK_INFO-MLEAF',
            'FB_INFO-AB','FB_INFO-PAIRED','FB_INFO-PAIREDR',
            'VS_FORMAT-PVAL',
            'PT_INFO-FR','PT_INFO-SbPval','PT_INFO-BRF',
            'ST_INFO-RPB','ST_INFO-BQB','ST_INFO-MQSB','ST_INFO-MQ0F','ST_INFO-ICB','ST_INFO-HOB','SV_INFO-PV','SV_INFO-SP','SV_INFO-FS']

        donotNORM += ['IEVA-iPNC AA',
            'IEVA-iPNC AC',
            'IEVA-iPNC AG',
            'IEVA-iPNC AT',
            'IEVA-iPNC CA',
            'IEVA-iPNC CC',
            'IEVA-iPNC CG',
            'IEVA-iPNC CT',
            'IEVA-iPNC GA',
            'IEVA-iPNC GC',
            'IEVA-iPNC GG',
            'IEVA-iPNC GT',
            'IEVA-iPNC TA',
            'IEVA-iPNC TC',
            'IEVA-iPNC TG',
            'IEVA-iPNC TT']
        
        num_dataset = dataset.drop(dataset.select_dtypes(['object']),axis=1)

        for f in donotNORM:
            try:
                num_dataset = num_dataset.drop(f,axis=1)
                refill += [f]
            except:
                continue
        #scaler = StandardScaler()
        #scaler = MinMaxScaler()
        scaler = MaxAbsScaler()
        scaler.fit(num_dataset)
        I = scaler.transform(num_dataset)
        X = pd.DataFrame(I, columns=num_dataset.columns)
        X.set_index(dataset.index, inplace =True)
        
        dataset = pd.concat([X,dataset.select_dtypes(['object']),dataset[refill]], axis=1)

        OUTPUT += '.NORM'
        dataset_to_print=pd.concat([dataset,fulldataset['CLASS']],axis=1)
        try:
            dataset.set_index(['ID'], inplace =True)
            dataset_to_print.set_index(['ID'], inplace =True)
        except:
            pass
        dataset_to_print.to_csv(path_or_buf=OUTPUT+'.csv',sep='\t')
        
        print_stats(OUTPUT+'.csv')


# 3) Utilizzo un LabelBinarizer per splittare le features di tipo object in features binarie (solo FORMAT-GT)
    if 'LB' in pipe:
        lb = LabelBinarizer()
        #print dataset.index
        FORMATGT_lb = pd.DataFrame(lb.fit_transform(dataset['FORMAT-GT']), columns=['FORMAT-GT-0/0','FORMAT-GT-0/1','FORMAT-GT-1/1'], index= dataset.index)       
        dataset = dataset.drop(['FORMAT-GT'], axis=1)
        dataset = pd.concat([dataset, FORMATGT_lb], axis=1)


        IEVAISR_lb = pd.DataFrame(lb.fit_transform(dataset['IEVA-iSR']), columns=['IEVA-SR-0','IEVA-SR-SRS','IEVA-SR-HP'], index= dataset.index)
        dataset = dataset.drop(['IEVA-iSR'], axis=1)
        dataset = pd.concat([dataset,IEVAISR_lb ], axis=1)

        OUTPUT += '.LB'
        dataset_to_print=pd.concat([dataset,fulldataset['CLASS']],axis=1)
        try:
            dataset.set_index(['ID'], inplace =True)
            dataset_to_print.set_index(['ID'], inplace =True)
        except:
            pass
        dataset_to_print.to_csv(path_or_buf=OUTPUT+'.csv',sep='\t')
        print_stats(OUTPUT+'.csv')


# 
    if 'DROP' in pipe:

        TODROP = []

        #FEATURES COMPOSTE DA CONTEGGIO DI READS
        #TODROP += ['IEVA-iAD REF', 'IEVA-iAD ALT']
        # TODROP += ['IEVA-iAD REF']
        #TODROP += ['IEVA-iSBD RF', 'IEVA-iSBD RR', 'IEVA-iSBD AF', 'IEVA-iSBD AR'] # sono utilizzati nello SB
        # TODROP += ['IEVA-iAXS0 ALT', 'IEVA-iACR ALT', 'IEVA-iAMQ0 ALT']
        # TODROP += ['IEVA-iAXS0 REF', 'IEVA-iACR REF', 'IEVA-iAMQ0 REF']

        # #FEATURES DI TIPO CATEGORICO INUTILE AI FINI DEL ML
        TODROP += ['FB_INFO-TYPE','IEVA-iSRU','IEVA-iVC'] #sequenza ripetuta nei SR e tipo di variante (snv o indel)

        # #FEATURES SEMPRE 0
        TODROP += ['IEVA-iUnMap', 'IEVA-iSA', 'IEVA-iNP','IEVA-iAMQ0-REF','IEVA-iAMQ0-ALT']
        TODROP += ['IEVA-iACR', 'IEVA-iAXS', 'IEVA-iAXS0', 'IEVA-iAMQ0']
        TODROP += ['IEVA-iQual','IEVA-iAMQ','IEVA-iAAS']


        # #FEATURES STRETTAMENTE LEGATE AD ALTRE
        # TODROP += ['FORMAT-GT-0/0', 'FORMAT-GT-0/1', 'FORMAT-GT-1/1'] #strettamente legata a FORMAT-PL
        #TODROP += ['IEVA-iDP', 'IEVA-iFREQ']
        # #TODROP += ['FORMAT-DP', 'FORMAT-AF']

        # #FEATURES CON INFO GAIN CIRCA 0
        # #TODROP += ['IEVA-SR-HP', 'IEVA-SR-0']

        
        # TODROP += ['INFO-MLEAC','INFO-MLEAF']
        # TODROP += ['IEVA-iMQ','IEVA-iQual REF','IEVA-iAMMQ REF','IEVA-iANRP','IEVA-iAS','IEVA-iBQVA','']

        # #LE TENGO E VEDO DOPO SE FILTRARLE
        # #TODROP += ['IEVA-iAMMQ REF', 'IEVA-iQual REF', 'IEVA-iAAS REF', 'IEVA-iAXS REF'] # missing values > 15%
        # #TODROP += ['IEVA-iAMMQ ALT', 'IEVA-iQual ALT', 'IEVA-iAAS ALT', 'IEVA-iAXS ALT'] # perche abbiamo tolto i ref?
        # #TODROP += ['IEVA-iCRT', 'IEVA-iQRT', 'IEVA-iMRT', 'IEVA-iPRT']  # missing values > 15%


        #FEATURES DA ELIMINARE DA CLUSTERMAP

        # TODROP += ['IEVA-iUnMap', 'IEVA-iSA', 'IEVA-iNP', 'IEVA-iAS','IEVA-iRM',]
        # TODROP += ['IEVA-iACR', 'IEVA-iAXS', 'IEVA-iAXS0', 'IEVA-iAMQ0']
        # TODROP += ['IEVA-iMIU', 'IEVA-iNPA', 'IEVA-iMQ0', 'IEVA-iAMQ','IEVA-iDUP', 'IEVA-iNPP']
        # TODROP += ['IEVA-iAD REF','IEVA-iQual REF']
        # TODROP += ['IEVA-iAAS REF','IEVA-iAAS ALT']
        # TODROP += ['IEVA-iSBD RF', 'IEVA-iSBD RR', 'IEVA-iSBD AF', 'IEVA-iSBD AR'] # sono utilizzati nello SB
       # TODROP += ['IEVA-iAXS0 ALT','IEVA-iAXS ALT', 'IEVA-iACR ALT', 'IEVA-iAMQ0 ALT','IEVA-iAMMQ ALT']
       # TODROP += ['IEVA-iAXS0 REF','IEVA-iAXS REF', 'IEVA-iACR REF', 'IEVA-iAMQ0 REF','IEVA-iAMMQ REF']
        # TODROP += ['IEVA-iDP', 'IEVA-iFREQ','IEVA-iMQ']
        # TODROP += ['FORMAT-RO','IEVA-iAD ALT']

        if opts.droplist is not None:
            droplist += [(f.rstrip() for f in open(opts.droplist,'r').readlines())]
            dataset = dataset.drop(droplist, axis=1)
        else:
            # try:
            #     dataset = drop_pseudonucleotidi(dataset)
            # except:
            #     pass
            for f in TODROP:
                try:
                    dataset = dataset.drop(f,axis=1)
                except:
                    continue

        dataset_to_print=pd.concat([dataset,fulldataset['CLASS']],axis=1)
        #OUTPUT += '.DROP.NOPSNC'
        OUTPUT += '.DROP'
        try:
            dataset.set_index(['ID'], inplace =True)
            dataset_to_print.set_index(['ID'], inplace =True)
        except:
            pass
        dataset_to_print.to_csv(path_or_buf=OUTPUT+'.csv',sep='\t')
        print_stats(OUTPUT+'.csv')


# 4) Elimina tutte le features di ieva cosi da avere un dataset no ieva
    if 'NOIEVA' in pipe:
        dataset=pd.concat([dataset,fulldataset['CLASS']],axis=1)
        try:
            dataset.set_index(['ID'], inplace =True)
        except:
            pass
        train_noieva = dataset.copy()

        for f in dataset.columns.values.tolist():
            if f.startswith('IEVA'):
                train_noieva = train_noieva.drop([f],axis=1)

        dataset.to_csv(path_or_buf=OUTPUT+'.csv',sep='\t')
        train_noieva.to_csv(path_or_buf=OUTPUT+'.NOIEVA.csv',sep='\t')

# 5)PCA E TSNE per vedere come si distribuiscono le classi sulle componenti principali:
    if 'PCA' in pipe:
        pca = PCA(n_components=0.9999,random_state=64)
        CLASS = fulldataset['CLASS'].reset_index
        X_pca = pca.fit_transform(dataset)
        X_pca = pd.concat([pd.DataFrame(X_pca),pd.DataFrame(dataset.index)], axis=1)
        X_pca.set_index(['ID'], inplace =True)
        df_pca = pd.concat([X_pca.copy(),fulldataset['CLASS']],axis=1)
        df_pca['PCA 0'] = df_pca[0]
        df_pca['PCA 1'] = df_pca[1]
        print ggplot(df_pca,aes(x='PCA 0',y='PCA 1',color='CLASS')) + geom_point() + ggtitle("PCA IEVA")
    
    # 5.1) calcolo i pesi della PCA per vedere quali features sono piu informative
        # for comp in pca.components_:
        #     co =[]
        #     for c in comp:
        #         co += [str(c)]
        #     print '\t'.join(co)

        out = open(OUTPUT +'.PCA-FeaturesVariance.STATS.csv','w')
        out.write('\t'.join(['FEATURE','PCA0','PCA1','PCA2','PCA3','PCA4']) + '\n')

        for f in dataset.columns.values.tolist():
            co=[]
            for comp in pca.components_:
                co += [str(comp[dataset.columns.values.tolist().index(f)])]
            out.write('\t'.join([f]+co) + '\n')
        out.close()

    # 6.2) calcolo del numero di Principal Components necessari ad ottenere il 98% della varianza
        cum_var_exp = np.cumsum(pca.explained_variance_ratio_)
        pca_var = pd.DataFrame(cum_var_exp, columns=['% Variance'])
        pca_var = pd.concat([pd.DataFrame(range(1,len(pca.components_)+1),columns=['Num PC']),pca_var],axis=1)

        print ggplot(pca_var,aes(y='% Variance', x='Num PC', label='% Variance')) + geom_point() \
            + geom_line(color = 'red', alpha = 0.50, size = 2.5) \
            + ggtitle("cumulative explained variance") \
            + geom_text(position='jitter',size=6) \
            #+ geom_text(aes(label = if('% Variance'>=0.95:'% Variance')))
            #+ geom_text(label=ifelse('% Variance'>=0.95 and '% Variance'<0.955),hjust=0,vjust=0)
            

    if 'TSNE' in pipe:
        tsne = TSNE(n_components=2)
        X_tsne = tsne.fit_transform(dataset)
        X_tsne = pd.concat([pd.DataFrame(X_tsne),fulldataset['CLASS']], axis=1)

        df_tsne = pd.DataFrame(X_tsne.copy())
        df_tsne['TSNE 0'] = X_tsne[0]
        df_tsne['TSNE 1'] = X_tsne[1]
        print ggplot(df_tsne,aes(x='TSNE 0',y='TSNE 1',color='CLASS')) + geom_point() + ggtitle("T-SNE IEVA") #+ facet_wrap("CLASS")

    
    if 'CLUSTERING' in pipe:
        print dataset.columns[dataset.nunique() == 1]
        class_map = {'PASS': 1, 'FILTER': 0}
        Y_dataset = fulldataset['CLASS'].map(class_map)

        lut = dict(zip(Y_dataset.unique(), "rbg"))
        row_colors = Y_dataset.map(lut)

        ds = pd.concat([dataset,Y_dataset],axis=1)

        CLUSTERMAP = seaborn.clustermap(dataset,method='average', metric='correlation', row_colors=row_colors, xticklabels=1, cmap='bwr')

        ax = CLUSTERMAP.ax_heatmap
        plt.show()

# 6) Scegliere le features piu informative con INFOGAIN
    if 'IG' in pipe:
    # 6.1) calcolo dell' Info Gain
        out = open(OUTPUT +'.MUTUAL-INFOGAIN','w')
        res = dict(zip(dataset.columns.values.tolist(), mutual_info_classif(dataset, fulldataset['CLASS'], discrete_features=True)))
        
        for r in res.keys():
            out.write(r+'\t'+str(res[r])+'\n')

    if 'IGFilter' in pipe:
        #7/A) FACOLTATIVO: Filtro le features in base a una soglia di varianza: 0.9999
        selector = VarianceThreshold(0.0)
        selector.fit_transform(dataset,fulldataset['CLASS'])
        ds = dataset.loc[:,selector.get_support(indices=False)]
        dataset  = pd.concat([ds,fulldataset['CLASS']], axis=1)

        OUTPUT += '.IGF'
        dataset.to_csv(path_or_buf=OUTPUT+'.csv',sep='\t')
        print_stats(OUTPUT+'.csv')
    
    
    # 8)Calcolo della matrice di correlazione

    if 'CORR' in pipe:

        corr_matrix = dataset.corr()
        labels = dataset.columns.values.tolist()
        fig, ax = plt.subplots()

        cax = ax.matshow(corr_matrix)
        fig.colorbar(cax)
        ticks = np.arange(0,len(dataset.columns),1)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticklabels(labels)
        
        plt.title('Correlation Matrix', y=1.08)
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params(labelsize=7)
        plt.xticks(rotation=90)
        plt.show()
        plt.savefig(OUTPUT + 'CORR-MATRIX.png')

        print_full(corr_matrix.to_csv)

        writer = pd.ExcelWriter(OUTPUT+'.CORR-MATRIX.xlsx')
        corr_matrix.to_excel(writer,'Sheet1')
        writer.save()

        # 8.1) scatter matrix delle prime 10 features
        #scatter_matrix(dataset[labels[:10]],figsize=(12, 8))
        #plt.show()


    #9) Recursive Feature Elimination

    if 'RFE' in pipe:

        out=open(OUTPUT + '.RFE-RESULTS.csv','w')
        out.write('\t'.join(['FEATURE ELIMINATA','TN','FP','FN','TP','AUC','BALANCED ACC','RECALL','PRECISION','F1 SCORE'])+'\n')
        
        class_map = {'PASS': 1, 'FILTER': 0}
        Y_dataset = fulldataset['CLASS'].map(class_map)

        ROC_AUC=[]
        PRECISION = []
        RECALL = []
        BACC = []
        F1 = []
        
        cv = 4

        for feat in dataset.columns:
            ds = dataset.drop(feat,axis=1).reset_index()
            X_dataset = ds.drop(['ID'], axis=1)
            forest_clf = RandomForestClassifier(class_weight = 'balanced', random_state=42)
            #skfolds = StratifiedKFold(n_splits=3, random_state=42)
            y_train_pred = cross_val_predict(forest_clf, X_dataset, Y_dataset, cv=cv)
            tn, fp, fn, tp = confusion_matrix(Y_dataset, y_train_pred).ravel()

            scores = cross_validate(forest_clf, X_dataset, Y_dataset, cv=cv, scoring=['roc_auc','precision','recall','balanced_accuracy','f1_weighted'])
            
            out.write(feat +'\t'+ str(tn) +'\t'+str(fp) +'\t'+ str(fn) +'\t'+ str(tp)+'\t')
            out.write(
                str(scores['test_roc_auc'].mean()) + '\t'+
                str(scores['test_precision'].mean()) +'\t'+
                str(scores['test_recall'].mean()) + '\t'+
                str(scores['test_balanced_accuracy'].mean()) + '\t'+ 
                str(scores['test_f1_weighted'].mean()) + '\n')

            ROC_AUC += [scores['test_roc_auc'].mean()]
            PRECISION += [scores['test_precision'].mean()]
            RECALL += [scores['test_recall'].mean()]
            BACC += [scores['test_balanced_accuracy'].mean()]
            F1 += [scores['test_f1_weighted'].mean()]

        fig, ax = plt.subplots()
        ax.plot(BACC, label='Balanced Accuracy')
        ax.plot(ROC_AUC, label='ROC auc')
        ax.plot(RECALL, label='Recall')
        ax.plot(PRECISION, label='Precision')
        ax.set(xlabel='Eliminated Feature',
            title='Recursive Feature Elimination')
        xticks = np.arange(0,len(dataset.columns),1)
        yticks = np.arange(0.95,1.000,0.001)
        
        ax.set_xticks(xticks)
        ax.set_xticklabels(dataset.columns.tolist())
        plt.xticks(rotation=90)

        ax.set_yticks(yticks)

        ax.grid()
        plt.legend()

        plt.show()
        out.close()
        plt.savefig(OUTPUT + 'RFE.png')

        out=open(OUTPUT + '.RFE-RANKING-RESULTS.csv','w')
        forest_clf = RandomForestClassifier(class_weight = 'balanced', random_state=42)
        selector = RFE(forest_clf, step=1)
        selector = selector.fit(dataset, Y_dataset)
        for f in dataset.columns.tolist():
            out.write(f +'\t'+ str(selector.ranking_[ dataset.columns.tolist().index(f)])+'\n')


    if 'RFECV' in pipe:
        
        X_dataset = dataset.copy()
        class_map = {'PASS': 1, 'FILTER': 0}
        Y_dataset = fulldataset['CLASS'].map(class_map)
        # Create the RFE object and compute a cross-validated score.
        forest_clf = RandomForestClassifier(class_weight = 'balanced')
        # The "accuracy" scoring is proportional to the number of correct
        # classifications
        rfecv = RFECV(estimator=forest_clf, step=1, cv=StratifiedKFold(4),
                      scoring='balanced_accuracy')
        rfecv.fit(X_dataset, Y_dataset)

        

        # Plot number of features VS. cross-validation scores
        fig, ax = plt.subplots()
        ax.set(xlabel="Number of features selected",
            ylabel= "Cross validation score (Balanced Accuracy)",
            title='Cross Validation RFE')

        ax.set_xticks(np.arange(1, len(rfecv.grid_scores_) + 1,1))
        ax.set_yticks(np.arange(min(rfecv.grid_scores_)-0.002,1.000,0.001))
        ax.plot(rfecv.grid_scores_)
        plt.grid()
        plt.show()

        plt.savefig(OUTPUT + 'RFECV.png')

        out=open(OUTPUT + '.RFECV-RANKING-RESULTS.csv','w')
        out.write("Optimal number of features :\t %d \n" % rfecv.n_features_)
        for f in dataset.columns.tolist():
            out.write(f +'\t'+ str(rfecv.ranking_[ dataset.columns.tolist().index(f)])+'\n')

    
    #9)  k-fold cross validation utilizzando la random forest  

    if 'TUNING' in pipe:

        X_dataset = dataset.copy()
        class_map = {'PASS': 1, 'FILTER': 0}
        Y_dataset = fulldataset['CLASS'].map(class_map)

        ROC_AUC=[]
        PRECISION = []
        RECALL = []
        BACC = []
        F1 = []    
        
        out=open(OUTPUT + '.TUNING-RESULTS.csv','w')
        #out.write('NUM TREE\n')
        out.write('\t'.join(['N ESTIMATORS','TN','FP','FN','TP','AUC','BALANCED ACC','RECALL','PRECISION','F1 SCORE'])+'\n')
        cv = 4

        for n_tree in range(421)[10::30]:
            forest_clf = RandomForestClassifier(class_weight = 'balanced',n_estimators = n_tree,  random_state=42)
            
            y_train_pred = cross_val_predict(forest_clf, X_dataset, Y_dataset, cv=cv)
            tn, fp, fn, tp = confusion_matrix(Y_dataset, y_train_pred).ravel()
          
            scores = cross_validate(forest_clf, X_dataset, Y_dataset, cv=cv, scoring=['roc_auc','precision','recall','balanced_accuracy','f1_weighted'])
        
            out.write(str(n_tree)+'\t'+ str(tn) +'\t'+str(fp) +'\t'+ str(fn) +'\t'+ str(tp)+'\t'+
                str(scores['test_roc_auc'].mean()) + '\t'+
                str(scores['test_precision'].mean()) +'\t'+
                str(scores['test_recall'].mean()) + '\t'+
                str(scores['test_balanced_accuracy'].mean()) + '\t'+ 
                str(scores['test_f1_weighted'].mean()) + '\n')

            ROC_AUC += [scores['test_roc_auc'].mean()]
            PRECISION += [scores['test_precision'].mean()]
            RECALL += [scores['test_recall'].mean()]
            BACC += [scores['test_balanced_accuracy'].mean()]
            F1 += [scores['test_f1_weighted'].mean()]

        fig, ax = plt.subplots()
        ax.plot(BACC, label='Balanced Accuracy')
        ax.plot(ROC_AUC, label='ROC auc')
        ax.plot(RECALL, label='Recall')
        ax.plot(PRECISION, label='Precision')
        ax.set(xlabel='Num Estimators',
            title='Number of Estimators Tuning')

        xticks = np.arange(0,len(range(421)[10::30]),1)
        yticks = np.arange(0.96,1.000,0.001)

        ax.set_xticks(xticks)
        ax.set_xticklabels(range(421)[10::30])
        plt.xticks(rotation=90)

        ax.set_yticks(yticks)
        ax.grid()
        plt.legend()

        plt.show()
        plt.savefig(OUTPUT + '.N_tree-TUNING.png')


        # Number of trees in random forest
        n_estimators = [int(x) for x in np.linspace(start = 20, stop = 400, num = 20)]
        #n_estimators = [280]
        # Number of features to consider at every split
        max_features = ['auto', 'sqrt']
        # Maximum number of levels in tree
        max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
        max_depth.append(None)
        # Minimum number of samples required to split a node
        min_samples_split = [2, 5, 10]
        # Minimum number of samples required at each leaf node
        min_samples_leaf = [1, 2, 4]
        # Method of selecting samples for training each tree
        bootstrap = [True, False]

        # Create the random grid
        random_grid = {'n_estimators': n_estimators,
                       'max_features': max_features,
                       'max_depth': max_depth,
                       'min_samples_split': min_samples_split,
                       'min_samples_leaf': min_samples_leaf,
                       'bootstrap': bootstrap}


        rf = RandomForestClassifier()
        rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_grid, n_iter = 100, cv = 4, verbose=0, random_state=42, n_jobs = -1)
        rf_random.fit(X_dataset, Y_dataset)
        print rf_random.best_params_
        BP =rf_random.best_params_

        #grid_search = GridSearchCV(estimator = rf, param_grid = random_grid, cv = 4, n_jobs = -1, verbose = 0)
        #grid_search.fit(X_dataset, Y_dataset)
        #print grid_search.best_params_

        out.write('BEST HYPERPARAMETERS (RANDOM GRID):'+'\n'+
                        'n_estimators' +'\t' + str(BP['n_estimators']) +'\n'
                        'max_features' +'\t' + str(BP['max_features']) +'\n'
                        'max_depth' +'\t' + str(BP['max_depth']) +'\n'
                        'min_samples_split' +'\t' + str(BP['min_samples_split']) +'\n'
                        'min_samples_leaf' +'\t' + str(BP['min_samples_leaf']) +'\n'
                        'bootstrap' +'\t' + str(BP['bootstrap']))
        
        #grid_search = GridSearchCV(estimator = rf, param_grid = random_grid, cv = 4, n_jobs = -1, verbose = 0)
        #grid_search.fit(X_dataset, Y_dataset)
        #print grid_search.best_params_


    if 'KFOLD' in pipe:
        X_dataset = dataset.copy()
        class_map = {'PASS': 1, 'FILTER': 0}
        Y_dataset = fulldataset['CLASS'].map(class_map)

        # parametri tunati per dataset IEVA
        # n_estimators = 100
        # max_features = 'auto'
        # max_depth = None
        # min_samples_split = 5
        # min_samples_leaf = 1
        # bootstrap = True

        # ## parametri tunati per dataset NOIEVA
        n_estimators = 120
        max_features = 'sqrt'
        max_depth = 70
        min_samples_split = 5
        min_samples_leaf = 4
        bootstrap = False

        cv = 4
        out=open(OUTPUT + '.KFOLD-RESULTS.csv','w')

        RF_ARRAY = []
        BACC_ARRAY = []
        for j in range(10):
            forest_clf = RandomForestClassifier(class_weight = 'balanced',
                                                n_estimators = n_estimators,
                                                max_features = max_features,
                                                max_depth = max_depth,
                                                min_samples_split = min_samples_split,
                                                min_samples_leaf = min_samples_leaf,
                                                bootstrap = bootstrap)

            out.write('\t'.join(['N REPETITION','TN','FP','FN','TP','AUC','BALANCED ACC','RECALL','PRECISION','F1 SCORE'])+'\n')
            y_train_pred = cross_val_predict(forest_clf, X_dataset, Y_dataset, cv=cv)
            scores = cross_validate(forest_clf, X_dataset, Y_dataset, cv=cv, scoring=['roc_auc','precision','recall','balanced_accuracy','f1_weighted'])
            
            tn, fp, fn, tp = confusion_matrix(Y_dataset, y_train_pred).ravel()

            out.write(str(j)+'\t'+ str(tn) +'\t'+str(fp) +'\t'+ str(fn) +'\t'+ str(tp)+'\t'+
                str(scores['test_roc_auc'].mean()) + '\t'+
                str(scores['test_precision'].mean()) +'\t'+
                str(scores['test_recall'].mean()) + '\t'+
                str(scores['test_balanced_accuracy'].mean()) + '\t'+ 
                str(scores['test_f1_weighted'].mean()) + '\n')

            result = pd.concat([pd.DataFrame(Y_dataset).reset_index(), pd.DataFrame(y_train_pred,columns=['PRED'])], axis=1)

            FN = result[result['CLASS'] == 1]
            FN = FN[FN['PRED'] == 0]
            FP = result[result['CLASS'] == 0]
            FP = FP[FP['PRED'] == 1]
            #print FP.apply('\t'.join)
            FN = re.sub('\n\t','\n',re.sub(' ','\t',re.sub("[\[\]\']","",np.array2string(FN.values))))
            FP = re.sub('\n\t','\n',re.sub(' ','\t',re.sub("[\[\]\']","",np.array2string(FP.values))))

            RF_ARRAY += [forest_clf]
            BACC_ARRAY += [scores['test_balanced_accuracy'].mean()]
            #out.write(np.array2string(FN.values) +'\n')
            out.write('ID' +'\t' +'CLASS' + '\t' + 'PRED'+'\n')
            out.write(FN+'\n'+ FP+'\n')

            dump(RF_ARRAY[BACC_ARRAY.index(max(BACC_ARRAY))], OUTPUT+'.BESTCV_RFclassifier.joblib') 


    if opts.train is not None:

        TRAIN=pd.read_csv(opts.train,sep='\t', index_col=['ID'])
        X_train = TRAIN.drop(['CLASS'], axis=1)
        y_train = TRAIN['CLASS']
        forest_clf = RandomForestClassifier(class_weight = 'balanced',n_estimators = 50)
        forest_clf.fit(X_train, y_train)
        dump(forest_clf, OUTPUT+'.RFclassifier.joblib')

        if opts.test is not None:

            TEST=pd.read_csv(opts.test,sep='\t', index_col=['ID'])
            X_test = TEST.drop(['CLASS'], axis=1)
            y_test = TEST['CLASS']
            y_pred = forest_clf.predict(X_test)

            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
            out=open(OUTPUT + '.TEST-EVALUATION-RESULTS.csv','w')
            out.write('CONFUSION MATRIX TEST SET:'+'\n')
            out.write('\t' +str(tn) +'\t'+str(fp) +'\n\t'+ str(fn) +'\t'+ str(tp)+'\n')

    if opts.model is not None:

        clf = load(opts.classifier)

        TEST=pd.read_csv(opts.test,sep='\t', index_col=['ID'])
        X_test = TEST.drop(['CLASS'], axis=1)
        y_test = TEST['CLASS']
        y_pred = clf.predict(X_test)

        tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
        out=open(OUTPUT + '.TEST-EVALUATION-RESULTS.csv','w')
        out.write('CONFUSION MATRIX TEST SET:'+'\n')
        out.write(str(tn) +'\t'+str(fp) +'\n\t'+ str(fn) +'\t'+ str(tp)+'\n')
