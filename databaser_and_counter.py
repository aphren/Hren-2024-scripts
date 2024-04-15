# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 09:18:53 2022

@author: andre
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import math
import os 
import statistics

def flip_spacer(row):
    spacer = Seq(str(row.Spacer)).reverse_complement()
    return spacer

def database_maker(spacer_dict):
    database = []
    for key in spacer_dict.keys():
        name = '>' + key
        database.append(name)
        if key[-1] == 'r':
            seq_rev = seq1 + str(spacer_dict[key]) + seq2
            database.append(seq_rev)
        elif key[-1] == 'f':
            seq_fwd = seq3 + str(spacer_dict[key]) + seq4
            database.append(seq_fwd)
    return database

def value_counter(sample):
    # Count reads basically
    
    # If the file is >50million reaads you might get a memory error and 
    # will have to split the file
    filename = './complete/' + sample + '_aligned_id98.tsv'
    print('Processing ',filename)
    instances = []

    with open(filename) as f:
        for line in f:
            entry = line.split()
            instances.append(entry[9][:-1])
        df = pd.Series(instances)
        counts = pd.DataFrame(df.value_counts())
        counts.reset_index(inplace=True)
        counts.columns = ['ID','counts']
        name = sample + '_counts.csv'
        counts.to_csv(name)
            
def freq(row):
    frequency = int(row.counts) / int(row.total)
    return frequency

def lfc(row):
    lfc = math.log(row.freq_y/row.freq_x,2)
    return lfc

def strip_locus(row):
    info = row.locus.split('_')[1]
    return info

def strip_locus2(row):
    info = row.ID.split('_')[0]
    return info

## Creating putida database from Jacob's library csv



"""
## This code reads all variants of my 7002 library, then arranges the 
## spacer IDs and spacers into a fasta database for global alignment (vsearch)

variants = pd.read_csv('variants.csv',low_memory=False)
variants['reverse_spacer'] = variants.apply(flip_spacer,axis=1)

ids = variants['ID'].to_list()
spacers_fwd = variants['Spacer'].to_list()
spacers_rev = variants['reverse_spacer'].to_list()
spacer_dict = {}
for x in range(len(ids)):
    rev_id = ids[x] + 'r'
    fwd_id = ids[x] + 'f'
    spacer_dict[rev_id] = str(spacers_rev[x])
    spacer_dict[fwd_id] = str(spacers_fwd[x])

database = database_maker(spacer_dict)

with open('database2.fasta','w') as fp:
    fp.write("\n".join(database))
   

## This code reads the tsv file output from alignment and counts unique instances

samples = []
for x in os.listdir('./complete'):
    if x[-4:] == '.tsv':
        samples.append(x[:-17])


for sample in samples:
    value_counter(sample)





## Preparing files for DESeq2 input
samples = ['37B1','37B2','dCasLib1','dCasLib2']
dataframes = []
combined = pd.DataFrame()
for sample in samples:
    file = './counts/' + sample + '_counts.csv'
    data = pd.read_csv(file,names=['ID','counts'],header=0)
    count_name = 'counts_' + sample
    data.rename(columns={'counts':count_name},inplace=True)
    if len(combined) == 0:
        combined = data.copy()
    else:
        combined = pd.merge(combined,data,on='ID',how='outer')
combined.fillna(value=0,inplace=True)
combined.to_csv('cts_37B.csv')


## combine deseq2 and crisphiermix outputs into final data with descriptions
data = pd.read_csv('./condition_files/37W/condition_library_vs_37W.csv')
data2 = pd.read_csv('./condition_files/37W/crisphiermix_37W_output.csv')
full = pd.merge(data,data2,how='outer',on='gene')
full = full.rename(columns={'gene':'locus'})
gene_data = pd.read_csv('gene_descriptions.csv',usecols=['locus','product','note','gene'])   
gene_data['locus'] = gene_data.apply(strip_locus,axis=1)
full = pd.merge(full,gene_data,how='left',on='locus')
full.to_csv('37W_final.csv')


## create average LFC for each locus
## then output for comparisons
data = pd.read_csv('./condition_files/22B/condition_library_vs_22B.csv')
data = data.rename(columns={'Unnamed: 0':'ID'})
IDs = data['ID'].to_list()
genes = []

for ID in IDs:
    locus = ID.split('_')[0]
    genes.append(locus)
genes = list(set(genes))
average_lfcs = []
for gene in genes:    
    relevant = data[data['gene'] == gene]
    avg = relevant['LFC'].mean() 
    average_lfcs.append(avg)

full = pd.DataFrame(genes)
full['avg_lfc'] = average_lfcs
full.to_csv('22B_avg_lfcs.csv')


## combine avg lfcs from two conditions for graphing purposes
data1 = pd.read_csv('22B_avg_lfcs.csv')
data1.rename(columns={'avg_lfc':'avg_lfc_22B'},inplace=True)
data2 = pd.read_csv('./condition_files/22R/22R_avg_lfcs.csv')
data2.rename(columns={'avg_lfc':'avg_lfc_22R'},inplace=True)
full = pd.merge(data1,data2,on='locus')
full.to_csv('22B vs 22R avg lfcs.csv')

"""
## create a master file with all LFC data for each guide
conditions = ['37W','37R','37B','22W','22R','22B','37di']
passes = 0
for condition in conditions:
    path = './condition_files/' + condition + '/' + condition + '_final.csv'
    if passes == 0:
        full = pd.read_csv(path,header=0,usecols=['ID','LFC'])
        new_LFC = condition + '_LFC'
        full.rename(columns={'LFC':new_LFC},inplace=True)
    else:
        temp = pd.read_csv(path,header=0,usecols=['ID','LFC']) 
        new_LFC = condition + '_LFC'
        temp.rename(columns={'LFC':new_LFC},inplace=True)
        full = pd.merge(full,temp,on='ID')
    passes += 1

full['locus'] = full.apply(strip_locus2,axis=1)
gene_data = pd.read_csv('full_descriptions.csv',usecols=['locus','product','note','gene'])   
gene_data['locus'] = gene_data.apply(strip_locus,axis=1)
full = pd.merge(full,gene_data,how='left',on='locus')
full.to_csv('full_LFC_data.csv')

"""
## ADDING TRNAS TO GENE DESCRIPTIONS
gene_data = pd.read_csv('gene_descriptions.csv',usecols=['locus','product','note','gene'])   
#gene_data['locus'] = gene_data.apply(strip_locus,axis=1)
idk = pd.read_csv('trna_data.txt',sep='\t')
idk.rename(columns={'Aliases':'locus','description':'product'},inplace=True)
nontrna = gene_data['locus'].to_list()
trna = idk['locus'].to_list()
missing = set(trna).difference(set(nontrna))
idk2 = idk[['locus','product']].copy()

descriptions = pd.merge(gene_data,idk2,on=['locus','product'],how='outer')
#descriptions.to_csv('full_descriptions.csv')
"""