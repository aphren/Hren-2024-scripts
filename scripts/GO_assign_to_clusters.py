# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 14:21:30 2023

@author: andre
"""

import pandas as pd
import re

data = pd.read_excel('uniprot_data.xlsx')
data['locus'] = data['Gene Names (ordered locus)'].str.split('_').str[1]
data['gene'] = data['Gene Names'].str.split(' ').str[0]
data['Gene Ontology (biological process)'] = data['Gene Ontology (biological process)'].astype(str)
unique_GO = data['Gene Ontology (biological process)'].str.cat(sep=';').split(';')
unique_GO = [entry.lstrip() for entry in unique_GO]
GO_loci = pd.DataFrame(zip(data['locus'],data['Gene Ontology (GO)'],data['Gene Ontology (biological process)']),columns=['locus','GO1','GO2'])

#unique_GO = list(set(unique_GO))  
    
cluster_data = pd.read_excel('2023_07_19 clusters.xlsx')    
merged = pd.merge(cluster_data,GO_loci,on='locus',how='left')
merged.to_csv('GO_clusters.csv')    


output = pd.DataFrame.from_dict(final, orient='index').reset_index()
output['gene_count'] = output.count(axis=1)

# Move column 'C' to the second position
column_name = 'gene_count'
second_position = 1

# Extract the column to be moved
column_to_move = output[column_name]

# Drop the column from its current position
output = output.drop(column_name, axis=1)

# Insert the column at the desired position
output.insert(second_position, column_name, column_to_move)

output = output.drop(output[output['index'] == 'nan'].index)
output = output.dropna(axis=1, how='all')

#output.to_csv('PCC7002_GO_terms.csv')