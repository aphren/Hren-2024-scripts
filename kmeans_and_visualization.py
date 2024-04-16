# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 15:05:56 2023

@author: andre
"""

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler,MinMaxScaler
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

data_22b = pd.read_excel("22B_final.xlsx")
#data_22b = data_22b.drop(data_22b[data_22b['padj'] > 0.05].index)
data_22b = data_22b.drop(data_22b[data_22b['locus'] == 'nt-gRNA'].index)
data_22b = data_22b.drop(data_22b[data_22b['locus'] == 'dup'].index)

data_22r = pd.read_excel("22R_final.xlsx")
#data_22r = data_22r.drop(data_22r[data_22r['padj'] > 0.05].index)
data_22r = data_22r.drop(data_22r[data_22r['locus'] == 'nt-gRNA'].index)
data_22r = data_22r.drop(data_22r[data_22r['locus'] == 'dup'].index)

data_37w = pd.read_excel("37W_final.xlsx")
#data_22r = data_22r.drop(data_22r[data_22r['padj'] > 0.05].index)
data_37w = data_37w.drop(data_37w[data_37w['locus'] == 'nt-gRNA'].index)
data_37w = data_37w.drop(data_37w[data_37w['locus'] == 'dup'].index)

"""
X = 5
counts = data_22b['locus'].value_counts()
mask = data_22b['locus'].isin(counts[counts >= X].index)
data_22b = data_22b[mask]

X = 5
counts = data_22r['locus'].value_counts()
mask = data_22r['locus'].isin(counts[counts >= X].index)
data_22r = data_22r[mask]
"""

### metric generator




data_22b = data_22b[['gRNA','LFC','locus','padj']]
data_22b = data_22b.rename(columns={'LFC': 'LFC_22B','gRNA':'ID','padj':'padj_22B'})

data_22r = data_22r[['ID','LFC','locus','padj']]
data_22r = data_22r.rename(columns={'LFC': 'LFC_22R','padj':'padj_22R'})

data_37w = data_37w[['gRNA','LFC','locus','padj']]
data_37w = data_37w.rename(columns={'LFC': 'LFC_37W','gRNA':'ID','padj':'padj_37W'})


data = pd.merge(data_22b,data_22r,on=['ID','locus'],how='outer')
data = pd.merge(data,data_37w,on=['ID','locus'],how='outer')
data = data[['ID','locus','LFC_22B','LFC_22R','LFC_37W','padj_22B','padj_22R','padj_37W']]
data = data[~((data['padj_22B'] > 0.05) & (data['padj_22R'] > 0.05) & (data['padj_37W'] > 0.05))]

data = data.dropna()
X = 5
counts = data['locus'].value_counts()
mask = data['locus'].isin(counts[counts >= X].index)
data = data[mask]


#data = data_22b
loci = list(set(data['locus'].to_list()))

metrics = 2
X = np.zeros((len(loci),metrics*3))

w37_avg_lfc = []
r22_avg_lfc = []
b22_avg_lfc = []

row = 0
for locus in loci:
    current = data[data['locus'] == locus]
    #sort by 37W order
    #current = current.sort_values('37W_LFC',ascending=True)
    for p in range(3):
        ## Access a single condition (e.g. 37W)
        #if p == 0:
            #x_data = current.iloc[:,p+2]
            #X[row,0] = (x_data.mean())**2
        
        relevant = current.iloc[:,p+2]
        #X[row,0+p*metrics] = (relevant.nlargest(2).mean())**2             #= linregress(x_data,relevant).slope
        #X[row,0+p*metrics] = relevant.std()
        X[row,0+p*metrics] = (relevant > 1).sum()
        X[row,1+p*metrics] = relevant.mean()
        #X[row,2+p*metrics] = (relevant.nlargest(2).mean() - relevant.nsmallest(2).mean())**2
        
        if p==0:
            b22_avg_lfc.append(relevant.mean())
        elif p==1:
            r22_avg_lfc.append(relevant.mean())
        elif p==2:
            w37_avg_lfc.append(relevant.mean())
        
    row += 1
    
locus_means = pd.read_excel('full_locus_means.xlsx')
relevant_means = pd.DataFrame(loci,columns=['locus'])
locus_means = pd.merge(relevant_means,locus_means,on='locus',how='left')
## normalize for color mapping
columns_to_norm = ['37W_mean','22R_mean','22B_mean']
scaler = MinMaxScaler()
locus_means_out = locus_means.copy()

locus_means[columns_to_norm] = scaler.fit_transform(locus_means[columns_to_norm])



    
output = pd.DataFrame(X)
output['locus'] = loci
scale = StandardScaler()
X_norm = scale.fit_transform(X)
     
cluster_count = 5

km = KMeans(n_clusters=cluster_count, n_init=200, algorithm="elkan",random_state=42)  
km.fit(X_norm)
prediction = km.predict(X_norm)
df_km = pd.DataFrame(prediction)
df_km['locus'] = loci
df_km.rename(columns={0:'cluster'},inplace=True)

locus_means_out['cluster'] = prediction
locus_means_out.to_csv('clusters_and_means.csv')


cluster_sizes = df_km.groupby('cluster').count()
print(cluster_sizes)
# Create a bar graph
plt.figure()
plt.bar(cluster_sizes.index, cluster_sizes['locus'])
plt.xlabel('Cluster')
plt.ylabel('Size')
plt.title('Cluster Sizes')
plt.show()


### input variables for 
perp = 35
ee = 4


tsne = TSNE(n_components=2, random_state=42, perplexity=perp, early_exaggeration=ee)#, method='exact')
X_tsne = tsne.fit_transform(X_norm)
df_km['x'] = X_tsne[:, 0]
df_km['y'] = X_tsne[:, 1]

desc = pd.read_excel('uniprot_gene_descriptions.xlsx')
desc = desc[['locus','Gene Names (primary)']]
mixed = pd.merge(df_km,desc,on='locus',how='left')

#mixed[['v1','v2','v3','v4','v5','v6','v7','v8','v9']] = X

mixed = mixed.assign(W37=w37_avg_lfc, R22=r22_avg_lfc,B22=b22_avg_lfc)


#mixed.to_csv('clusters.csv',index=False)

# Map the cluster assignments back to the original data
color_dict = {0:'#FFA040',1:'#F94040',2:'#5757F9',3:'#B856D7',4:'#A0FFA0'}vv
colors = [color_dict[category] for category in prediction]
plt.figure(figsize=(10, 8))
for category, color in color_dict.items():
    mask = (prediction == category)
    plt.scatter(X_tsne[:, 0][mask], X_tsne[:, 1][mask], c=color)

plt.colorbar(ticks=range(cluster_count))

plt.savefig('tsne_clusters_perp' + str(perp) + '_ee' + str(ee) + '.png',dpi=1200,transparent=True,bbox_inches='tight')
#plt.title('PCC 7002 Gene Clusters')
#plt.xlabel('tSNE1')
#plt.ylabel('tSNE2')
#plt.savefig('t-SNE-figure.png',dpi=1200)


# Modified tSNE with color mapping based on gene means
#plt.show()

# Plotting
plots = ['37W_mean','22R_mean','22B_mean']
for x in range(3):
    cmap = plt.get_cmap('inferno')
    #norm = plt.Normalize(0, 1)
    
    
    plt.figure(figsize=(10, 8))
    plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=locus_means[plots[x]], cmap=cmap)#,norm=norm)
    plt.colorbar(ticks=range(cluster_count))
    plt.savefig('tSNE_' + plots[x] + '_perp' + str(perp) + '_ee' + str(ee) + '.png',dpi=1200,transparent=True,bbox_inches='tight')
    #plt.title('PCC 7002 Gene Clusters')
    #plt.xlabel('tSNE1')
    #plt.ylabel('tSNE2')



## Creating a visualizaiton using PCA (2D)
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_norm)
pca_data = pd.DataFrame(X_pca, columns=['PCA1','PCA2'])
pca_data['Cluster'] = prediction

plt.figure(figsize=(8,6))
for cluster_label in pca_data['Cluster'].unique():
    cluster_data = pca_data[pca_data['Cluster'] == cluster_label]
    plt.scatter(cluster_data['PCA1'],cluster_data['PCA2'],label=f'Cluster {cluster_label}')
plt.xlabel('PCA Component 1')
plt.ylabel('PCA Component 2')
plt.title('Data with Assigned Clusters')
plt.legend()
#plt.show()
#plt.savefig('Cluster_PCA.png',dpi=1200)