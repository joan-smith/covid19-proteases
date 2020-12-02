#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 18:37:28 2020

@author: Joan Smith
Airway Atlas: https://www.genomique.eu/cellbrowser/HCA/
"""

#%%
import pandas as pd
from matplotlib import pyplot as plt
import scanpy as sc
import matplotlib as mpl
import numpy as np

import os

#%% Set scanpy and matplotlib settings
sc.settings.figdir = 'data/protease-expression/cluster-plots/airway-atlas/'
sc.settings.verbosity = 4

mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams.update({
    'font.sans-serif': 'Arial',
    'font.family': 'sans-serif',
    'axes.titlesize': 18,
    'axes.labelsize': 18,
    'font.size': 18,
    })


#%%
df = pd.read_csv('data/raw-data/airway-atlas/exprMatrix.tsv', sep='\t', index_col=0,)

#%%
adata = sc.AnnData(df.T)
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.highly_variable_genes(adata)

sc.tl.pca(adata)
sc.pp.neighbors(adata)

#%%
LEARNING_RATE = 1000
EARLY_EXAGGERATION = 12
RESOLUTION = 1.5
PERPLEXITY=130

sc.tl.tsne(adata, learning_rate=LEARNING_RATE, n_jobs=8, early_exaggeration=EARLY_EXAGGERATION, perplexity=PERPLEXITY)
sc.tl.leiden(adata, resolution=RESOLUTION)
#%%
params = {'learning_rate': LEARNING_RATE,
          'early_exaggeration':EARLY_EXAGGERATION,
          'resolution': RESOLUTION,
          'perplexity': PERPLEXITY,
          'genes': 'all'}
pd.Series(params).to_csv(os.path.join(sc.settings.figdir, 'params.txt'))
adata.write(os.path.join(sc.settings.figdir, 'adata.h5ad'))

#%%
markers = pd.read_csv('data/highlighted_genes.csv', header=None, names=['gene', 'cluster'])
markers['gene'] = markers['gene'].str.upper()
markers = markers[markers['gene'].isin(adata.var.index)]

markers['title'] = markers['gene'] + '+: ' + markers['cluster']
markers = markers.set_index('gene')
markers.loc['PTPRC']['title'] = 'PTPRC (CD45)+: Immune Cells'

markers.loc['leiden'] = ['Leiden', 'Clusters']
#%%
addl_genes = pd.read_csv(os.path.join(sc.settings.figdir, 'additional_genes.csv'), header=None)
addl_genes['title'] = addl_genes[0]
addl_genes = addl_genes.set_index(0)
markers = markers.append(addl_genes)

#%%
sc.pl.tsne(adata, color='leiden',
       title='Leiden',
       color_map='plasma',
       size=25,
       save='_' + 'leiden' + '_all.pdf',
       show=False)

#%%
for i, g in markers.iterrows():
    sc.pl.tsne(adata, color=i,
           title=g['title'],
           color_map='plasma',
           size=25,
           save='_' + i + '_all.pdf',
           show=False)
    
#%%
adata.obs['TMPRS11F'] = 0
proteases = ['TMPRSS2', 'TMPRSS13', 'TMPRSS11D', 'TMPRSS11E']
for i in proteases:
    sc.pl.tsne(adata, color=i,
           title=i,
           color_map='plasma',
           size=25,
           save='_' + i + '_all.pdf',
           show=False)
sc.pl.tsne(adata, color='TMPRS11F',
   title='TMPRS11F',
   color_map='plasma',
   size=25,
   save='_' + 'TMPRS11F' + '_all.pdf',
   vmin=0,
   vmax=1,
   show=False)

#%%
any_tmprss_meta = ((adata.to_df()[proteases] > 0).sum(axis=1) > 0)
any_tmprss_meta = any_tmprss_meta.map({True: 'Present', False: 'Absent'}).astype('category').cat.reorder_categories(['Present', 'Absent'])
adata.obs['any_tmprss_metagene'] = any_tmprss_meta
fig = sc.pl.tsne(adata, color='any_tmprss_metagene',
   title='TTSP-expressing cells',   sort_order=True,
   palette=['#1fb5f0', '#DCDCDC'],
   size=25,
   groups=['Present'],
   save='_' + 'tmprss_metagene' + '_all.png',
   return_fig = True,
   show=False)


legend_elements = [plt.Line2D([0], [0], marker='o', color='#1fb5f0', label='TTSP-positive', linestyle='None'),
                   plt.Line2D([0],[0], marker='o', color='#DCDCDC', label='TTSP-negative', markerfacecolor='#DCDCDC', linestyle='None')]
fig.axes[0].legend(handles=legend_elements, frameon=False,
                loc='center left',
                bbox_to_anchor=(1, 0.5),)
plt.show()
fig.savefig(os.path.join(sc.settings.figdir, 'tsne_tmprss_metagene_all-ttsp.pdf'), bbox_inches='tight')

#%%
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', n_genes=20)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).to_csv(os.path.join(sc.settings.figdir, 'rank_groups.csv'))


#%%

count_fractions = (adata.to_df()[proteases + ['ACE2']] > 0).sum() / adata.shape[0]
ace2_tmprss2 = ((adata.to_df()[['TMPRSS2', 'ACE2']] > 0).sum(axis=1) >= 2).sum() / adata.shape[0]
count_fractions['ACE2+TMPRSS2'] = ace2_tmprss2

other_proteases = adata.to_df()[proteases[1:]]

ace2_other_protease = (other_proteases > 0).sum(axis=1).rename('other_tpmrss')
count_fractions['ACE2+other TMPRSS'] = ((adata.to_df()[['ACE2']].join(ace2_other_protease) > 0).sum(axis=1) >= 2).sum() / adata.shape[0]

count_fractions.to_csv(os.path.join(sc.settings.figdir, 'cell_count_fractions.csv'))

#%%

cluster_labels = pd.read_csv(os.path.join(sc.settings.figdir, 'key.csv'), header=None, names=['cluster', 'label'], dtype='str')
cluster_labels.loc[cluster_labels['label'].isna(),'label'] = cluster_labels.loc[cluster_labels['label'].isna(), 'cluster']
cluster_labels_dict = cluster_labels.set_index('cluster')['label'].to_dict()
adata.obs['Cell Types'] = adata.obs.leiden.map(cluster_labels_dict).astype('category')

# change secretory cell color
adata.uns['Cell Types_colors'][15] = '#9370DB' #smg serous
adata.uns['Cell Types_colors'][16] = '#FF4A46' #secretory
adata.uns['Cell Types_colors'][12] = '#809693' #nasal suprabasal
adata.uns['Cell Types_colors'][3] = '#004D43' #deuterosomal

ax = sc.pl.tsne(adata, color=['Cell Types'],
           size=25,
           title='Human Airway',
           save='_' + 'labeled_clusters.pdf')
adata.write(os.path.join(sc.settings.figdir, 'adata.h5ad'))

#%%
trackplot_markers = pd.read_csv(os.path.join(sc.settings.figdir, 'track_plot.csv'), header=None)
trackplot_markers = trackplot_markers.drop(4, axis=0)

axes_list = sc.pl.tracksplot(adata, var_names=trackplot_markers[0].values, groupby='Cell Types', figsize=(26,10), show=False)
[i.yaxis.set_ticks([]) for i in axes_list]
ax = plt.gca()

# shift fibroblasts and plasma cells so they have non-overlapping labels
fig = plt.gcf()
x_labels_axis = fig.get_axes()[-2]
xticks = [i for i in x_labels_axis.get_xticks()]
xticks[5] = xticks[5] - 400
xticks[6] = xticks[6] + 50
xticks[12] = xticks[12] - 300
xticks[14] = xticks[14] + 300
x_labels_axis.set_xticks(xticks)

ax.set_xlabel('')


plt.savefig(os.path.join(sc.settings.figdir, 'tracksplot_cell_types.pdf'), bbox_inches='tight', format='pdf')

#%%
adata = sc.read_h5ad(os.path.join(sc.settings.figdir, 'adata.h5ad'))
#%%

sc.tl.rank_genes_groups(adata, 'Cell Types', method='t-test', n_genes=20)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).to_csv(os.path.join(sc.settings.figdir, 'labeled_clusters_rank_groups.csv'))
adata.write(os.path.join(sc.settings.figdir, 'adata.h5ad'))

#%% Nov 2020 revisions

sc.pl.tsne(adata, color='EPCAM', size=25, color_map='plasma', save='_epcam.pdf')
           
epcam = adata.obs[['Cell Types']].join(adata.to_df()[['EPCAM']])
epcam_rank = epcam.groupby('Cell Types').apply(np.average).sort_values(ascending=False)
epcam_rank.to_csv(os.path.join(sc.settings.figdir, 'epcam_ranks.csv'))
#%%
epithelial_cell_types = list(epcam_rank[epcam_rank > .3].index)
adata.obs['epithelial'] = adata.obs['Cell Types'].isin(epithelial_cell_types)
adata.obs['epithelial'] = adata.obs['epithelial'].where(adata.obs['epithelial'], 'Non-epithelial')
adata.obs['epithelial'] = adata.obs['epithelial'].where(adata.obs['epithelial'] == 'Non-epithelial', 'Epithelial').astype(pd.CategoricalDtype())
adata.uns['epithelial_colors'] = ['#f76f68', '#bdbdbd']
fig, ax = plt.subplots()
sc.pl.tsne(adata, color='epithelial', size=25, ax=ax, title='Epithelial cells')
plt.savefig(os.path.join(sc.settings.figdir, 'epithelial_highlighted.pdf'), bbox_inches='tight') 
#nasal and above highlighted
#highlight within the epithlelial the proteases

#%%
adata.uns['any_tmprss_metagene_colors'] = ['#1fb5f0', '#bdbdbd']
sc.pl.tsne(adata, color='any_tmprss_metagene', size=25, title='TTSP-expressing cells', groups=['Present', 'Absent'])
#%%
adata.obs['epith+protease'] = adata.obs['epithelial'].copy().astype('str')
adata.obs.loc[(adata.obs['epithelial'] == 'Epithelial') & (adata.obs['any_tmprss_metagene'] == 'Present'), 'epith+protease'] = 'Epithelial and TTSP-expressing'
adata.obs.loc[(adata.obs['epithelial'] == 'Epithelial') & (adata.obs['any_tmprss_metagene'] == 'Absent'), 'epith+protease'] = 'Epithelial and TTSP-negative'
adata.obs.loc[(adata.obs['epithelial'] == 'Non-epithelial') & (adata.obs['any_tmprss_metagene'] == 'Present'), 'epith+protease'] = 'Non-epithelial and TTSP-expressing'
adata.obs['epith+protease'] = adata.obs['epith+protease'].astype(pd.CategoricalDtype(['Non-epithelial', 'Non-epithelial and TTSP-expressing', 'Epithelial and TTSP-negative', 'Epithelial and TTSP-expressing'], ordered=True))
adata.uns['epith+protease_colors'] = [ '#bdbdbd', '#000000', '#f76f68','#1fb5f0']
fig, ax = plt.subplots()
sc.pl.tsne(adata, color='epith+protease', size=25, ax=ax, save='_epith_ttsp_status')



