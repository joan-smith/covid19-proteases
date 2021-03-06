adat#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 2020

@author: Joan Smith

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122960
"""
#%%
import pandas as pd
import os

import scanpy as sc
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
#%%

sc.settings.figdir = 'data/protease-expression/cluster-plots/ipf_gse1229560/'
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
adata = sc.read_h5ad(os.path.join(sc.settings.figdir, 'adata.h5ad'))
params = pd.read_csv(os.path.join(sc.settings.figdir, 'params.txt'), index_col=0).to_dict('dict')['0']

#%%

proteases = ['TMPRSS2', 'TMPRSS13', 'TMPRSS11D', 'TMPRSS11E', 'TMPRSS11F']
adata.obs['any_tmprss_metagene'] = any_tmprss_meta

for i in proteases + ['ACE2']:
    sc.pl.tsne(adata, color=i,
           title=i,
           color_map='plasma',
           size=25,
           save='_' + i + '_all.pdf',
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
fig.savefig(os.path.join(sc.settings.figdir, 'tsne_any_tmprss_metagene_all.pdf'), bbox_inches='tight')

#%%
sc.pl.tsne(adata, color='TMPRSS11F',
   title='TMPSRSS11F',
   color_map='plasma',
   size=25,
   save='_' + 'TMPRSS11F' + '_all.pdf',
   vmin=0,
   vmax=1,
   show=False)


#%%

count_fractions = (adata.to_df()[proteases + ['ACE2']] > 0).sum() / adata.shape[0]
ace2_tmprss2 = ((adata.to_df()[['TMPRSS2', 'ACE2']] > 0).sum(axis=1) >= 2).sum() / adata.shape[0]
count_fractions['ACE2+TMPRSS2'] = ace2_tmprss2

other_proteases = adata.to_df()[proteases[1:]]

ace2_other_protease = (other_proteases > 0).sum(axis=1).rename('other_tpmrss')
count_fractions['ACE2+other TMPRSS'] = ((adata.to_df()[['ACE2']].join(ace2_other_protease) > 0).sum(axis=1) >= 2).sum() / adata.shape[0]

count_fractions.to_csv(os.path.join(sc.settings.figdir, 'cell_count_fractions.csv'))

#%%

adata.obs['Cell Types'] = adata.obs['Cell Types'].replace('Goblet/Club cells', 'Secretory cells')
ax = sc.pl.tsne(adata, color=['Cell Types'],
           size=25,
           title='Human Lung',
           save='_' + 'labeled_clusters.pdf')

#%%
trackplot_markers = pd.read_csv(os.path.join(sc.settings.figdir, 'trackplot_human.csv'), header=None)
trackplot_markers = trackplot_markers.drop(4, axis=0)
axes_list = sc.pl.tracksplot(adata, var_names=trackplot_markers[0].values, groupby='Cell Types', figsize=(20,6), show=False)
[i.yaxis.set_ticks([]) for i in axes_list]
ax = plt.gca()
ax.set_xlabel('')
plt.savefig(os.path.join(sc.settings.figdir, 'tracksplot_cell_types.tiff'), bbox_inches='tight', dpi=1000)

#%%
# Nov 2020 Revisions

sc.pl.tsne(adata, color='EPCAM', size=25, color_map='plasma', save='_epcam.pdf')
           
epcam = adata.obs[['Cell Types']].join(adata.to_df()[['EPCAM']])
epcam_rank = epcam.groupby('Cell Types').apply(np.average).sort_values(ascending=False)
epcam_rank.to_csv(os.path.join(sc.settings.figdir, 'epcam_ranks.csv'))
#%%
epithelial_cell_types = list(epcam_rank[epcam_rank > 1].index)
adata.obs['epithelial'] = adata.obs['Cell Types'].isin(epithelial_cell_types)
adata.obs['epithelial'] = adata.obs['epithelial'].where(adata.obs['epithelial'], 'Non-epithelial')
adata.obs['epithelial'] = adata.obs['epithelial'].where(adata.obs['epithelial'] == 'Non-epithelial', 'Epithelial').astype(pd.CategoricalDtype())
adata.uns['epithelial_colors'] = ['#f76f68', '#bdbdbd']
fig, ax = plt.subplots()
sc.pl.tsne(adata, color='epithelial', size=25, ax=ax, title='Epithelial cells')
plt.savefig(os.path.join(sc.settings.figdir, 'epithelial_highlighted.pdf'), bbox_inches='tight') 
