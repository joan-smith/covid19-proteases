#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 19:26:43 2020

@author: Joan Smith
"""

#%%
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import os

gtex = pd.read_csv("data/raw-data/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", sep="\t", header=2, index_col=1)
gtex = gtex.drop('Name', axis=1).astype(float)
attributes = pd.read_csv("data/raw-data/GTEX/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep='\t', index_col=0, dtype=None)

#%%
def set_labels(ax, labels):
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)

#%% Plot All Tissue
proteases = ['TMPRSS2', 'TMPRSS13', 'TMPRSS11D', 'TMPRSS11E', 'TMPRSS11F', 'ACE2']

for p in proteases: 
    protease_tissue = gtex.loc[[p]].T.join(attributes)
    protease_tissue[p] = np.log2(protease_tissue[p] + 1)
    
    fig, ax = plt.subplots(figsize=(12,6))
    plt.title(p + ' by tissue')
    
    order = protease_tissue.groupby('SMTS')[p].apply(np.mean).sort_values()
    print(order)
    g = {i[0]: i[1][p].dropna() for i in protease_tissue.groupby('SMTS')}
    ordered_g= [(k, g[k]) for k in order.index]
    labels, split = zip(*ordered_g)
    ax.violinplot(split, showmeans=True)
    set_labels(ax, labels)
    plt.xticks(rotation=90)
    plt.savefig(os.path.join('data/protease-expression', 'gtex_' + p + '_all_tissues.png'), bbox_inches='tight')
    
    pivot = protease_tissue[[p]].join(attributes).reset_index().pivot_table(columns='SMTS', values=p, index='index')
    pivot[order.index].to_csv(
        os.path.join('data/protease-expression', 'gtex_'+p+'_by_tissue.csv'))

