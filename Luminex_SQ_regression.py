# -*- coding: utf-8 -*-
"""
Created on Mon May 11 18:00:34 2020

@author: kosta
"""

import pandas as pd
from Exe_SQanalysis import protein_df
import seaborn as sns
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

def rescale(master_df, col):
    
    scaled_df = pd.DataFrame(columns = master_df.columns)
    proteins = master_df.protein.unique()
    
    for prot in proteins:
        temp = master_df[master_df.protein == prot]
        m = max(temp[col])
        temp[col] = np.array(temp[col])/m
        scaled_df = scaled_df.append(temp)
    
    return scaled_df

#########SQ vs Luminex############
patients = ['02A-001', '02A-002', '02A-005', '02A-010', '02A-013', '02A-021']
df_sq = protein_df[protein_df.patient.isin(patients)]
df_sq = df_sq[df_sq.protein != 'ELNE']

df_lum = pd.read_excel('Luminex final data.xlsx')
df_lum = df_lum.drop(0)
df_lum = df_lum.drop('Unnamed: 1', axis=1)

timep = []
for i in range(0,len(df_lum), 7):
    ind = 1
    for k in range(7):
        timep.append(ind)
        ind +=1

df_lum.insert(1, 'timepoint', timep)
df_lum = df_lum.rename(columns={'Values are in pg/ml':'patient',
                                'TNF-a':'TNFA',
                                'IL-1b':'IL1B',
                                'Collagen':'CO1A1',
                                'Fibronectin':'FINC',
                                'S100A8':'S10A8',
                                'S100A9':'S10A9'})

df_lum = pd.melt(df_lum, id_vars=['patient', 'timepoint'])
df_lum = df_lum.rename(columns={'variable':'protein'})
df_lum = df_lum.fillna(0)

df = df_sq.merge(df_lum, on=['patient','timepoint','protein'])
df = df.rename(columns={'value':'luminex','total_area_protein':'surequant'})

#df = rescale(df,'luminex')
#df = rescale(df,'surequant')
df = df.sort_values('sample')

df_mel = df.melt(id_vars=['patient','sample','condition','type','protein','timepoint'],
                 value_vars=['luminex', 'surequant'], var_name='method')


colors = sns.color_palette('muted')
models = []
proteins = list(df.protein.unique())
fig, ax = plt.subplots(4,2,figsize = (8,10))
proteins.pop(4)
c=0
for x in range(4):
    for y in range(2):
        
        if c==7:
            continue
        protein = proteins[c]
        df_prot = df[df.protein == protein]
        lum_vals = list(df_prot.luminex)
        sq_vals = list(df_prot.surequant)
        
        peru = 90
        perl = 0
        perulum = np.percentile(lum_vals, peru)
        perusq = np.percentile(sq_vals, peru)
        perllum = np.percentile(lum_vals, perl)
        perlsq = np.percentile(sq_vals, perl)
        for i in range(len(lum_vals)-1,0,-1):
            if lum_vals[i] == 0:
                lum_vals.pop(i)
                sq_vals.pop(i)
                continue
            if lum_vals[i] > perulum or lum_vals[i] < perllum:
                lum_vals.pop(i)
                sq_vals.pop(i)
            elif sq_vals[i] > perusq or sq_vals[i] < perlsq:
                lum_vals.pop(i)
                sq_vals.pop(i)
        
        lv,sv = zip(*sorted(zip(lum_vals, sq_vals)))
        
        X = np.array(sv).reshape(-1, 1)
        Y = np.array(lv)
        reg = LinearRegression().fit(X, Y)
        
        predx = np.arange(min(sv), max(sv), (max(sv)-min(sv))/100).reshape(-1,1)
        predy = reg.predict(predx)
        
        sns.set_style('ticks')
        
        sns.scatterplot(sv, lv, color=colors[c], ax=ax[x][y])
        sns.lineplot(predx.reshape(-1),predy, color=colors[c], ax=ax[x][y])
        ax[x][y].set_xlabel('Surequant')
        ax[x][y].set_ylabel('Luminex')
        ax[x][y].set_title(protein)
        ax[x][y].text(max(sv)*0.82, min(lv)*1.2, s=f'R^2 = {round(reg.score(X,Y),3)}')
        sns.despine()
        c += 1

fig.delaxes(ax[3][1])
plt.tight_layout()
plt.savefig('Reg_sqlum.png', dpi=400, format='png')