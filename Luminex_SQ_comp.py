# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 10:44:28 2020

@author: kosta
"""

import pandas as pd
from Exe_SQanalysis import protein_df
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns
import numpy as np

def rescale(master_df, col):
    
    scaled_df = pd.DataFrame(columns = master_df.columns)
    proteins = master_df.protein.unique()
    
    for prot in proteins:
        temp = master_df[master_df.protein == prot]
        m = max(temp[col])
        temp[col] = np.array(temp[col])/m
        scaled_df = scaled_df.append(temp)
    
    return scaled_df

def boxswarmls(df):
    
    proteins = df.protein.unique()
    fig, ax = plt.subplots(4,2, figsize = (12,16))
    ind = 0
    for x in range(4):
        for y in range(2):
            pdata = df[df.protein == proteins[ind]]
    
            sns.boxplot(x="timepoint", y="value", hue="method", data=pdata, palette="vlag", ax=ax[x][y])
            #sns.swarmplot(x="timepoint", y="value", hue="method", data=pdata, edgecolors='b', size = 6, linewidth=1, palette="vlag",ax=ax[x][y])
            sns.despine()
            
            ax[x][y].set_title(proteins[ind])
            if x == 0:
                ax[x][y].set_ylabel('Normalized area')
            else:
                ax[x,y].set_ylabel('')
            if y == 1:
                ax[x][y].set_xlabel('Timepoint')
            else:
                ax[x,y].set_xlabel('')
            if x == 0 and y == 1:
                ax[x,y].legend(loc='center left', bbox_to_anchor=(1.02, 1))
            else:
                ax[x,y].get_legend().set_visible(False)
            ind += 1
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig('Comp_luminex_sq.svg',format='svg')
    plt.savefig('Comp_luminex_sq.png',format='png', dpi=300) 

def boxswarmps(df):
    
    proteins = df.protein.unique()
    fig, ax = plt.subplots(4,2, figsize = (12,16))
    ind = 0
    for x in range(4):
        for y in range(2):
            if ind == 7:
                continue
            pdata = df[df.protein == proteins[ind]]
    
            sns.boxplot(x="timepoint", y="value", hue="method", data=pdata, palette="vlag", ax=ax[x][y])
            #sns.swarmplot(x="timepoint", y="value", hue="method", data=pdata, edgecolors='b', size = 6, linewidth=1, palette="vlag",ax=ax[x][y])
            sns.despine()
            
            ax[x][y].set_title(proteins[ind])
            if x == 0:
                ax[x][y].set_ylabel('Normalized area')
            else:
                ax[x,y].set_ylabel('')
            if y == 1:
                ax[x][y].set_xlabel('Timepoint')
            else:
                ax[x,y].set_xlabel('')
            if x == 0 and y == 1:
                ax[x,y].legend(loc='center left', bbox_to_anchor=(1.02, 1))
            else:
                ax[x,y].get_legend().set_visible(False)
            ind += 1
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig('Comp_prm_sq.svg',format='svg')
    plt.savefig('Comp_prm_sq.png',format='png', dpi=300)

             
def boxswarmpls(df):
    
    proteins = df.protein.unique()
    fig, ax = plt.subplots(3,2, figsize = (12,16))
    ind = 0
    for x in range(3):
        for y in range(2):
            pdata = df[df.protein == proteins[ind]]
    
            sns.boxplot(x="timepoint", y="value", hue="method", data=pdata, palette="vlag", ax=ax[x][y])
            #sns.swarmplot(x="timepoint", y="value", hue="method", data=pdata, edgecolors='b', size = 6, linewidth=1, palette="vlag",ax=ax[x][y])
            sns.despine()
            
            ax[x][y].set_title(proteins[ind])
            if x == 0:
                ax[x][y].set_ylabel('Normalized area')
            else:
                ax[x,y].set_ylabel('')
            if y == 1:
                ax[x][y].set_xlabel('Timepoint')
            else:
                ax[x,y].set_xlabel('')
            if x == 0 and y == 1:
                ax[x,y].legend(loc='center left', bbox_to_anchor=(1.02, 1))
            else:
                ax[x,y].get_legend().set_visible(False)
            ind += 1
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig('Comp_all.svg',format='svg')
    plt.savefig('Comp_all.png',format='png', dpi=300)
    
    
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

df = rescale(df,'luminex')
df = rescale(df,'surequant')
df = df.sort_values('sample')

df_mel = df.melt(id_vars=['patient','sample','condition','type','protein','timepoint'],
                 value_vars=['luminex', 'surequant'], var_name='method')

boxswarmls(df_mel)

#######PRM vs SQ#########
df_sq = protein_df[protein_df.patient.isin(patients)]
df_sq = df_sq[df_sq.protein != 'TNFA']
df_sq = df_sq[df_sq.protein != 'IL1B']
df_prm = pd.read_excel('PRM final data.xlsx')
df_prm = df_prm.drop(0)
df_prm = df_prm.drop('Unnamed: 1', axis=1)
df_prm.insert(1, 'timepoint', timep)
df_prm = df_prm.rename(columns={'Unnamed: 0':'patient',
                                'Collagen':'CO1A1',
                                'Fibronectin':'FINC',
                                'S100A8':'S10A8',
                                'S100A9':'S10A9',
                                'Neutrophil Elastase':'ELNE'})
df_prm = pd.melt(df_prm, id_vars=['patient', 'timepoint'])
df_prm = df_prm.rename(columns={'variable':'protein'})
df = df_sq.merge(df_prm, on=['patient','timepoint','protein'])
df = df.rename(columns={'value':'prm','total_area_protein':'surequant'})
df = rescale(df,'prm')
df = rescale(df,'surequant')
df = df.sort_values('sample')
df_mel = df.melt(id_vars=['patient','sample','condition','type','protein','timepoint'],
                 value_vars=['prm', 'surequant'], var_name='method')
df_mel['value'] = df_mel['value'].astype(float)
boxswarmps(df_mel)

#########ALL#########
df_sq = protein_df[protein_df.patient.isin(patients)]
df_sq = df_sq[df_sq.protein != 'TNFA']
df_sq = df_sq[df_sq.protein != 'IL1B']
df_sq = df_sq[df_sq.protein != 'ELNE']
df = df_sq.merge(df_prm, on=['patient','timepoint','protein'])
df = df.rename(columns={'value':'prm','total_area_protein':'surequant'})
df = df.merge(df_lum,  on=['patient','timepoint','protein'])
df = df.rename(columns={'value':'luminex'})
df = rescale(df,'prm')
df = rescale(df,'surequant')
df = rescale(df,'luminex')
df_mel = df.melt(id_vars=['patient','sample','condition','type','protein','timepoint'],
                 value_vars=['prm', 'surequant','luminex'], var_name='method')
df_mel['value'] = df_mel['value'].astype(float)
boxswarmpls(df_mel)