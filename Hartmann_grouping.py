# -*- coding: utf-8 -*-
"""
Created on Wed May 13 16:19:33 2020

@author: kosta
"""

from Exe_SQanalysis import protein_df, rescaled_df
import pandas as pd
from scipy import stats
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats import multitest

def lineplot2(rescaled_df, x_target, y_target, hue, fname):
    sns.set(style="ticks")
    colors = sns.color_palette('muted')
    proteins = rescaled_df.protein.unique()
    fig, ax = plt.subplots(3,3,figsize = (12,12))
    
    x = 0
    for i in range(0, len(proteins), 3):
        
        y = 0
        sns.lineplot(x=x_target, y=y_target, hue=hue, data=rescaled_df[rescaled_df.protein == proteins[i]],
                     color = colors[i], ax = ax[x,y])
        
        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')
        ax[x,y].set_title(proteins[i])
        
        if y == 0:
            ax[x,y].set_ylabel('total_area_protein', fontsize = 12)
        if x == 2 and y < 3:
            ax[x,y].set_xlabel('Sampling timepoint')
        
        y += 1
        sns.lineplot(x=x_target, y=y_target, hue=hue, data=rescaled_df[rescaled_df.protein == proteins[i+1]],
                     color = colors[i], ax = ax[x,y])
    
        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')    
        ax[x,y].set_title(proteins[i + 1])
        
        if x == 2 and y < 3:
            ax[x,y].set_xlabel('Sampling timepoint')
        
        y += 1
        sns.lineplot(x=x_target, y=y_target, hue=hue, data=rescaled_df[rescaled_df.protein == proteins[i+2]], 
                color = colors[i], ax = ax[x,y])
        
        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')    
        ax[x,y].set_title(proteins[i + 2])
        
        if x == 2 and y < 3:
            ax[x,y].set_xlabel('Sampling timepoint')
        x += 1
        
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    fig.suptitle(f'Proteins, impaired healing')
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(f'{fname}.png', dpi = 300, format = 'png')
    plt.savefig(f'{fname}.svg', format = 'svg')
    plt.close()

    
def boxplot2(df, x_target, y_target, hue, fname):
    sns.set(style="ticks")
    colors = sns.color_palette('muted')
    proteins = df.protein.unique()
    fig, ax = plt.subplots(3,3,figsize = (17,12))
    c = 0
    
    for x in range(3):
        for y in range(3):
            protein = proteins[c]
            data= df[df.protein == protein]
            
            sns.boxplot(x=x_target, y=y_target,
             hue=hue, data=data, color = colors[c], ax = ax[x,y])
            
            ax[x,y].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            
            if x == 2:
                ax[x,y].set_xlabel(x_target, fontsize = 12)
            if y == 2:
                ax[x,y].set_ylabel(y_target, fontsize = 12)
            
            ax[x,y].set_title(protein)
            c += 1
            
    fig.suptitle(f'Proteins, impaired healing')
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(f'{fname}.png', dpi = 300, format = 'png')
    plt.savefig(f'{fname}.svg', format = 'svg')
    plt.close()

    
'''Categorization on mode'''
 
#read and cleanup
# =============================================================================
# df_annot = pd.read_excel('Hartmann_annot.xlsx')
# df_annot = df_annot.iloc[1:,:]
# 
# df_annot = df_annot[df_annot.VISIT != 'Screening']
# df_annot['timepoint'] = df_annot.Visit
# df_annot['patient'] = df_annot.apply(lambda row: str(row.SITE) + '-' + str(row.PATIENT), axis=1)
# df_annot = df_annot[df_annot.patient.isin(rescaled_df.patient.unique())]
# df_annot = df_annot[['patient', 'timepoint', 'Mean', df_annot.columns[-4]]]
# 
# df_annot = df_annot[~((df_annot.timepoint==7) & (df_annot.patient=='03-007'))]
# df_annot = df_annot[~((df_annot.timepoint==7) & (df_annot.patient=='08-007'))]
# df_annot = df_annot[~((df_annot.timepoint==1) & (df_annot.patient=='02A-020'))]
# df_annot = df_annot[~((df_annot.timepoint==2) & (df_annot.patient=='02A-020'))]
# 
# #change df_annot to reflect visits
# df_annot['timepoint'] = df_annot.apply(lambda row: row.timepoint + 1, axis=1)
# df_annot = df_annot[((~df_annot.patient.str.contains('B')) & (df_annot.timepoint <=7)) | ((df_annot.patient.str.contains('B')) & (df_annot.timepoint <=8))]
# 
# #more cleanup and merge
# df_annot.timepoint = df_annot.timepoint.astype(int)
# df_annot = df_annot.rename(columns={'Mean':'mean',df_annot.columns[-1]:'trajectory'})
# rescaled_df = rescaled_df[rescaled_df.timepoint > 1]
# fdf = rescaled_df.merge(df_annot, on=['timepoint', 'patient'], how='right')
# fdf = fdf.sort_values('patient')
# 
# fdf['outcome'] = 0
# for p in fdf.patient.unique():
#     arr = np.array(fdf[fdf.patient==p].trajectory)
#     idx = fdf[fdf.patient==p].index
#     m = stats.mode(arr)[0][0]
#     if m == 1:
#         fdf.loc[idx,'outcome'] = 'Improved'
#     elif m == 0:
#         fdf.loc[idx,'outcome'] = 'Neutral'
#     elif m == -1:
#         fdf.loc[idx,'outcome'] = 'Worsened'
# 
# 
# #log transformation
# log_df = fdf.copy()
# log_df['log2_area'] = np.log2(log_df['total_area_protein'])
# log_df.to_excel('Dataframe_annot.xlsx', index=False)
# log_df = log_df[log_df.condition=='Impaired']
# with pd.option_context('mode.use_inf_as_na', True):
#     log_df = log_df.dropna(subset=['log2_area'], how='all')
# 
# =============================================================================
# =============================================================================
# lineplot2(log_df, 'timepoint', 'total_area_protein', 'outcome', 'Modal_lineplot_annot')
# boxplot2(log_df, 'timepoint', 'total_area_protein', 'outcome', 'Modal_boxplot_annot')
# 
# g = sns.catplot(x='timepoint', y='total_area_protein', col='protein', hue='outcome',
#             col_wrap = 3, data=log_df, kind='bar', sharey=False)
# g.fig.subplots_adjust(top=0.9)
# g.fig.suptitle('Proteins, impaired healing')
# g.savefig('Modal_barplot_annot.png', dpi = 300, format = 'png')
# g.savefig('Modal_barplot_annot.svg', format = 'svg')
# 
# =============================================================================


'''Categorization based on paper'''

# =============================================================================
# df_annot = pd.read_excel('Hartmann_annot.xlsx')
# df_annot = df_annot.iloc[1:,:]
# df_annot = df_annot[df_annot.VISIT != 'Screening']
# df_annot['timepoint'] = df_annot.Visit
# df_annot['patient'] = df_annot.apply(lambda row: str(row.SITE) + '-' + str(row.PATIENT), axis=1)
# df_annot = df_annot[df_annot.patient.isin(rescaled_df.patient.unique())]
# df_annot = df_annot[['patient', 'timepoint', 'Mean']]
# df_annot = df_annot[~((df_annot.timepoint==2) & (df_annot.patient=='02A-020'))]
# with pd.option_context('mode.use_inf_as_na', True):
#     df_annot = df_annot.dropna(subset=['Mean'], how='all')
# df_annot.timepoint = df_annot.timepoint.astype(int)
# df_annot = df_annot.sort_values(['patient', 'timepoint'])
# 
# new_annot = df_annot.copy()
# new_annot['Trajectory'] = 0
# for p in new_annot.patient.unique():
#     temp = new_annot[new_annot.patient==p]
#     idx = temp.index
#     for i in range(1,len(idx) -1):
#         if temp.iloc[i-1]['Mean'] > temp.iloc[i]['Mean']and temp.iloc[i+1]['Mean'] < temp.iloc[i]['Mean']:
#             new_annot.loc[idx[i],'Trajectory'] = 'Improving'
#         elif temp.iloc[i-1]['Mean'] > temp.iloc[i]['Mean'] and temp.iloc[i+1]['Mean'] > temp.iloc[i]['Mean']:
#             new_annot.loc[idx[i],'Trajectory'] = 'Intermediate'
#         elif temp.iloc[i-1]['Mean'] < temp.iloc[i]['Mean'] and temp.iloc[i+1]['Mean'] < temp.iloc[i]['Mean']:
#             new_annot.loc[idx[i],'Trajectory'] = 'Intermediate'
#         elif temp.iloc[i-1]['Mean'] < temp.iloc[i]['Mean']and temp.iloc[i+1]['Mean'] > temp.iloc[i]['Mean']:
#             new_annot.loc[idx[i],'Trajectory'] = 'Non healing'
#         else:
#             new_annot.loc[idx[i],'Trajectory'] = 'Error'
# 
# new_annot = new_annot[new_annot.Trajectory != 0]
# fdf = rescaled_df.merge(new_annot, on=['timepoint', 'patient'], how='right')
# fdf = fdf.sort_values('patient')
# fdf = fdf[fdf.condition=='Impaired']
# fdf = fdf.sort_values(['patient', 'timepoint'])
# =============================================================================

# =============================================================================
# lineplot2(fdf, 'timepoint', 'total_area_protein', 'Trajectory', 'Group_lineplot_annot' )
# boxplot2(fdf, 'timepoint', 'total_area_protein', 'Trajectory', 'Group_boxplot_annot')
# 
# col = sns.color_palette('muted')
# b = col.pop(0)
# col.insert(2, b)
# g = sns.catplot(x='timepoint', y='total_area_protein', col='protein', hue='Trajectory',
#             col_wrap = 3, data=fdf, kind='bar', sharey=False, palette=col)
# g.fig.subplots_adjust(top=0.9)
# g.fig.suptitle('Proteins, impaired healing')
# g.savefig('Group_barplot_annot.png', dpi = 300, format = 'png')
# g.savefig('Group_barplot_annot.svg', format = 'svg')
# 
# =============================================================================

#All timepoints grouped
# =============================================================================
# g = sns.catplot(x='Trajectory', y='total_area_protein', col='protein',
#             col_wrap = 3, data=fdf, kind='box', sharey=False, palette=col)
# g.fig.subplots_adjust(top=0.9)
# g.fig.suptitle('Proteins, impaired healing')
# g.savefig('Mesh_boxplot_annot.png', dpi = 300, format = 'png')
# g.savefig('Mesh_boxplot_annot.svg', format = 'svg')
# =============================================================================


'''Categorization based on mean'''

# =============================================================================
# df_annot = pd.read_excel('Hartmann_annot.xlsx')
# df_annot = df_annot.iloc[1:,:]
# df_annot = df_annot[df_annot.VISIT != 'Screening']
# df_annot['timepoint'] = df_annot.Visit
# df_annot['patient'] = df_annot.apply(lambda row: str(row.SITE) + '-' + str(row.PATIENT), axis=1)
# df_annot = df_annot[df_annot.patient.isin(rescaled_df.patient.unique())]
# df_annot = df_annot[['patient', 'timepoint', 'Mean']]
# df_annot = df_annot[~((df_annot.timepoint==2) & (df_annot.patient=='02A-020'))]
# with pd.option_context('mode.use_inf_as_na', True):
#     df_annot = df_annot.dropna(subset=['Mean'], how='all')
# df_annot.timepoint = df_annot.timepoint.astype(int)
# df_annot = df_annot.sort_values(['patient', 'timepoint'])
# 
# new_annot = df_annot.copy()
# new_annot['Trajectory'] = 0
# for p in new_annot.patient.unique():
#     temp = new_annot[new_annot.patient==p]
#     temp = temp.sort_values('timepoint')
#     idx = temp.index
#     valst = temp.loc[idx[0], 'Mean']
#     valend = temp.loc[idx[-1], 'Mean']
#     if valst < valend:
#         new_annot.loc[idx,'Trajectory'] = 'Chronic'
#     elif valend <= valst:
#         if valend < valst/2:
#             new_annot.loc[idx,'Trajectory'] = 'Improved'
#         elif valend > valst/2:
#             new_annot.loc[idx,'Trajectory'] = 'Neutral'
#     
# new_annot = new_annot[new_annot.Trajectory != 0]
# fdf = rescaled_df.merge(new_annot, on=['timepoint', 'patient'], how='right')
# fdf = fdf.sort_values('patient')
# fdf = fdf[fdf.condition=='Impaired']
# fdf = fdf.sort_values(['patient', 'timepoint'])
# =============================================================================

# =============================================================================
# lineplot2(fdf, 'timepoint', 'total_area_protein', 'Trajectory', 'Mean_lineplot_annot' )
# boxplot2(fdf, 'timepoint', 'total_area_protein', 'Trajectory', 'Mean_boxplot_annot')
# 
# col = sns.color_palette('muted')
# b = col.pop(0)
# col.insert(2, b)
# g = sns.catplot(x='timepoint', y='total_area_protein', col='protein', hue='Trajectory',
#             col_wrap = 3, data=fdf, kind='bar', sharey=False, palette=col)
# g.fig.subplots_adjust(top=0.9)
# g.fig.suptitle('Proteins, impaired healing')
# g.savefig('Mean_barplot_annot.png', dpi = 300, format = 'png')
# g.savefig('Mean_barplot_annot.svg', format = 'svg')
# =============================================================================



''' END OF ANALYSES '''

# =============================================================================
# #statistics for groups
# ind_l = []
# time_l = []
# p_l = []
# fh = open(f'Group_Kruskal_allgroups.txt', 'w')
# fhraw = open(f'Group_Kruskal_allgroups_raw.txt', 'w')
# for i in fdf.protein.unique():
#     for t in fdf.timepoint.unique():
#         test = fdf[(fdf.protein== i) & (fdf.timepoint == t)]
#         h, p = stats.kruskal(*[group["total_area_protein"].values for name, group in test.groupby("Trajectory")])
#         p_l.append(p)
#         ind_l.append(i)
#         time_l.append(t)
# reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.05, method = 'bonferroni')    
# for i in range(len(ind_l)):    
#     fh.write(f'{ind_l[i]}-{time_l[i]}\t{p_lnew[i]}\n')
#     fhraw.write(f'{ind_l[i]}-{time_l[i]}\t{p_l[i]}\n')
#     
# fh.close()
# fhraw.close()
# 
# #Statistics ranksums for condition comparing timepoints
# ind_l = []
# p_l = []
# t_l = [] 
# t2_l = []
# 
# for k in rescaled_df.protein.unique():
#     s = 2
#     c = 3
#     
#     while s <= 6:
#         end = 7
#         for i in range(c,end):
#             sample1 = fdf[(fdf.protein == k) & (fdf.Trajectory == 'Improving') 
#                                   & (fdf.timepoint == s)]
#             sample2 = fdf[(fdf.protein == k) & (fdf.Trajectory == 'Non healing') 
#                                   & (fdf.timepoint == i)]
#             u, p = stats.ranksums(sample1.total_area_protein, sample2.total_area_protein)
#             ind_l.append(k)
#             t_l.append(s)
#             t2_l.append(i)
#             p_l.append(p)
#         
#         c += 1
#         s += 1
#     
# reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.05, method = 'bonferroni')
# 
# fh = open(f'Group_ranksum_2group_prot.txt', 'w')
# 
# for i in range(len(reject)):
#     if reject[i]:
#         fh.write(f'{ind_l[i]}\t{t_l[i]}+{t2_l[i]} - {p_lnew[i]}\n')
# fh.close()
# 
# with open(f'Group_ranksum_2group_prot_raw.txt','w') as fhraw:
#     for i,p in enumerate(p_l):
#         fhraw.write(f'{ind_l[i]}\t{t_l[i]}+{t2_l[i]} - {p}\n')
# 
# =============================================================================

#statistics for mode analysis
# =============================================================================
# ind_l = []
# time_l = []
# p_l = []
# fh = open(f'Modal_Kruskal_allgroups.txt', 'w')
# fhraw = open(f'Modal_Kruskal_allgroups_raw.txt', 'w')
# for i in log_df.protein.unique():
#     for t in log_df.timepoint.unique():
#         test = log_df[(log_df.protein== i) & (log_df.timepoint == t)]
#         h, p = stats.kruskal(*[group["total_area_protein"].values for name, group in test.groupby("outcome")])
#         p_l.append(p)
#         ind_l.append(i)
#         time_l.append(t)
# reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.05, method = 'bonferroni')    
# for i in range(len(ind_l)):    
#     fh.write(f'{ind_l[i]}-{time_l[i]}\t{p_lnew[i]}\n')
#     fhraw.write(f'{ind_l[i]}-{time_l[i]}\t{p_l[i]}\n')
#     
# fh.close()
# fhraw.close()
# 
# #Statistics ranksums for condition comparing timepoints
# ind_l = []
# p_l = []
# t_l = [] 
# t2_l = []
# 
# for k in rescaled_df.protein.unique():
#     s = 2
#     c = 3
#     
#     while s <= 6:
#         end = 7
#         for i in range(c,end):
#             sample1 = log_df[(log_df.protein == k) & (log_df.outcome == 'Improved') 
#                                   & (log_df.timepoint == s)]
#             sample2 = log_df[(log_df.protein == k) & (log_df.outcome == 'Worsened') 
#                                   & (log_df.timepoint == i)]
#             u, p = stats.ranksums(sample1.total_area_protein, sample2.total_area_protein)
#             ind_l.append(k)
#             t_l.append(s)
#             t2_l.append(i)
#             p_l.append(p)
#         
#         c += 1
#         s += 1
#     
# reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.05, method = 'bonferroni')
# 
# fh = open(f'Modal_ranksum_2group_prot.txt', 'w')
# 
# for i in range(len(reject)):
#     if reject[i]:
#         fh.write(f'{ind_l[i]}\t{t_l[i]}+{t2_l[i]} - {p_lnew[i]}\n')
# fh.close()
# 
# with open(f'Modal_ranksum_2group_prot_raw.txt','w') as fhraw:
#     for i,p in enumerate(p_l):
#         fhraw.write(f'{ind_l[i]}\t{t_l[i]}+{t2_l[i]} - {p}\n')
# 
# =============================================================================


#statistics for mean analysis
# =============================================================================
# ind_l = []
# time_l = []
# p_l = []
# fh = open(f'Mean_Kruskal_allgroups.txt', 'w')
# fhraw = open(f'Mean_Kruskal_allgroups_raw.txt', 'w')
# for i in fdf.protein.unique():
#     for t in fdf.timepoint.unique():
#         test = fdf[(fdf.protein== i) & (fdf.timepoint == t)]
#         h, p = stats.kruskal(*[group["total_area_protein"].values for name, group in test.groupby("Trajectory")])
#         p_l.append(p)
#         ind_l.append(i)
#         time_l.append(t)
# reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.05, method = 'bonferroni')    
# for i in range(len(ind_l)):    
#     fh.write(f'{ind_l[i]}-{time_l[i]}\t{p_lnew[i]}\n')
#     fhraw.write(f'{ind_l[i]}-{time_l[i]}\t{p_l[i]}\n')
#     
# fh.close()
# fhraw.close()
# 
# #Statistics ranksums for condition comparing timepoints
# ind_l = []
# p_l = []
# t_l = [] 
# t2_l = []
# 
# for k in rescaled_df.protein.unique():
#     s = 2
#     c = 3
#     
#     while s <= 6:
#         end = 7
#         for i in range(c,end):
#             sample1 = fdf[(fdf.protein == k) & (fdf.Trajectory == 'Improved') 
#                                   & (fdf.timepoint == s)]
#             sample2 = fdf[(fdf.protein == k) & (fdf.Trajectory == 'Chronic') 
#                                   & (fdf.timepoint == i)]
#             u, p = stats.ranksums(sample1.total_area_protein, sample2.total_area_protein)
#             ind_l.append(k)
#             t_l.append(s)
#             t2_l.append(i)
#             p_l.append(p)
#         
#         c += 1
#         s += 1
#     
# reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.05, method = 'bonferroni')
# 
# fh = open(f'Mean_ranksum_2group_prot.txt', 'w')
# 
# for i in range(len(reject)):
#     if reject[i]:
#         fh.write(f'{ind_l[i]}\t{t_l[i]}+{t2_l[i]} - {p_lnew[i]}\n')
# fh.close()
# 
# with open(f'Mean_ranksum_2group_prot_raw.txt','w') as fhraw:
#     for i,p in enumerate(p_l):
#         fhraw.write(f'{ind_l[i]}\t{t_l[i]}+{t2_l[i]} - {p}\n')
# 
# =============================================================================

#statistics for mesh analysis
# =============================================================================
# ind_l = []
# p_l = []
# fh = open(f'Mesh_Kruskal_allgroups.txt', 'w')
# fhraw = open(f'Mesh_Kruskal_allgroups_raw.txt', 'w')
# for i in fdf.protein.unique():
#     test = fdf[fdf.protein== i]
#     h, p = stats.kruskal(*[group["total_area_protein"].values for name, group in test.groupby("Trajectory")])
#     p_l.append(p)
#     ind_l.append(i)
# reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.05, method = 'bonferroni')    
# for i in range(len(ind_l)):    
#     fh.write(f'{ind_l[i]}\t{p_lnew[i]}\n')
#     fhraw.write(f'{ind_l[i]}\t{p_l[i]}\n')
#     
# fh.close()
# fhraw.close()
# 
# 
# ind_l = []
# time_l = []
# p_l = []
# fh = open(f'Mesh_ranksum_2groups.txt', 'w')
# fhraw = open(f'Mesh_ranksum_2groups_raw.txt', 'w')
# for i in fdf.protein.unique():
#     sample1 = fdf[(fdf.protein == i) & (fdf.Trajectory == 'Improving')]
#     sample2 = fdf[(fdf.protein == i) & (fdf.Trajectory == 'Non healing')]
#     u, p = stats.ranksums(sample1.total_area_protein, sample2.total_area_protein)
#     p_l.append(p)
#     ind_l.append(i)
# reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.05, method = 'bonferroni')    
# for i in range(len(ind_l)):    
#     fh.write(f'{ind_l[i]}\t{p_lnew[i]}\n')
#     fhraw.write(f'{ind_l[i]}\t{p_l[i]}\n')
#     
# fh.close()
# fhraw.close()
# =============================================================================
