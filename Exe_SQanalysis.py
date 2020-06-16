# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 12:50:59 2020

@author: konka
"""

#import functions
from AnalyzeSureQuant import *

#RESULTS 2

#process data
file_pre = 'ResultSQ2_plate'
dataframe = pd.DataFrame()

for i in range(1,5):    
    filename = file_pre + str(i) + '.csv'
    plate = read_clean_plate(filename)      
    dataframe = dataframe.append(plate, ignore_index = True)

factor = ReadChrom()    
df = ChromNorm(dataframe, factor)
protein_df = get_protein_values2(df)

#remove MEC
#protein_df = protein_df[protein_df.timepoint != 1]

rescaled_df = rescale(protein_df)    
fold_df = fold_change(rescaled_df)

#log transformation
log_df = fold_df.copy()
log_df['total_area_protein'] = np.log2(log_df['total_area_protein'])

# =============================================================================
# #Rename timepoints, remove no8 from acute
# final_df = rescaled_df.copy()
# lst = []
# for row in range(len(final_df)):
#     lst.append(final_df.iloc[row]['timepoint'] -1 )
# final_df.timepoint = lst
# =============================================================================

#final_df = final_df[final_df.timepoint < 7]

# =============================================================================
# #plot
# #remember to switch column from fold_change to total_protein_area in the functions
# lineplot2(log_df)
# plot_all(log_df)
# clustermap(log_df)
# final_box(log_df, 'total_area_protein')
# box2(log_df,'Acute')
# swarmbox(log_df)
# =============================================================================



# =============================================================================
# #NORMALITY PLOTS
# from statsmodels.graphics.gofplots import qqplot
# test = final_df[(final_df.protein == 'MMP9') & (final_df.condition == 'Impaired') & 
#                 (final_df.timepoint == 3)].total_area_protein
# test = test.reset_index(drop=True)
# 
# lst_lag = []
# for i in range(len(test)):
#     if i == len(test) - 1:
#         lst_lag.append(test.iloc[0])
#     else:
#         lst_lag.append(test.iloc[i +1])
#     
# fig, ax = plt.subplots(2,2,figsize=(8,8))
# ax[0,0].plot(test.index, test)
# ax[0,1].scatter(np.array(test), lst_lag, s=10)
# ax[1,0].hist(test) 
# qqplot(test, line='s', ax=ax[1,1])
# plt.show()
# =============================================================================


# =============================================================================
# #Statistics kruskal mann whitney for condition comparing timepoints
# for cond in ['Acute', 'Impaired']:
#     
#     fh = open(f'Kruskal_allprot_{cond}_woMC.txt', 'w')
#     for i in rescaled_df.protein.unique():
#         test = rescaled_df[(rescaled_df.protein == i) & (rescaled_df.condition == cond)]
#         h, p = stats.kruskal(*[group["total_area_protein"].values for name, group in test.groupby("timepoint")])
#         print(i,'\t',p,'\n', file = fh)
#     fh.close()
#     
#     ind_l = []
#     p_l = []
#     t_l = [] 
#     t2_l = []
#     
#     for k in rescaled_df.protein.unique():
#         s = 1
#         c = 2
#         
#         while s <= 7:
#             if cond == 'Acute':
#                 end = 9
#             elif cond == 'Impaired':
#                 end = 8
#             for i in range(c,end):
#                 sample1 = rescaled_df[(rescaled_df.protein == k) & (rescaled_df.condition == cond) 
#                                       & (rescaled_df.timepoint == s)]
#                 sample2 = rescaled_df[(rescaled_df.protein == k) & (rescaled_df.condition == cond) 
#                                       & (rescaled_df.timepoint == i)]
#                 u, p = stats.mannwhitneyu(sample1.total_area_protein, sample2.total_area_protein)
#                 ind_l.append(k)
#                 t_l.append(s)
#                 t2_l.append(i)
#                 p_l.append(p)
#             
#             c += 1
#             s += 1
#         
#     reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.01, method = 'bonferroni')
#     
#     fh = open(f'MannwhitneyU_res_proteinpertimepoint_{cond}_woMC.txt', 'w')
#     
#     for i in range(len(reject)):
#         if reject[i]:
#             fh.write(f'{ind_l[i]}\t{t_l[i]}+{t2_l[i]} - {p_lnew[i]}\n')
#     fh.close()
#     
#     with open(f'MannWhitney_{cond}_raw.txt','w') as fhraw:
#         for i,p in enumerate(p_l):
#             fhraw.write(f'{ind_l[i]}\t{t_l[i]}+{t2_l[i]} - {p}\n')
# =============================================================================


# =============================================================================
# #Statistics kruskal wilcoxon for condition comparing timepoints
# for cond in ['Acute', 'Impaired']:
#     
#     fh = open(f'Kruskal_allprot_{cond}.txt', 'w')
#     for i in rescaled_df.protein.unique():
#         test = rescaled_df[(rescaled_df.protein == i) & (rescaled_df.condition == cond)]
#         h, p = stats.kruskal(*[group["total_area_protein"].values for name, group in test.groupby("timepoint")])
#         print(i,'\t',p,'\n', file = fh)
#     fh.close()
#     
#     ind_l = []
#     p_l = []
#     t_l = [] 
#     t2_l = []
#     
#     for k in rescaled_df.protein.unique():
#         s = 1
#         c = 2
#         
#         while s <= 7:
#             if cond == 'Acute':
#                 end = 9
#             elif cond == 'Impaired':
#                 end = 8
#             for i in range(c,end):
#                 sample1 = rescaled_df[(rescaled_df.protein == k) & (rescaled_df.condition == cond) 
#                                       & (rescaled_df.timepoint == s)]
#                 sample2 = rescaled_df[(rescaled_df.protein == k) & (rescaled_df.condition == cond) 
#                                       & (rescaled_df.timepoint == i)]
#                 u, p = stats.wilcoxon(sample1.total_area_protein, sample2.total_area_protein)
#                 ind_l.append(k)
#                 t_l.append(s)
#                 t2_l.append(i)
#                 p_l.append(p)
#             
#             c += 1
#             s += 1
#         
#     reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.01, method = 'bonferroni')
#     
#     fh = open(f'Wilcoxon_res_proteinpertimepoint_{cond}.txt', 'w')
#     
#     for i in range(len(reject)):
#         if reject[i]:
#             fh.write(f'{ind_l[i]}\t{t_l[i]}+{t2_l[i]} - {p_lnew[i]}\n')
#     fh.close()
#     
#     with open(f'Wilcoxon_{cond}_raw.txt','w') as fhraw:
#         for i,p in enumerate(p_l):
#             fhraw.write(f'{ind_l[i]}\t{t_l[i]}+{t2_l[i]} - {p}\n')
# =============================================================================

# =============================================================================
# protl = []
# pl = []
# 
# for prot in fold_df.protein.unique():
#     dat = fold_df[fold_df.protein == prot]
#     dat = dat[dat.condition == 'Acute']
#     dat = dat.pivot('patient', 'timepoint', 'total_area_protein')
#     data = []
#     for i in dat.columns:
#         data.append(np.array(dat[i]))
#     d1, d2, d3, d4, d5, d6, d7, d8 = [i for i in data]
#     res = stats.friedmanchisquare(d1, d2, d3, d4, d5, d6, d7)
#     print(prot,res[1])
#     protl.append(prot)
#     pl.append(res[1])
# 
# reject, p_lnew, asid, abon = multitest.multipletests(pl,alpha = 0.01, method = 'bonferroni')
# fh = open('Friedman_acute.txt','w')
# for i,p in enumerate(protl):
#     fh.write(p+'\t'+str(p_lnew[i])+'\n')
# fh.close()
# =============================================================================

        
'''Obscure'''
# =============================================================================
# #boxplots all
# sns.set(context='notebook', style='whitegrid', palette = 'deep', font= 'Helvetica')
# #final_df = final_df[final_df.timepoint < 8]
# 
# for prot in final_df.protein.unique():
#     sns.boxplot(x = 'timepoint', y = 'fold_change', hue = 'condition', 
#                 data = final_df[final_df.protein == prot], 
#                 showfliers = False, whis = 1.2)
#     plt.title(prot)
#     plt.show()
#     plt.savefig(f'Boxplot_condition{prot}_woMEC.png', format = 'png')
#     plt.close()
# =============================================================================

# =============================================================================
# sns.catplot(x = 'timepoint', y = 'fold_change', col = 'protein', kind = 'box', hue = 'condition',
#             data = final_df, showfliers = False, col_wrap = 3, height = 3, aspect = 1)            
# =============================================================================

# =============================================================================
# #RESULTS 1
# #process data
# file_pre = 'ResultsSQ_plate'
# dataframes = []
# dics_innernorm = []
# for i in range(1,5):    
#     filename = file_pre + str(i) + '.csv'
#     plate = read_clean_plate(filename)    
#     plate_df, lst_norm = normalize_innerplate(plate)    
#     dataframes.append(plate_df)
#     dics_innernorm.append(lst_norm)
# 
# peptide_final_df = normalize_plates(dataframes, dics_innernorm)    
# protein_df = get_protein_values(peptide_final_df)
# 
# #Omit MEC
# #protein_df = protein_df[protein_df.timepoint != 1]
# 
# treated_df = remove_outliers(protein_df, flag = False, perc = 0.10)
# rescaled_df = rescale(treated_df)
# 
# '''plot function calling'''
# sns.set(context='notebook', style='whitegrid', palette = 'deep', font= 'Helvetica')
# =============================================================================

# =============================================================================
# #hue timepoint or type
# boxplot_all(rescaled_df, hue = 'type')
# boxplot_subplots(rescaled_df, hue = 'type')
# boxplot_all(rescaled_df, hue = 'timepoint')
# boxplot_subplots(rescaled_df, hue = 'timepoint')
# 
# #pretty_plots(rescaled_df)
# 
# boxplot_type(rescaled_df, hue = 'condition')
# 
# lineplot(rescaled_df)
# 
# heatmap(rescaled_df)
# 
# clustermap(rescaled_df)
# 
# #streamflow(rescaled_df, 'Acute')
# #streamflow(rescaled_df, 'Impaired')
# 
# plot_all(rescaled_df)
# 
# =============================================================================

# =============================================================================
# #with MEC removed, rename timepoint and plot figs
# final_df = rescaled_df.copy()
# lst = []
# for row in range(len(final_df)):
#     lst.append(final_df.iloc[row]['timepoint'] -1 )
# final_df.timepoint = lst
# 
# final_df = final_df[final_df.timepoint < 7]
# =============================================================================

#lineplot(rescaled_df)

#clustermap(final_df)

# =============================================================================
# #mfuzz tables
# proteins = rescaled_df.protein.unique()
# for i in range(len(proteins)):
#     prot_out = rescaled_df[rescaled_df.protein == proteins[i]]
#     prot_out = prot_out[prot_out.timepoint < 8]
#     prot_out = prot_out.pivot('patient', 'timepoint', 'total_area_protein')
#     prot_out.to_csv(f'{proteins[i]}_pivot.csv')
# =============================================================================
