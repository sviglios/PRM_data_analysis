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
final_df = fold_change(rescaled_df)    
    
# =============================================================================
# #Rename timepoints, remove no8 from acute
# final_df = rescaled_df.copy()
# lst = []
# for row in range(len(final_df)):
#     lst.append(final_df.iloc[row]['timepoint'] -1 )
# final_df.timepoint = lst
# 
# final_df = final_df[final_df.timepoint < 7]
# =============================================================================

# =============================================================================
# #plot
# lineplot2(rescaled_df)
# plot_all(rescaled_df)
# clustermap(rescaled_df)
# =============================================================================

#boxplots all
# =============================================================================
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
