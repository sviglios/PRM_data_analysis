# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 12:50:59 2020

@author: konka
"""

#import functions
from AnalyzeSureQuant import *

#process data
file_pre = 'ResultsSQ_plate'
dataframes = []
dics_innernorm = []
for i in range(1,5):    
    filename = file_pre + str(i) + '.csv'
    plate = read_clean_plate(filename)    
    plate_df, lst_norm = normalize_innerplate(plate)    
    dataframes.append(plate_df)
    dics_innernorm.append(lst_norm)

peptide_final_df = normalize_plates(dataframes, dics_innernorm)    
protein_df = get_protein_values(peptide_final_df)
treated_df = remove_outliers(protein_df, flag = False, perc = 0.10)
rescaled_df = rescale(treated_df)

'''plot function calling'''
# =============================================================================
# sns.set(context='notebook', style='whitegrid', palette = 'deep', font= 'Helvetica')
# 
# #hue timepoint or type
# boxplot_all(rescaled_df, hue = 'type')
# boxplot_subplots(rescaled_df, hue = 'type')
# 
# pretty_plots(rescaled_df)
# 
# boxplot_type(rescaled_df, hue = 'condition')
# 
# lineplot(rescaled_df)
# 
# heatmap(rescaled_df)
# 
# clustermap(rescaled_df)
# =============================================================================



        
        
        
        
        
        
        