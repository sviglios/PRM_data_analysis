# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:06:13 2020

@author: konka
"""

'''
Use 3 apostrophes for longer comments, # for smaller comments.

Every task should be executed in functions for better flow and understanding
of the code. Comment function in short after every function and job you write.
No need for classes.

Would be nice to include links if code comes from somewhere, so both can read the source.

In Spyder, Ctrl+4 after highlighting the code comments it out. Do that if you want to 
replace or keep code but you don't want to delete it. We can clean up later. Placing the cursor
inside the commented code, Ctrl+5 uncomments it. F5 is running the script, F9 is running selection.
Switch to dark mode, monokai theme imo. Tools -> Preferences for customization.
'''

import pandas as pd
import statistics as stat
import numpy as np
from sklearn.preprocessing import MinMaxScaler

import matplotlib.pyplot as plt
import seaborn as sns

'''
sum ion intensity -> peptide signal
normalize for heavy pep signal sample wise
normalize for heave pep signal plate wise
include peptides function for each protein. That is on the avg normalization on the peptide lvl
exclude samples function on the protein signal lvl

SOME PEPTIDES ARE ONLY GOOD IN SOME PLATES. FIGURE OUT HOW TO EXCLUDE THEM FOR SOME SAMPLES
SOME PEPTIDES DO NOT HAVE VALUES IN THE LIGHT VERSION IN SOME SAMPLES
FEWER PEPTIDES DO NOT HAVE VALUES FOR EITHER LIGHT OR HEAVY IN A FEW SAMPLES
'''


def read_clean_plate(plate_file):
    '''Read plate, convert columns to an easier to select format,
    remove NA values, create column to store the normalized values'''
    
    plate = pd.read_csv(plate_file, sep=';')
    plate.columns = plate.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
    plate = plate.sort_values(by='replicate_name')
    plate = plate.fillna(0)
    plate['total_area_fragment_norm'] = 0
    
    return plate


def normalize_innerplate(plate):
    '''Normalization of light peptide fragment area values
    based on the heavy peptide area. The median of the heavy peptide
    values is used to rescale all light peptide values'''
    
    #make plate intermediate dataframe after inner plate normalization
    plate_df = pd.DataFrame(columns = plate.columns) 

    seqs = plate.peptide_sequence.unique()
    
    #transform to a loop that goes over all peptides later
    #pepdf = plate[plate['peptide_sequence'] == seqs[3]]
    
    dic_inter_norm = {}
    
    for seq in seqs:
        
        pepdf = plate[plate['peptide_sequence'] == seq]
        
        #get masses of heavy and light
        prec_masses = pepdf.precursor_mz.unique()
        light_mass, heavy_mass = min(prec_masses), max(prec_masses)
        
        pepdf = pepdf[pepdf.fragment_ion == 'precursor']
        
        #get heavy and light peptide values for all samples, single peptide
        pepdf_heavy, pepdf_light = pepdf[pepdf.precursor_mz == heavy_mass], pepdf[pepdf.precursor_mz == light_mass]
        heavy_lst, light_lst = np.array(pepdf_heavy.total_area_fragment), np.array(pepdf_light.total_area_fragment)
          
        #normalize light values to heavy values
        light_lst = light_lst * (heavy_lst / stat.median(heavy_lst))
        inter_plate = stat.median(heavy_lst)
        
        pepdf_light['total_area_fragment_norm'] = light_lst
            
        #append to plate df, keep values of heavy median in a tuple for inter plate normalization
        plate_df = plate_df.append(pepdf_light)
        dic_inter_norm[seq] = inter_plate

    return plate_df, dic_inter_norm

    
def normalize_plates(dataframes, dics_innernorm):
    '''Normalize on plate level. Compare medians of heavy peptides,
    for each peptide, for each plate. Normalize between them, refactor
    all light peptide values with that. Use only peptides that are 
    present in all plates'''
    
    #this is really crappy way to do this, and function in general.
    
    #find sequences present in dataframes
    seqset = []
    for i in range(len(dataframes)):    
        seqs_df = dataframes[i].peptide_sequence.unique() 
        for s in seqs_df:
            seqset.append(s)    
        
    seqset = set(seqset)
    seqset = list(seqset)
    
    #get values for heavy peptide median for all plates, default 0 if seq not in plate
    seq_lst = []
    median_lst = []
    for i in range(len(seqset)):    
        seq = seqset[i]   
        seq_lst.append(seq)
        median_lst.append([])    
        for d in range(len(dataframes)):        
            df = dataframes[d]
            median = 0        
            if seq in list(df.peptide_sequence.unique()):
                median = dics_internorm[d][seq]        
            median_lst[i].append(median)
    
    #remove seqs not present in all plates
    ind_del = []
    for s in range(len(seq_lst)):
        tempmed = median_lst[s]
        if 0 in tempmed:
            ind_del.append(s)
    for ind in sorted(ind_del, reverse=True):
        del seq_lst[ind]
        del median_lst[ind]
    
    #create final dataframe    
    master_df = pd.DataFrame(columns = dataframes[0].columns)
    
    #create factor list
    for s in range(len(seq_lst)):
        medls = median_lst[s]
        med = stat.median(medls)
        median_lst[s] = np.array(median_lst[s]) / med
    
    for s in range(len(seq_lst)):
        
        seq = seq_lst[s]
        
        for d in range(len(dataframes)):
            
            factor = median_lst[s][d]
            df = dataframes[d] 
            df_seq = df[df.peptide_sequence == seq]
            df_seq['total_area_plate_norm'] = list(factor * df_seq.total_area_fragment_norm)
            df_seq['plate_id'] = d + 1
            
            master_df = master_df.append(df_seq)
    
    return master_df


def get_protein_values(master_df):      
    '''Average values of all peptides for a given protein, across all samples.
    Keep only the columns we need for analysis and clean up a bit. Rescale to 
    a 0-1000 scale as well'''
    
    #average proteins based on peptides, keep columns of interest
    protein_df = master_df.groupby(['protein_name','replicate_name', 'sample', 'condition',
           'timepoint', 'type'], as_index = False)['total_area_plate_norm'].mean()
    
    protein_df = protein_df.rename(columns = {'total_area_plate_norm':'total_area_protein'})
    protein_df['protein'] = protein_df.protein_name.str[10:-6]
    
    return protein_df


def rescale(master_df):
    '''Rescale values from 0 - 100, for each protein individually'''
    
    scaled_df = pd.DataFrame(columns = protein_df.columns)
    proteins = master_df.protein.unique()
    
    for prot in proteins:
        temp = master_df[master_df.protein == prot]
        
        #https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MinMaxScaler.html
        #https://stackoverflow.com/questions/24645153/pandas-dataframe-columns-scaling-with-sklearn
        scaler = MinMaxScaler(feature_range = (0,100))    
        temp[['total_area_protein']] = scaler.fit_transform(temp[['total_area_protein']])
        
        scaled_df = scaled_df.append(temp)
    
    return scaled_df


def remove_outliers(master_df, flag, perc):
    '''Remove outliers, if flag is True'''
    
    if flag:
        
        clean_df = pd.DataFrame(columns = master_df.columns)
        proteins = master_df.protein.unique()
        
        for prot in proteins:
            
            temp = master_df[master_df.protein == prot]
            
            qh = temp.total_area_protein.quantile(1-perc)
            ql = temp.total_area_protein.quantile(perc)
            temp = temp[(temp.total_area_protein > ql) & (temp.total_area_protein < qh)]
            clean_df = clean_df.append(temp)
            
        return clean_df
            
    else:
        return master_df



#run things under here
      
file_pre = 'ResultsSQ_plate'
dataframes = []
dics_internorm = []

for i in range(1,5):
    
    filename = file_pre + str(i) + '.csv'
    plate = read_clean_plate(filename)
    
    plate_df, lst_norm = normalize_innerplate(plate)
    
    dataframes.append(plate_df)
    dics_internorm.append(lst_norm)

peptide_final_df = normalize_plates(dataframes, dics_internorm)    

protein_df = get_protein_values(peptide_final_df)
#protein_df.to_excel('Norm_prot.xlsx')

treated_df = remove_outliers(protein_df, flag = True, perc = 0.10)

rescaled_df = rescale(treated_df)
#rescaled_df.to_excel('Scaled_prot.xlsx')


# =============================================================================
# sns.set(context='notebook', style='whitegrid', palette = 'deep', font= 'Helvetica')
# 
# fig = plt.figure(figsize = (9,6))
# 
# ax = sns.boxplot(x='protein', y='total_area_protein', hue='timepoint', data= rescaled_df[rescaled_df.condition == 'Acute'], palette='vlag')
# # =============================================================================
# # for tick in ax.xaxis.get_major_ticks():
# #     tick.label.set_fontsize(10)
# # ax.tick_params( rotation=30)
# # =============================================================================
# #ax.set_ylim(0,100)
# ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
# #ax.set_title('Acute healing, all timepoints')
# ax.set_ylabel('Protein relative concentration')
# ax.set_xlabel('Protein')
# plt.tight_layout()
# =============================================================================
    
# =============================================================================
# plate1 = pd.read_csv('ResultsSQ_plate1.csv', sep=';')
# 
# #fix column names so it's easier to select
# plate1.columns = plate1.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
# 
# #for later, pick peptides for quantification
# #plate1 = plate1[plate1.peptide_note == 'q']
# 
# plate1 = plate1.sort_values(by='replicate_name')
# plate1 = plate1.fillna(0)
# 
# samples = plate1.replicate_name.unique()
# seqs = plate1.peptide_sequence.unique()
# 
# #transform to a loop that goes over all peptides later
# pepdf = plate1[plate1['peptide_sequence'] == seqs[3]]
# 
# #get masses of heavy and light
# prec_masses = pepdf.precursor_mz.unique()
# light_mass, heavy_mass = min(prec_masses), max(prec_masses)
# 
# # =============================================================================
# # #check areas and if they match
# # sample_pep = pepdf[pepdf.replicate_name == samples[0]]
# # sample_pep_light = sample_pep[sample_pep.precursor_mz == light_mass]
# # 
# # sample_pep_light_prec = sample_pep_light[sample_pep.fragment_ion == 'precursor']
# # sample_pep_light_ion = sample_pep_light[sample_pep_light.fragment_ion != 'precursor']
# # 
# # print(int(sample_pep_light_prec.total_area_fragment), int(sum(sample_pep_light_ion.area)))
# # #they are the same, we continue only with precursor, values in column total_area_fragment
# # =============================================================================
# 
# pepdf = pepdf[pepdf.fragment_ion == 'precursor']
# 
# #get heavy and light peptide values for all samples, single peptide
# pepdf_heavy, pepdf_light = pepdf[pepdf.precursor_mz == heavy_mass], pepdf[pepdf.precursor_mz == light_mass]
# 
# sample_lst = np.array(pepdf_heavy.replicate_name)
# heavy_lst, light_lst = np.array(pepdf_heavy.total_area_fragment), np.array(pepdf_light.total_area_fragment)
# 
# # =============================================================================
# # #check if everything is correct
# # print(sample_lst[5], heavy_lst[5], light_lst[5])
# # print(pepdf[pepdf.replicate_name == 6].total_area_fragment)
# # =============================================================================
# 
# #normalize light values to heavy values
# light_lst = light_lst * (heavy_lst / stat.median(heavy_lst))
# intra_plate = stat.median(heavy_lst)
# 
# pepdf_light['total_area_fragment_norm'] = light_lst
# =============================================================================



# =============================================================================
# plate = pd.read_csv('ResultsSQ_plate1.csv', sep=';')
# 
# plate.columns = plate.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
# plate = plate.sort_values(by='replicate_name')
# plate = plate.fillna(0)
# plate['total_area_fragment_norm'] = 0
# 
# #make plate intermediate dataframe after inner plate normalization
# plate_df = pd.DataFrame(columns = plate.columns) 
# 
# samples = plate.replicate_name.unique()
# seqs = plate.peptide_sequence.unique()
# 
# #transform to a loop that goes over all peptides later
# #pepdf = plate[plate['peptide_sequence'] == seqs[3]]
# 
# lst_intra_norm = []
# 
# for seq in seqs:
#     
#     pepdf = plate[plate['peptide_sequence'] == seq]
#     
#     #get masses of heavy and light
#     prec_masses = pepdf.precursor_mz.unique()
#     light_mass, heavy_mass = min(prec_masses), max(prec_masses)
#     
#     pepdf = pepdf[pepdf.fragment_ion == 'precursor']
#     
#     #get heavy and light peptide values for all samples, single peptide
#     pepdf_heavy, pepdf_light = pepdf[pepdf.precursor_mz == heavy_mass], pepdf[pepdf.precursor_mz == light_mass]
#     heavy_lst, light_lst = np.array(pepdf_heavy.total_area_fragment), np.array(pepdf_light.total_area_fragment)
#       
#     #normalize light values to heavy values
#     light_lst = light_lst * (heavy_lst / stat.median(heavy_lst))
#     intra_plate = stat.median(heavy_lst)
#     
#     pepdf_light['total_area_fragment_norm'] = light_lst
#         
#     #append to plate df, keep values of heavy median in a tuple for intra plate normalization
#     plate_df = plate_df.append(pepdf_light)
#     lst_intra_norm.append((seq,intra_plate))
# =============================================================================