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
    
    plate = pd.read_csv(plate_file, sep=';')
    plate.columns = plate.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
    plate = plate.sort_values(by='replicate_name')
    plate = plate.fillna(0)
    plate['total_area_fragment_norm'] = 0
    
    return plate

def normalize_innerplate(plate):
    
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

final_df = normalize_plates(dataframes, dics_internorm)    
        
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