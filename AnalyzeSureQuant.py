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
from scipy import stats
from statsmodels.stats import multitest
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.integrate as integrate
from decimal import Decimal
import os 


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


'''
PROCESSING AND CLEANING UP FUNCTIONS
'''


def read_clean_plate(plate_file):
    '''Read plate, convert columns to an easier to select format,
    remove NA values, create column to store the normalized values'''
    
    plate = pd.read_csv(plate_file, sep=',')
    plate.columns = plate.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
    plate = plate.sort_values(by='replicate_name')
    plate = plate.fillna(0)
    plate['total_area_fragment_norm'] = 0
    
    return plate


def ReadChrom():
    '''Read all chromatograms, integrate area and calculate factor to scale values'''
    
    samples = range(1,366)
    df = pd.DataFrame(columns=['RetentionTime','Intensity'])
    folder = r'Chrom'
    
    for s in samples:
        filename = os.path.join(folder,f'20200218_SS_MK_EVO_WF_{s}.raw_Ms_TIC_chromatogram.txt')
        temp = pd.read_csv(filename, sep='\t')
    
        df = df.join(temp, lsuffix='', rsuffix=str(s), how='outer')
    
    df = df.drop(['RetentionTime','Intensity'], axis=1)
    
    lst_area = []
    for s in samples:
        lst_area.append(integrate.simps(df[f'Intensity{str(s)}'].dropna(), 
                                        df[f'RetentionTime{str(s)}'].dropna())) 
    factor = []
    median = stat.median(lst_area)
    for area in lst_area:
        factor.append(median / area)
    
    return factor


def ChromNorm(dataframe, factor):
    '''Normalize based on chromatogram factor values'''
    
    df = dataframe[dataframe.fragment_ion == 'precursor']
    light_df = pd.DataFrame(columns = df.columns)
    
    seqs = df.peptide_sequence.unique()
    
    for seq in seqs:    
        pepdf = df[df['peptide_sequence'] == seq]
        prec_masses = pepdf.precursor_mz.unique()
        light_mass = min(prec_masses)
        pepdf_light = pepdf[pepdf.precursor_mz == light_mass]
        light_df = light_df.append(pepdf_light)
    
    norm_df = pd.DataFrame(columns = light_df.columns)
    for i, val in enumerate(factor):
        rep_df = light_df[light_df.replicate_name == i+1]
        col = np.array(rep_df.total_area_fragment)
        rep_df['total_area_fragment_norm'] =  col * val
        norm_df = norm_df.append(rep_df)
    
    return norm_df


def get_protein_values2(master_df):      
    '''Average values of all peptides for a given protein, across all samples.
    Keep only the columns we need for analysis and clean up a bit.'''
    
    #average proteins based on peptides, keep columns of interest
    protein_df = master_df.groupby(['protein_name','replicate_name', 'sample', 'condition',
           'timepoint', 'type'], as_index = False)['total_area_fragment_norm'].mean()
    
    #clean columns fro later visualization steps
    protein_df = protein_df.rename(columns = {'total_area_fragment_norm':'total_area_protein'})
    protein_df['protein'] = protein_df.protein_name.str[10:-6]
    protein_df['patient'] = protein_df['sample'].str[:-3]
    
    return protein_df


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
        #light_lst = light_lst * (heavy_lst / stat.median(heavy_lst))
        light_lst = light_lst * (stat.median(heavy_lst) / heavy_lst)
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
                median = dics_innernorm[d][seq]        
            median_lst[i].append(median)
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
        #median_lst[s] = np.array(median_lst[s]) / med
        median_lst[s] = med / np.array(median_lst[s])
    
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
    Keep only the columns we need for analysis and clean up a bit.'''
    
    #average proteins based on peptides, keep columns of interest
    protein_df = master_df.groupby(['protein_name','replicate_name', 'sample', 'condition',
           'timepoint', 'type'], as_index = False)['total_area_plate_norm'].mean()
    
    #clean columns fro later visualization steps
    protein_df = protein_df.rename(columns = {'total_area_plate_norm':'total_area_protein'})
    protein_df['protein'] = protein_df.protein_name.str[10:-6]
    protein_df['patient'] = protein_df['sample'].str[:-3]
    
    return protein_df


def rescale(master_df):
    '''Rescale values from 0 - 100, for each protein individually'''
    
    scaled_df = pd.DataFrame(columns = master_df.columns)
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
        
        #remove 0s are below the limit of detection
        #master_df = master_df[master_df['total_area_protein'] != 0]

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


def fold_change(rescaled_df):
    '''Creates a column and calculates fold change from timepoint 1'''
    
    rescaled_df['fold_change'] = 0

    for prot in rescaled_df.protein.unique():
        
        for pat in rescaled_df.patient.unique():
            
            ind = rescaled_df[(rescaled_df.protein == prot) & (rescaled_df.patient == pat)].index
            init_val = rescaled_df.loc[ind][rescaled_df.loc[ind].timepoint == 1].total_area_protein.iloc[0]
            #print(rescaled_df.loc[ind][rescaled_df.loc[ind].timepoint == 1][['protein','patient','timepoint']])
            
            for i in ind:
                cell_val = rescaled_df.loc[i,'total_area_protein']
                if init_val == 0:
                    rescaled_df.loc[i,'fold_change'] =  0
                    init_val = cell_val
                else:
                    rescaled_df.loc[i,'fold_change'] =  (cell_val / init_val)
                        
    return rescaled_df

          
'''
VISUALIZATION FUNCTIONS
'''


def pretty_plots(pdata):
    '''Creates 4 plots. Two boxplots of impaired vs acute and MEC vs HYC, 
    as well as two line plots, one for acute and one for impaired samples over time'''
    sns.set(style="ticks")

    # initialize the figure with a logarithmic x axis
    f, ax = plt.subplots(figsize=(7, 6))
    ax.set_xscale("log")    
    sns.boxplot(x="total_area_protein", y="protein", hue="condition", data=pdata, palette="vlag")
    
    # tweak the visual presentation
    ax.xaxis.grid(True)
    ax.set(ylabel="")
    ax.set(xlim = (0.000001,110))
    ax.legend(loc='center left', bbox_to_anchor=(1.08, 0.95), ncol=1)
    sns.despine(trim=True, left=True)
    
    # start new plot
    plt.figure()
    
    # lineplot for acute
    ac_line = sns.lineplot(x="timepoint", y="total_area_protein",
                 hue="protein", data=pdata[pdata.condition == "Acute"])
    ac_line.legend(loc='center left', bbox_to_anchor=(1.08, 0.8), ncol=1)
    ac_line.set_title("Acute")

    # start new plot
    plt.figure()
    
    # lineplot for impaired
    imp_line = sns.lineplot(x="timepoint", y="total_area_protein",
                 hue="protein", data=pdata[pdata.condition == "Impaired"])
    imp_line.legend(loc='center left', bbox_to_anchor=(1.08, 0.8), ncol=1)
    imp_line.set_title("Impaired")
    
    # initialize the figure with a logarithmic x axis
    f, ax = plt.subplots(figsize=(7, 6))
    ax.set_xscale("log")    
    sns.boxplot(x="total_area_protein", y="protein", hue="type", data=pdata, palette="vlag")
    
    # tweak the visual presentation
    ax.xaxis.grid(True)
    ax.set(ylabel="")
    ax.set(xlim = (0.000001,110))
    ax.legend(loc='center left', bbox_to_anchor=(1.08, 0.95), ncol=1)
    sns.despine(trim=True, left=True)
    
    
def boxplot_all(rescaled_df, hue):
    '''Boxplots, all acute/impaired, based on hue'''
    
    fig = plt.figure(figsize = (9,6))
    ax = sns.boxplot(x='protein', y='total_area_protein', hue=hue, 
                     data= rescaled_df[rescaled_df.condition == 'Acute'], palette='vlag', showfliers = False)
    # =============================================================================
    # for tick in ax.xaxis.get_major_ticks():
    #     tick.label.set_fontsize(10)
    # ax.tick_params( rotation=30)
    # =============================================================================
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_title('Acute healing, all timepoints')
    ax.set_ylabel('Protein relative concentration')
    ax.set_xlabel('Protein')
    plt.tight_layout()
    plt.savefig('Acute_'+hue+'.png', dpi = 300, format = 'png')
    plt.close()
    
    fig = plt.figure(figsize = (9,6))
    ax = sns.boxplot(x='protein', y='total_area_protein', hue=hue, 
                     data= rescaled_df[rescaled_df.condition == 'Impaired'], palette='vlag', showfliers = False)
    # =============================================================================
    # for tick in ax.xaxis.get_major_ticks():
    #     tick.label.set_fontsize(10)
    # ax.tick_params( rotation=30)
    # =============================================================================
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_title('Impaired healing, all timepoints')
    ax.set_ylabel('Protein relative concentration')
    ax.set_xlabel('Protein')
    plt.tight_layout()
    plt.savefig('Impaired_'+hue+'.png', dpi = 300, format = 'png')
    plt.close()


def boxplot_subplots(rescaled_df, hue):
    '''Boxplot subplots for all proteins, all acute/impaired, based on hue'''
    
    proteins = rescaled_df.protein.unique()
    fig, ax = plt.subplots(3,2,figsize = (8,12))
    
    x = 0
    for i in range(0, len(proteins), 2):
        y = 0
        
        sns.boxplot(x='protein', y='total_area_protein', hue=hue, 
          data= rescaled_df[(rescaled_df.condition == 'Acute') & (rescaled_df.protein == proteins[i])], 
          palette='vlag', showfliers = False, ax = ax[x,y]) 
   
        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')
    
        if x == 1 and y == 0:
            ax[x,y].set_ylabel('Normalized rel. concentration')
        elif x == 2 and y == 0:
            ax[x,y].set_xlabel('Protein')
            ax[x,y].xaxis.set_label_coords(1.05, -0.1)
            
        y += 1
        sns.boxplot(x='protein', y='total_area_protein', hue=hue, 
          data= rescaled_df[(rescaled_df.condition == 'Acute') & (rescaled_df.protein == proteins[i+1])], 
          palette='vlag', showfliers = False, ax = ax[x,y]) 

        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')    
            
        x += 1
    
    fig.suptitle('Acute healing, all timepoints')
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.savefig('Subplot_'+hue+'_acute_all.png', dpi = 300, format = 'png')
    plt.close()
    
    #second fig
    proteins = rescaled_df.protein.unique()
    fig, ax = plt.subplots(3,2,figsize = (8,12))
    
    x = 0
    for i in range(0, len(proteins), 2):
        y = 0
        
        sns.boxplot(x='protein', y='total_area_protein', hue=hue, 
          data= rescaled_df[(rescaled_df.condition == 'Impaired') & (rescaled_df.protein == proteins[i])], 
          palette='vlag', showfliers = False, ax = ax[x,y]) 
  
        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')
        
        if x == 1 and y == 0:
            ax[x,y].set_ylabel('Normalized rel. concentration')
        elif x == 2 and y == 0:
            ax[x,y].set_xlabel('Protein')
            ax[x,y].xaxis.set_label_coords(1.05, -0.1)
        
        y += 1
        sns.boxplot(x='protein', y='total_area_protein', hue=hue, 
          data= rescaled_df[(rescaled_df.condition == 'Impaired') & (rescaled_df.protein == proteins[i+1])], 
          palette='vlag', showfliers = False, ax = ax[x,y]) 

        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')    
        
        x += 1
    
    fig.suptitle('Impaired healing, all timepoints')
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.savefig('Subplot_'+hue+'_impaired_all.png', dpi = 300, format = 'png')
    plt.close()


def boxplot_type(rescaled_df, hue):
    
    proteins = rescaled_df.protein.unique()
    fig, ax = plt.subplots(3,2,figsize = (8,12))
    
    x = 0
    for i in range(0, len(proteins), 2):
        y = 0
        
        sns.boxplot(x='type', y='total_area_protein', hue=hue, 
            data= rescaled_df[rescaled_df.protein == proteins[i]], 
            palette='vlag', showfliers = False, ax = ax[x,y]) 
 
        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')
        ax[x,y].set_title(proteins[i])
        
        if x == 1 and y == 0:
            ax[x,y].set_ylabel('Normalized rel. concentration')
        elif x == 2 and y == 0:
            ax[x,y].set_xlabel('Protein')
            ax[x,y].xaxis.set_label_coords(1.05, -0.1)
        
        y += 1
        sns.boxplot(x='type', y='total_area_protein', hue=hue, 
          data= rescaled_df[rescaled_df.protein == proteins[i+1]], 
          palette='vlag', showfliers = False, ax = ax[x,y]) 

        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')    
        ax[x,y].set_title(proteins[i + 1])
        
        if x == 2 and y == 1:
            ax[x,y].legend()
        x += 1
    
    fig.suptitle(f'All timepoints, {hue}')
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig('Subplot_'+hue+'type_all.png', dpi = 300, format = 'png')
    plt.close()

def lineplot(rescaled_df):
    sns.set(style="ticks")
    colors = sns.color_palette('muted')
    proteins = rescaled_df.protein.unique()
    proteins = ['TNFA', 'IL1B', 'MMP2', 'S10A8', 'S10A9', 'MMP9']
    fig, ax = plt.subplots(2,3,figsize = (12,8))
    
    x = 0
    for i in range(0, len(proteins), 3):
        
        y = 0
        sns.lineplot(x="timepoint", y="total_area_protein",
             hue="condition", data=rescaled_df[rescaled_df.protein == proteins[i]], color = colors[i], ax = ax[x,y])
        
        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')
        ax[x,y].set_title(proteins[i])
        
        if y == 0:
            ax[x,y].set_ylabel('Normalized rel. concentration', fontsize = 14)
        if x == 1 and y < 3:
            ax[x,y].set_xlabel('Sampling timepoint')
        #elif x == 1 and y == 0:
            #ax[x,y].set_xlabel('Protein')
            #ax[x,y].xaxis.set_label_coords(1.05, -0.1)
        
        y += 1
        sns.lineplot(x="timepoint", y="total_area_protein",
                hue="condition", data=rescaled_df[rescaled_df.protein == proteins[i+1]], color = colors[i], ax = ax[x,y])
    
        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')    
        ax[x,y].set_title(proteins[i + 1])
        
        if x == 1 and y < 3:
            ax[x,y].set_xlabel('Sampling timepoint')
        
        y += 1
        sns.lineplot(x="timepoint", y="total_area_protein",
                hue="condition", data=rescaled_df[rescaled_df.protein == proteins[i+2]], color = colors[i], ax = ax[x,y])
        
        if x == 0: 
            ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')    
        ax[x,y].set_title(proteins[i + 2])
        
        if x == 1 and y < 3:
            ax[x,y].set_xlabel('Sampling timepoint')
        x += 1
        
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    fig.suptitle(f'Proteins')
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig('Subplot_lineplot_all1.png', dpi = 300, format = 'png')
    plt.close()


def heatmap(master_df):  
    '''Heatmap for each individual protein, patient-timepoint'''
    
    cmaps = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']

    proteins = master_df.protein.unique()
    
    for i in range(len(proteins)):
        
        cbarval = False
        if (i +1) %2 == 0:
            cbarval = True
        
        prot_frame = master_df[master_df.protein == proteins[i]]
        #https://seaborn.pydata.org/generated/seaborn.heatmap.html
        heat_frame = prot_frame.pivot('patient','timepoint', 'total_area_protein')
        
        plt.figure(figsize=(6,8))
        ax = plt.axes()
        sns.heatmap(heat_frame, vmin = 0, vmax = 25, cmap = cmaps[i], ax = ax, cbar= cbarval)
        ax.set_title(proteins[i])
        plt.savefig(f'Heatmap_{proteins[i]}.png', dpi = 300, format = 'png')
        plt.close()


def clustermap(master_df):
    '''Clustermap for each individual protein, patient-timepoint'''
    
    cmaps = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']

    proteins = master_df.protein.unique()
    
    for i in range(len(proteins)):
        
        prot_frame = master_df[master_df.protein == proteins[i]]
        heat_frame = prot_frame.pivot('patient','timepoint', 'total_area_protein')
        cluster_frame = heat_frame[list(heat_frame.columns)[:6]]
        #sample 325 is missing
        
        for col in cluster_frame.columns:
            cluster_frame[col] = np.log(cluster_frame[col])
        
        cluster_frame = cluster_frame.replace([np.inf, -np.inf], np.nan)
        cluster_frame = cluster_frame.fillna(0)
        
        sns.clustermap(cluster_frame, cmap = cmaps[i], linewidths=.5, yticklabels=True, figsize =(8, 12))
        
        plt.suptitle(f'{proteins[i]}')
        plt.savefig(f'Clustermap_{proteins[i]}.png', dpi = 300, format = 'png')
        plt.close()


def streamflow(master_df, cond):
    
    cond = master_df[master_df.condition == cond]
    cond = cond.groupby(['protein', 'timepoint'], as_index = False).mean()
    cond = cond.pivot('timepoint', 'protein', 'total_area_protein')
    
    plt.figure(figsize=(9,6))
    plt.stackplot(cond.index,cond.T, labels = cond.columns, baseline = 'weighted_wiggle')
    #plt.ylim([-18,18])
    plt.title(f'All proteins, acute healing')
    plt.xlabel('Timepoint')
    plt.ylabel('Weighted protein amount')
    plt.subplots_adjust(right=0.82)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    plt.savefig('Streamflow_acute.png', dpi = 300, format = 'png')


def plot_all(master_df):
    
    for patient in master_df.patient.unique():
        
        test = master_df[(master_df['patient'] == patient)]
        sns.set(style = 'ticks')
        g = sns.relplot(x='timepoint', y='total_area_protein', col ='protein', data = test, col_wrap= 3, 
                    height=3, aspect=1, kind = 'line')
        
        plt.subplots_adjust(top=0.92)
        g.fig.suptitle(f'Patient {patient}')
        plt.savefig(f'{patient}.png', dpi = 300, format = 'png')
        plt.close()
        
        
def lineplot2(rescaled_df):
    sns.set(style="ticks")
    colors = sns.color_palette('muted')
    proteins = rescaled_df.protein.unique()
    fig, ax = plt.subplots(3,3,figsize = (12,12))
    
    x = 0
    for i in range(0, len(proteins), 3):
        
        y = 0
        sns.lineplot(x="timepoint", y="total_area_protein",
             hue="condition", data=rescaled_df[rescaled_df.protein == proteins[i]], color = colors[i], ax = ax[x,y])
        
        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')
        ax[x,y].set_title(proteins[i])
        
        if y == 0:
            ax[x,y].set_ylabel('Rel. concentration', fontsize = 14)
        if x == 2 and y < 3:
            ax[x,y].set_xlabel('Sampling timepoint')
        
        y += 1
        sns.lineplot(x="timepoint", y="total_area_protein",
                hue="condition", data=rescaled_df[rescaled_df.protein == proteins[i+1]], color = colors[i], ax = ax[x,y])
    
        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')    
        ax[x,y].set_title(proteins[i + 1])
        
        if x == 2 and y < 3:
            ax[x,y].set_xlabel('Sampling timepoint')
        
        y += 1
        sns.lineplot(x="timepoint", y="total_area_protein",
                hue="condition", data=rescaled_df[rescaled_df.protein == proteins[i+2]], color = colors[i], ax = ax[x,y])
        
        ax[x,y].get_legend().set_visible(False)
        ax[x,y].set_ylabel('')
        ax[x,y].set_xlabel('')    
        ax[x,y].set_title(proteins[i + 2])
        
        if x == 2 and y < 3:
            ax[x,y].set_xlabel('Sampling timepoint')
        x += 1
        
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    fig.suptitle(f'Proteins')
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig('Subplot_lineplot_all.png', dpi = 300, format = 'png')
    plt.savefig('Subplot_lineplot_all.svg', format = 'svg')
    plt.close()


def final_box(final_df, col):
    sns.set(context='notebook', style='darkgrid', palette = 'deep', font= 'Helvetica')
    
    c = 0
    colors = ['blue', 'orange', 'green', 'red', 'purple', 'olive', 'cyan', 'brown', 'pink']
    titles = final_df.protein_name.unique()
    
    for prot in final_df.protein.unique():
        
        fig, ax = plt.subplots(2,1, figsize=(6,8))
              
        sns.boxplot(x='timepoint', y=col,
                    data = final_df[(final_df.protein == prot) & (final_df.condition == 'Impaired')], 
                    ax = ax[0], showfliers = False, whis = 1.5, color = colors[c])
        
        ax[0].set_ylabel('Relative concentration')
        ax[0].set_title('Impaired', size=14)
        ymin, ymax = ax[0].get_ylim()
        ax[0].set_ylim([ymin+(ymin*0.1), ymax+(ymax*0.1)])
        
        sns.boxplot(x='timepoint', y=col, 
                    data = final_df[(final_df.protein == prot) & (final_df.condition == 'Acute')], 
                    ax = ax[1], showfliers = False, whis = 1.5, color = 'tab:'+colors[c])
        
        ymin, ymax = ax[1].get_ylim()
        ax[1].set_ylim([ymin+(ymin*0.1), ymax+(ymax*0.1)])
        
        ax[1].set_ylabel('Relative concentration')
        ax[1].set_title('Acute', size= 14)
        
        acc = titles[c].split('|')[1]
        plt.suptitle(f'{prot} - {acc}', x = 0.24, size = 15)
        #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout(pad = 1.5)
        plt.subplots_adjust(top=0.92, left=0.1, right=0.95)
        plt.savefig(prot+f'_{col}box.svg', format='svg')
        
        c += 1
    
#run things under here
  
# =============================================================================
# file_pre = 'ResultsSQ_plate'
# dataframes = []
# dics_internorm = []
# 
# for i in range(1,5):
#     
#     filename = file_pre + str(i) + '.csv'
#     plate = read_clean_plate(filename)
#     
#     plate_df, lst_norm = normalize_innerplate(plate)
#     
#     dataframes.append(plate_df)
#     dics_internorm.append(lst_norm)
# 
# peptide_final_df = normalize_plates(dataframes, dics_internorm)    
# 
# protein_df = get_protein_values(peptide_final_df)
# #protein_df.to_excel('Norm_prot.xlsx')
# 
# treated_df = remove_outliers(protein_df, flag = False, perc = 0.10)
# 
# rescaled_df = rescale(treated_df)
# #rescaled_df.to_excel('Scaled_prot.xlsx')
# 
# sns.set(context='notebook', style='whitegrid', palette = 'deep', font= 'Helvetica')
# 
# '''plot function calling'''
# #hue timepoint or type
# #boxplot_all(rescaled_df, hue = 'type')
# #boxplot_subplots(rescaled_df, hue = 'type')
# #pretty_plots(rescaled_df)
# #boxplot_type(rescaled_df, hue = 'condition')
# #lineplot(rescaled_df)
# #heatmap(rescaled_df)
# =============================================================================


'''PLOTS'''
# =============================================================================
# temp = rescaled_df[(rescaled_df.protein == 'TNFA') & (rescaled_df.condition == 'Impaired')]
# sns.violinplot(x = 'timepoint', y = 'total_area_protein', scale="count", inner="quartile", 
#                data= temp, bw=.5)
# =============================================================================
        
# =============================================================================
# sns.catplot(x="timepoint", y="total_area_protein", hue = 'condition', col="protein",
#                 col_wrap = 2, data=rescaled_df, 
#                 kind="violin", split=False, scale="count", inner="quartile", 
#                 scale_hue=False, bw=.5, height=5, aspect=1);
# plt.savefig('violin_all.png', dpi = 300, format = 'png')
# =============================================================================
        
# =============================================================================
# sns.catplot(x="timepoint", y="total_area_protein", col="protein",
#                 col_wrap = 2, data=rescaled_df[rescaled_df.condition == 'Impaired'], 
#                 kind="violin", scale="count", inner="quartile", 
#                 scale_hue=False, bw=.2, height=5, aspect=1,);
# plt.suptitle('Impaired, all timepoints')
# plt.tight_layout()
# plt.savefig('violin_impaired_very_low_bw.png', dpi = 300, format = 'png')
# =============================================================================
    
       
# =============================================================================
# test = rescaled_df[(rescaled_df['patient'] == '02A-001')]
#     
# g = sns.catplot(x='timepoint', y='total_area_protein', col ='protein', data = test, col_wrap= 2, 
#             height=4, aspect=0.9, kind = 'bar', color = sns.color_palette('deep')[0])
# 
# plt.subplots_adjust(top=0.92)
# g.fig.suptitle('Patient 02A-001')
# =============================================================================

    
'''STATISTICS'''
# =============================================================================
# #CHECK IF SAMPLES ARE NORMALLY DISTRIBUTED
# mod_df = rescaled_df[rescaled_df.condition == 'Impaired']
# for i in mod_df.protein.unique():
#     for k in mod_df.timepoint.unique():
#         test = mod_df[(mod_df.protein == i) & (mod_df.timepoint == k)]
# 
#         k2, p = stats.normaltest(np.array(test.total_area_protein), nan_policy = 'omit')
#         alpha = 1e-3
#         if p > alpha:
#             print(i, k, p)
# 
# test = rescaled_df[(rescaled_df.protein == 'IL1B') & (rescaled_df.timepoint == 7) &
#                    (rescaled_df.condition == 'Impaired')]
# 
# plt.hist(test.total_area_protein)
# 
# for i in rescaled_df.protein.unique():
#     test = rescaled_df[(rescaled_df.protein == i) & (rescaled_df.condition == 'Impaired')]
#     h, p = stats.kruskal(*[group["total_area_protein"].values for name, group in test.groupby("timepoint")])
#     print(i,p)
#     
# #MEC vs HYC
# for k in rescaled_df.protein.unique():
#     prot_df = rescaled_df[rescaled_df.protein == k]
#     testM = prot_df[prot_df.type == 'MEC']
#     testH = prot_df[prot_df.type == 'HYC']
#     u, p = stats.mannwhitneyu(testM.total_area_protein, testH.total_area_protein)
#     print(k, p)
# 
# #All timepoints, acute vs impaired
# ind_l = []
# p_l = []
# t_l = []
# for k in rescaled_df.protein.unique():
#     prot_df = rescaled_df[rescaled_df.protein == k]
#     for i in range(1,8):
#         testA = prot_df[(prot_df.condition == 'Acute') & (rescaled_df.timepoint == i)]
#         testI = prot_df[(prot_df.condition == 'Impaired') & (rescaled_df.timepoint == i)]
#         u, p = stats.mannwhitneyu(testA.total_area_protein, testI.total_area_protein)
#         print(p, i, k)
#         ind_l.append(k)
#         t_l.append(i)
#         p_l.append(p)
# 
# 
# reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.01, method = 'bonferroni')
# 
# fh = open('MannwhitneyU_res_conditionpertimepoint.txt', 'w')
# 
# for i in range(len(reject)):
#     if reject[i]:
#         fh.write(f'{ind_l[i]}\t{t_l[i]} - {p_lnew[i]}\n')
# fh.close()
# 
# 
# ind_l = []
# p_l = []
# t_l = [] 
# t2_l = []
# 
# for k in rescaled_df.protein.unique():
#     s = 1
#     c = 2
#     
#     while s <= 7:
#         for i in range(c,8):
#             sample1 = rescaled_df[(rescaled_df.protein == k) & (rescaled_df.condition == 'Impaired') 
#                                   & (rescaled_df.timepoint == s)]
#             sample2 = rescaled_df[(rescaled_df.protein == k) & (rescaled_df.condition == 'Impaired') 
#                                   & (rescaled_df.timepoint == i)]
#             u, p = stats.mannwhitneyu(sample1.total_area_protein, sample2.total_area_protein)
#             ind_l.append(k)
#             t_l.append(s)
#             t2_l.append(i)
#             p_l.append(p)
#         
#         c += 1
#         s += 1
#     
# reject, p_lnew, asid, abon = multitest.multipletests(p_l,alpha = 0.01, method = 'bonferroni')
# 
# fh = open('MannwhitneyU_res_proteinpertimepoint_impaired.txt', 'w')
# 
# for i in range(len(reject)):
#     if reject[i]:
#         fh.write(f'{ind_l[i]}\t{t_l[i]}+{t2_l[i]} - {p_lnew[i]}\n')
# fh.close()
#     
# fh.close()
# =============================================================================

      
'''PROCESSING'''
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