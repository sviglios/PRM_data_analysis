# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 17:05:12 2020

@author: kosta
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.integrate as integrate
from decimal import Decimal
import numpy as np
import os 

samples = [1,10,11,12]
df = pd.DataFrame(columns=['RetentionTime','Intensity'])
df2 = pd.DataFrame(columns=['RetentionTime','Intensity'])
folder = r'C:\Users\kosta\Google Drive\MSdata\Chrom'

for s in samples:
    filename = os.path.join(folder,f'20200218_SS_MK_EVO_WF_{s}.raw_Ms_TIC_chromatogram.txt')
    temp = pd.read_csv(filename, sep='\t')

    df = df.join(temp, lsuffix='', rsuffix=str(s), how='outer')
    
    temp['rep'] = 'R'+str(s)
    temp['quant'] = s
    df2 = df2.append(temp, ignore_index=True, sort='False')

df = df.drop(['RetentionTime','Intensity'], axis=1)
df = df.dropna()

sns.set()
colors = sns.color_palette('muted')
plt.figure()

for i, s in enumerate(samples):
    sns.lineplot(df[f'RetentionTime{str(s)}'], df[f'Intensity{str(s)}'], 
            label =str(s), color = colors[i])
plt.legend()
plt.show()
plt.close()
for s in samples:
    print('%.4E' % Decimal(str(df[f'Intensity{str(s)}'].mean())), 
    '%.4E' % Decimal(str(integrate.simps(df[f'Intensity{str(s)}'], 
                                        df[f'RetentionTime{str(s)}']))))

plt.figure()
sns.lineplot('RetentionTime', 'Intensity', hue = 'rep', data = df2)  
plt.savefig('chro.svg', format='svg')
plt.show()
plt.close()
# =============================================================================
# print(df.Intensity.mean(),df2.Intensity.mean(),df3.Intensity.mean(),df4.Intensity.mean())
# print(df.Intensity.max(),df2.Intensity.max(),df3.Intensity.max(),df4.Intensity.max())
# 
# print('%.4E' % Decimal(str(integrate.simps(df.Intensity, df.RetentionTime))),
#       Decimal(integrate.simps(df2.Intensity, df2.RetentionTime)),
#       Decimal(integrate.simps(df3.Intensity, df3.RetentionTime)),
#       Decimal(integrate.simps(df4.Intensity, df4.RetentionTime)))
# =============================================================================
