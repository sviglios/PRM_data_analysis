# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:15:42 2020

@author: konka
"""

import pandas as pd
import numpy as np

df_in = pd.read_excel('Mapped_samples_SQruns.xlsx')

df_in['Condition'] = np.where(df_in['Sample name'].str.contains('B'), 'Acute', 'Impaired')

df_in['Timepoint'] = 0
for i in range(1,9):
    df_in.loc[df_in['Sample name'].str.contains(f'#{i}'), 'Timepoint'] = i
              
df_in['Type'] = np.where(df_in['Sample name'].str.contains('#1'), 'MEC', 'HYC')
     
df_in.to_excel('Annotation_samples.xlsx', index=False)