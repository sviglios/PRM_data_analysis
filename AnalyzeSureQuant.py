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
import matplotlib.pyplot as plt
import seaborn as sns

#sum ion intensity -> peptide signal
#normalize for heavy pep signal sample wise
#normalize for heave pep signal plate wise
#include peptides function for each protein. That is on the avg normalization on the peptide lvl
#exclude samples function on the protein signal lvl