#!/usr/bin/env python
# coding: utf-8

# In[5]:


# Imports
import numpy as np
import pandas as pd
from causalnex.structure.dynotears import from_pandas_dynamic
from IPython.display import Image 

# Data
data_all = pd.read_csv("sepsis_TTTF.csv")
data_all = data_all.loc[:, ["sep","co","temp","rh"]]
data_all.rename(columns={'sep': "Sepsis",
                           'co': "CO",
                         'temp': 'Temp',
                         'rh': 'Humid'}, inplace=True)
data_all
df = data_all.to_numpy()

# Run DYNOTEARS
g_learnt,w,a = from_pandas_dynamic(data_all,p=3,lambda_w=0.05,lambda_a=0.05,w_threshold=0.01)
pd.DataFrame(g_learnt.edges(data="weight"))


# In[ ]:




