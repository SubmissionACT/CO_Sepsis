#!/usr/bin/env python
# coding: utf-8

# In[23]:


# Imports
import numpy as np
import pandas as pd
from causalnex.structure.dynotears import from_pandas_dynamic
from IPython.display import Image 

from causalnex.plots import plot_structure, NODE_STYLE, EDGE_STYLE

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


# In[19]:


sm=g_learnt
graph_attributes = {
    "splines": "spline",  
    "ordering": "out",
    "ratio": "fill",  
    "size": "4.5,3!",
    "fontcolor": "#FFFFFFD9",
    "fontname": "Helvetica",
    "fontsize": 100,
    "labeljust": "l",
    "labelloc": "t",
    "pad": "1,1",
    "dpi": 500,
    "nodesep": 0.8,
    "ranksep": ".5 equally",
    "bgcolor":"transparent",
}

# Making all nodes hexagonal with black coloring
node_attributes = {
    node: {
        "shape": "oval",
        "width": 5,
        "height": 1.5,
        "fillcolor": "#E1EDF3",##C9C0BB
        "penwidth": "1",
        "color": "#000000",
        "fontsize": 52,
        "labelloc": "t",
        "fontcolor": "#000000",
    }
    for node in sm.nodes
}

# Target nodes are colored differently
for node in sm.nodes:
    if "Sepsis_lag0" in node:  
        node_attributes[node]["fillcolor"] = "#F57A49"#DF5F00

# Customising edges
edge_attributes = {
    (u, v): {
        "penwidth": abs(w) * 15 + 4,  # Setting edge thickness
        "weight": int(5 * abs(w)),  # Higher "weight"s mean shorter edges
        "arrowsize": 2.5 - 2.0 * abs(w),  # Avoid too large arrows
        "arrowtail": "dot",
        "color": "#F9B86D" #CF793A
    }
    for u, v, w in sm.edges(data="weight")
}

for u, v, w in sm.edges(data="weight"):
    if w<0:  
        edge_attributes[(u, v)]["color"] = "#2A7EBC"

viz = plot_structure(
    sm.get_largest_subgraph(),
    prog="dot",
    graph_attributes=graph_attributes,
    node_attributes=node_attributes,
    edge_attributes=edge_attributes,
)

Image(viz.draw(format='png'))


# In[ ]:




