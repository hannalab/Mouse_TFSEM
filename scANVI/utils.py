#####################
# Written by Dr. Noa Novershtern 
# https://github.com/hannalab/Mouse_TFSEM
#####################
import seaborn as sb
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

def plot_dist(adata,group_by='',query=''):
    
    items = adata.obs[group_by].unique().tolist()
    
    for d in (items): 
        sb.kdeplot(adata[adata.obs[group_by] == d].obs, x=query)
 
    plt.legend(items)
    plt.show()    
    
    
        



