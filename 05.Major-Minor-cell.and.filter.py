import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import seaborn as sns
import scanpy.external as sce
import matplotlib.pyplot as plt
sc.set_figure_params(figsize=(4,4), dpi_save=600, dpi=300)


# df_meta = pd.read_excel("df.meta.xlsx")
# ids = df_meta['Sample_id'].to_list()
# ids


if False:

    ids = ["Immune", "Stromal", "Epi_Tumor"]
    leiden=["leiden_0.6","leiden_0.4","leiden_0.2"]
    meta_data_sub_list = [None] * 3
    for i in range(0,3):
        id = ids[i]
        group = leiden[i]
        h5 = "pydata/04-ad-" + id + ".h5ad"
        adata = sc.read_h5ad(h5)
        

        cell_anno = pd.read_excel("cell_anno.xlsx", sheet_name= id + "-" + group)
        cell_type = pd.Series(cell_anno['Minor_cell'].values, index = cell_anno[group].astype('str')).to_dict()
        cell_type_1 = pd.Series(cell_anno['Minor_cell_1'].values, index = cell_anno[group].astype('str')).to_dict()
        
        adata.obs['Minor_cell'] = adata.obs[group].map(cell_type).astype('category')
        adata.obs['Minor_cell_1'] = adata.obs[group].map(cell_type_1).astype('category')
        
        del cell_anno
        del cell_type
        del cell_type_1
        
        sc.pl.umap(adata, color=['Major_cell_1', 'Major_cell','Minor_cell_1', 'Minor_cell'], 
                   save= "-" + id + '-CellTypes-raw.png', 
                   ncols=2,
                   legend_fontoutline=True,
                   legend_loc='on data', 
                   legend_fontsize=8)
        sc.pl.umap(adata, color=['Major_cell_1', 'Major_cell','Minor_cell_1', 'Minor_cell'], 
                   save= "-" + id + '-CellTypes-raw-1.png', 
                   ncols=2,
                   legend_fontoutline=True,
                   # legend_loc='on data', 
                   legend_fontsize=6)
        adata.write_h5ad("pydata/04-ad-" + id + "-Minor_cell-raw.h5ad")
        adata.obs.to_csv("pydata/04-ad-" + id + "-Minor_cell-raw.csv")
        
        meta_data_sub_list[i] = adata.obs










    ids = ["Immune", "Stromal", "Epi_Tumor"]
    leiden=["leiden_0.6","leiden_0.4","leiden_0.2"]
    meta_data_sub_list = [None] * 3
    for i in range(0,3):
        id = ids[i]
        df_tmp = pd.read_csv("pydata/04-ad-" + id + "-Minor_cell-raw.csv")
        meta_data_sub_list[i] = df_tmp
        
    ## concate meta_data_sub_list
    df_meta_data_sub = pd.concat(meta_data_sub_list[0:3])
    del meta_data_sub_list

    df_meta_data_sub.to_csv("pydata/df_meta_data_sub.csv")

    df_meta_data_sub = pd.read_csv("pydata/df_meta_data_sub.csv")
    ad_merge = sc.read_h5ad("pydata/03-ad_merge_Major_cell.h5ad")
    ad_merge = ad_merge[ad_merge.obs_names.isin(df_meta_data_sub['Cell_id'])]


    cell_type = pd.Series(df_meta_data_sub['Minor_cell'].values, index = df_meta_data_sub['Cell_id'].astype('str')).to_dict()
    cell_type_1 = pd.Series(df_meta_data_sub['Minor_cell_1'].values, index = df_meta_data_sub['Cell_id'].astype('str')).to_dict()

    ad_merge.obs['Minor_cell'] = ad_merge.obs_names.map(cell_type).astype('category')
    ad_merge.obs['Minor_cell_1'] = ad_merge.obs_names.map(cell_type_1).astype('category')

    ad_merge.write_h5ad("pydata/05-ad_merge_Major-Minor-raw.h5ad")
    ad_merge.obs.to_csv("pydata/05-ad_merge_Major-Minor-raw.csv")

    ad_merge.obs_names















