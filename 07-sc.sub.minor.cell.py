import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import seaborn as sns
import scanpy.external as sce
import matplotlib.pyplot as plt
sc.set_figure_params(figsize=(4,4), dpi_save=600, dpi=300)


os.chdir("/mnt/h/Work_Home/Jiwei_ICMC_data/OvarianCancer/AdiposeTissue/")


meta_data = pd.read_csv("./pydata/05-ad_merge_Major-Minor-filter.csv")

ids = ["Stromal", "Immune", "Epi_Tumor"]
for i in range(3):
    print(i)
    id=ids[i]
    print(id)

    adata = sc.read_h5ad("pydata/04-ad-" + id + "-Minor_cell-raw.h5ad")

    
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-unfilter-0.3.png",
           legend_fontsize=8, legend_fontoutline=True)
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-unfilter-0.3-legend.png",
           legend_fontsize=8, legend_fontoutline=True)

    sc.pl.umap(adata[adata.obs_names.isin(meta_data['Cell_id'])], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-0.3.png",
           legend_fontsize=8, legend_fontoutline=True)
    sc.pl.umap(adata[adata.obs_names.isin(meta_data['Cell_id'])], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-0.3-legend.png",
           legend_fontsize=8, legend_fontoutline=True)



    
    sc.tl.umap(adata, min_dist=0.2)
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-unfilter-0.2.png",
           legend_fontsize=8, legend_fontoutline=True)
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-unfilter-0.2-legend.png",
           legend_fontsize=8, legend_fontoutline=True)

    sc.pl.umap(adata[adata.obs_names.isin(meta_data['Cell_id'])], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-0.2.png",
           legend_fontsize=8, legend_fontoutline=True)
    sc.pl.umap(adata[adata.obs_names.isin(meta_data['Cell_id'])], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-0.2-legend.png",
           legend_fontsize=8, legend_fontoutline=True)


    sc.tl.umap(adata, min_dist=0.1)
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-unfilter-0.1.png",
           legend_fontsize=8, legend_fontoutline=True)
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-unfilter-0.1-legend.png",
           legend_fontsize=8, legend_fontoutline=True)

    sc.pl.umap(adata[adata.obs_names.isin(meta_data['Cell_id'])], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-0.1.png",
           legend_fontsize=8, legend_fontoutline=True)
    sc.pl.umap(adata[adata.obs_names.isin(meta_data['Cell_id'])], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-0.1-legend.png",
           legend_fontsize=8, legend_fontoutline=True)




    


    adata = adata[adata.obs_names.isin(meta_data['Cell_id'])]
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-reUMAP-0.3.png",
           legend_fontsize=8, legend_fontoutline=True)
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-reUMAP-0.3-legend.png",
           legend_fontsize=8, legend_fontoutline=True)

    sc.tl.umap(adata, min_dist=0.1)
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-reUMAP-0.1.png",
           legend_fontsize=8, legend_fontoutline=True)
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-reUMAP-0.1-legend.png",
           legend_fontsize=8, legend_fontoutline=True)

    sc.tl.umap(adata, min_dist=0.2)
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-reUMAP-0.2.png",
           legend_fontsize=8, legend_fontoutline=True)
    sc.pl.umap(adata, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_sub-"+ id + "-filter-reUMAP-0.2-legend.png",
           legend_fontsize=8, legend_fontoutline=True)




