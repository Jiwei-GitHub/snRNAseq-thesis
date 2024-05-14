import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import seaborn as sns
import scanpy.external as sce
import matplotlib.pyplot as plt
sc.set_figure_params(figsize=(4,4), dpi_save=600, dpi=300)


os.chdir("/mnt/h/Work_Home/Jiwei_ICMC_data/OvarianCancer/AdiposeTissue")



ad_merge = sc.read_h5ad("pydata/05-ad_merge_Major-Minor-raw.h5ad")
ad_merge.obs['Major_Minor'] = ad_merge.obs['Major_cell'].astype('str') + "+" + ad_merge.obs['Minor_cell'].astype('str')

df_kept = pd.read_excel("pydata/df.kept.xlsx")
Kept = pd.Series(df_kept['Kept'].values, index = df_kept['Major_Minor'].astype('str')).to_dict()
ad_merge.obs['Kept'] = ad_merge.obs['Major_Minor'].map(Kept).astype('category')



sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_merge-unfilter-0.3.png",
           legend_fontsize=8, legend_fontoutline=True)
sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_merge-unfilter-0.3-legend.png",
           legend_fontsize=8, legend_fontoutline=True)



sc.pl.umap(ad_merge[ad_merge.obs['Kept'] == "Yes"], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_merge-filter-0.3.png",
           legend_fontsize=8, legend_fontoutline=True)
sc.pl.umap(ad_merge[ad_merge.obs['Kept'] == "Yes"], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_merge-filter-0.3-legend.png",
           legend_fontsize=8, legend_fontoutline=True)






sc.tl.umap(ad_merge, min_dist=0.2)

sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_merge-unfilter-0.2.png",
           legend_fontsize=8, legend_fontoutline=True)
sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_merge-unfilter-0.2-legend.png",
           legend_fontsize=8, legend_fontoutline=True)


sc.pl.umap(ad_merge[ad_merge.obs['Kept'] == "Yes"], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_merge-filter-0.2.png",
           legend_fontsize=8, legend_fontoutline=True)
sc.pl.umap(ad_merge[ad_merge.obs['Kept'] == "Yes"], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_merge-filter-0.2-legend.png",
           legend_fontsize=8, legend_fontoutline=True)



sc.tl.umap(ad_merge, min_dist=0.1)

sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_merge-unfilter-0.1.png",
           legend_fontsize=8, legend_fontoutline=True)
sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_merge-unfilter-0.1-legend.png",
           legend_fontsize=8, legend_fontoutline=True)
           
           
sc.pl.umap(ad_merge[ad_merge.obs['Kept'] == "Yes"], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_merge-filter-0.1.png",
           legend_fontsize=8, legend_fontoutline=True)
sc.pl.umap(ad_merge[ad_merge.obs['Kept'] == "Yes"], color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_merge-filter-0.1-legend.png",
           legend_fontsize=8, legend_fontoutline=True)







ad_merge = ad_merge[ad_merge.obs['Kept'] == "Yes"]
ad_merge.write_h5ad("./pydata/05-ad_merge_Major-Minor-filter.h5ad")
ad_merge.obs.to_csv("./pydata/05-ad_merge_Major-Minor-filter.csv")




sc.tl.umap(ad_merge, min_dist=0.3)

sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_merge-filter-reUMAP_0.3.png",
           legend_fontsize=8, legend_fontoutline=True)
sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_merge-filter-reUMAP_0.3-legend.png",
           legend_fontsize=8, legend_fontoutline=True)



sc.tl.umap(ad_merge, min_dist=0.1)

sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_merge-filter-reUMAP_0.1.png",
           legend_fontsize=8, legend_fontoutline=True)
sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_merge-filter-reUMAP_0.1-legend.png",
           legend_fontsize=8, legend_fontoutline=True)



sc.tl.umap(ad_merge, min_dist=0.2)

sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_merge-filter-reUMAP_0.2.png",
           legend_fontsize=8, legend_fontoutline=True)
sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_merge-filter-reUMAP_0.2-legend.png",
           legend_fontsize=8, legend_fontoutline=True)


sc.tl.umap(ad_merge, min_dist=0.4)

sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, legend_loc = 'on data', 
           save="-ad_merge-filter-reUMAP_0.4.png",
           legend_fontsize=8, legend_fontoutline=True)
sc.pl.umap(ad_merge, color=["Major_cell_1", "Major_cell", "Minor_cell_1", "Minor_cell"], 
           ncols=2, 
           # legend_loc = 'on data', 
           save="-ad_merge-filter-reUMAP_0.4-legend.png",
           legend_fontsize=8, legend_fontoutline=True)