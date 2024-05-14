import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import seaborn as sns
import scanpy.external as sce
import matplotlib.pyplot as plt
sc.set_figure_params(figsize=(4,4), dpi_save=600, dpi=300)


df_meta = pd.read_excel("df.meta.xlsx")
ids = df_meta['Sample_id'].to_list()
ids

adata = sc.read_h5ad("./pydata/03-ad.h5ad")



## Read in Major_cells
cell_anno = pd.read_excel("cell_anno.xlsx", sheet_name='major_0.2')
cell_type = pd.Series(cell_anno['Major_cell'].values, index = cell_anno['leiden_0.2'].astype('str')).to_dict()
cell_type_1 = pd.Series(cell_anno['Major_cell_1'].values, index = cell_anno['leiden_0.2'].astype('str')).to_dict()

# adata = sc.read_h5ad("pydata/03-ad.h5ad")
adata.obs['Major_cell'] = adata.obs['leiden_0.2'].map(cell_type).astype('category')
adata.obs['Major_cell_1'] = adata.obs['leiden_0.2'].map(cell_type_1).astype('category')
adata.write_h5ad("pydata/03-ad_merge_Major_cell.h5ad")
adata.obs.to_csv("pydata/03-ad_merge_Major_cell.csv")


## split adata based on Major_cell_1 &&&  Read in Major_cells to raw.h5ad
## 
# meta_data = sc.read_h5ad("pydata/03-ad_merge_Major_cell.h5ad")
# meta_data = meta_data.obs
meta_data = adata.obs
adata = sc.read_h5ad("pydata/raw.h5ad")
adata.obs = meta_data


# split = adata.obs.groupby('Major_cell_1').indices
# adata_list = [adata[Major_cell_1] for Major_cell_1 in split.values()]




## ad_sub re-cluster to h5
ids = ["Stromal", "Immune", "Epi_Tumor"]
for i in range(3):
    print(i)
    id=ids[i]
    print(id)
    ad = adata[adata.obs.Major_cell_1 == id]
    ad.raw = ad
    # id = ad.obs.Major_cell_1.astype("str").to_list()[0]
    print(id)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, flavor='seurat', n_top_genes = 5000, subset=False, inplace=True, batch_key='orig.ident')
    
    sc.pp.regress_out(ad, ['nCount_RNA', 'MT_pct'])
    sc.pp.scale(ad, max_value=10)
    
    sc.tl.pca(ad, svd_solver='arpack')
    sce.pp.harmony_integrate(adata = ad,  max_iter_harmony=100,
                             key = "orig.ident", 
                             basis='X_pca', adjusted_basis='X_pca_harmony')
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=50, use_rep = "X_pca_harmony")
    sc.tl.umap(ad, min_dist=0.3)
    sc.tl.leiden(ad, resolution=0.05, key_added="leiden_0.05")
    sc.tl.leiden(ad, resolution=0.1, key_added="leiden_0.1")
    sc.tl.leiden(ad, resolution=0.2, key_added="leiden_0.2")
    sc.tl.leiden(ad, resolution=0.4, key_added="leiden_0.4")
    sc.tl.leiden(ad, resolution=0.6, key_added="leiden_0.6")
    sc.tl.leiden(ad, resolution=0.8, key_added="leiden_0.8")
    # sc.tl.leiden(ad, resolution=0.8, key_added="leiden_0.8")

    ad.write_h5ad("pydata/04-ad-" + id + ".h5ad")
    ad.obs.to_csv("pydata/04-ad-" + id + "-meta_data.csv")
    
    sc.pl.umap(ad, color=["leiden_0.05","leiden_0.1", "leiden_0.2",
                             "leiden_0.4" ,"leiden_0.6", "leiden_0.8",
                             'orig.ident', 'Tissue_Site','Sample_Type',
                             'MT_pct','Ribo_pct' ,'log2nCount'], 
               ncols=6, legend_loc = 'on data', 
               save="-"+id+".png",
               legend_fontsize=8, legend_fontoutline=True)
    sc.pl.umap(ad, color=["leiden_0.05","leiden_0.1", "leiden_0.2",
                             "leiden_0.4" ,"leiden_0.6", "leiden_0.8",
                             'orig.ident', 'Tissue_Site','Sample_Type',
                             'MT_pct','Ribo_pct' ,'log2nCount'], 
               ncols=6, 
               # legend_loc = 'on data', 
               save="-"+id+"-legend.png",
               legend_fontsize=8, legend_fontoutline=True)
    
    c_marker = pd.read_excel('./Markers_new.xlsx', sheet_name='Adipose')
    # intersection = set(adata.var_names,c_marker['gene_name'])
    c_marker = c_marker[c_marker['gene_name'].isin(ad.raw.var_names)]
    
    sc.pl.umap(ad, color=c_marker['gene_name'], 
    ncols=6, 
    # legend_loc = 'on data', 
    save="-"+id+"-marker0-legend.png",
    legend_fontsize=8, legend_fontoutline=True)
               
   
    
    ## dotplot for markers_1-Adipose
    c_marker = pd.read_excel('./Markers_new.xlsx', sheet_name='Adipose')
    # intersection = set(adata.var_names,c_marker['gene_name'])
    c_marker = c_marker[c_marker['gene_name'].isin(ad.raw.var_names)]

    # markers_1 = pd.Series(c_marker['gene_name'].values, index = c_marker['Cell_type'].astype('str'))
    markers_1 = c_marker.groupby('Cell_type').gene_name.agg(list).to_dict()
    markers_1

    sc.pl.dotplot(ad, markers_1, ['leiden_0.2'], use_raw=None, save="04-ad-" + id + "adipose.marker-leiden_0.2.png")
    sc.pl.dotplot(ad, markers_1, ['leiden_0.4'], use_raw=None, save="04-ad-" + id + "adipose.marker-leiden_0.4.png")
    sc.pl.dotplot(ad, markers_1, ['leiden_0.6'], use_raw=None, save="04-ad-" + id + "adipose.marker-leiden_0.6.png")
    sc.pl.dotplot(ad, markers_1, ['leiden_0.8'], use_raw=None, save="04-ad-" + id + "adipose.marker-leiden_0.8.png")


    ## dotplot for markers_2
    c_marker = pd.read_excel('./Markers_new.xlsx', sheet_name='Sheet2')
    c_marker = c_marker[c_marker['gene_name'].isin(ad.raw.var_names)]
    markers_2 = c_marker.groupby('Cell_type').gene_name.agg(list).to_dict()

    sc.pl.dotplot(ad, markers_2, ['leiden_0.2'], use_raw=None, save="04-ad-" + id + "adipose.marker_2-leiden_0.2.png")
    sc.pl.dotplot(ad, markers_2, ['leiden_0.4'], use_raw=None, save="04-ad-" + id + "adipose.marker_2-leiden_0.4.png")
    sc.pl.dotplot(ad, markers_2, ['leiden_0.6'], use_raw=None, save="04-ad-" + id + "adipose.marker_2-leiden_0.6.png")
    sc.pl.dotplot(ad, markers_2, ['leiden_0.8'], use_raw=None, save="04-ad-" + id + "adipose.marker_2-leiden_0.8.png")


    c_marker = pd.read_excel('./Markers_new.xlsx', sheet_name='Adipose_1')
    c_marker = c_marker[c_marker['gene_name'].isin(ad.raw.var_names)]
    markers_3 = c_marker.groupby('Cell_type').gene_name.agg(list).to_dict()

    sc.pl.dotplot(ad, markers_3, ['leiden_0.2'], use_raw=None, save="04-ad-" + id + "adipose.marker_3-leiden_0.2.png")
    sc.pl.dotplot(ad, markers_3, ['leiden_0.4'], use_raw=None, save="04-ad-" + id + "adipose.marker_3-leiden_0.4.png")
    sc.pl.dotplot(ad, markers_3, ['leiden_0.6'], use_raw=None, save="04-ad-" + id + "adipose.marker_3-leiden_0.6.png")
    sc.pl.dotplot(ad, markers_3, ['leiden_0.8'], use_raw=None, save="04-ad-" + id + "adipose.marker_3-leiden_0.8.png")
    
    
    group = 'leiden_0.2'
    sc.tl.rank_genes_groups(ad, groupby=group, method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(ad, n_genes=10, groupby=group,
                                   save="-deg-"+id+"-"+group+".png")
                                   
    group = 'leiden_0.4'
    sc.tl.rank_genes_groups(ad, groupby=group, method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(ad, n_genes=10, groupby=group,
                                   save="-deg-"+id+"-"+group+".png")    
    
    group = 'leiden_0.6'
    sc.tl.rank_genes_groups(ad, groupby=group, method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(ad, n_genes=10, groupby=group,
                                   save="-deg-"+id+"-"+group+".png") 

    
    del ad
    del group



































