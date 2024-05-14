import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import seaborn as sns
import scanpy.external as sce
import matplotlib.pyplot as plt
import scvelo as scv
sc.set_figure_params(figsize=(4,4), dpi_save=600, dpi=300)



os.chdir("/mnt/h/Work_Home/Jiwei_ICMC_data/OvarianCancer/AdiposeTissue/")

adata = sc.read_h5ad("velocity/Combined_raw-filter_with_meta.h5ad")

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata, n_jobs=20)
scv.tl.velocity(adata, mode='dynamical', n_jobs=18)
scv.tl.velocity_graph(adata, n_jobs=18)
adata.write_h5ad('velocity/ladata-all.h5ad')




lad_adipocyte = adata[adata.obs['Minor_cell'] == 'Adipocyte']
scv.pp.filter_and_normalize(lad_adipocyte, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(lad_adipocyte, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(lad_adipocyte, n_jobs=20)
scv.tl.velocity(lad_adipocyte, mode='dynamical', n_jobs=18)
scv.tl.velocity_graph(lad_adipocyte, n_jobs=18)
lad_adipocyte.write_h5ad('velocity/lad_adipocyte.h5ad')

del lad_adipocyte


lad_Stromal = adata[adata.obs['Major_cell_1'] == 'Stromal']
scv.pp.filter_and_normalize(lad_Stromal, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(lad_Stromal, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(lad_Stromal, n_jobs=20)
scv.tl.velocity(lad_Stromal, mode='dynamical', n_jobs=18)
scv.tl.velocity_graph(lad_Stromal, n_jobs=18)
lad_Stromal.write_h5ad('velocity/lad_Stromal.h5ad')

del lad_Stromal








