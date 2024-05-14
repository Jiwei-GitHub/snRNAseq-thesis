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

import loompy
df_meta = pd.read_excel("df.meta.xlsx")
ids = df_meta['Sample_id'].to_list()
ids
files = ["./velocity/looms/" + id + "/"+ id + ".loom" for id in ids]
loompy.combine(files=files,
               output_file="./velocity/Combined.loom")
ldata=scv.read("./velocity/Combined.loom")
ldata.write_h5ad("./velocity/Combined_raw.h5ad")
ldata.obs.to_csv("./velocity/Combined_raw.csv")
ldata

