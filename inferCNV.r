options(stringsAsFactors = F)
library(tidyverse)
library(Seurat)
library(patchwork)
# library(harmony)
# library(SeuratData)
library(future, warn.conflicts = T)
plan("multicore", workers = 6)
options(future.globals.maxSize = 18 * 1024 * 1024^2)

### Run inferCNV
library(infercnv)

options(scipen = 100)

if(T){library(ggsci)
  colorseq = c(pal_npg("nrc")(10),
               pal_jco("default", alpha = 0.6)(10),
               pal_futurama("planetexpress", alpha = 0.6)(10))
  mycolor.1 <- RColorBrewer::brewer.pal(8, 'Set2')
  mycolor.2 <- RColorBrewer::brewer.pal(8, 'Accent')
  mycolor.3 <- RColorBrewer::brewer.pal(8, 'Dark2')
  
  
  colorseq = c(
    '#1f77b4',
    '#ff7f0e',
    '#279e68',
    '#d62728',
    '#aa40fc',
    '#8c564b',
    '#e377c2',
    '#b5bd61',
    '#17becf',
    '#aec7e8',
    '#ffbb78',
    '#98df8a',
    '#ff9896',
    '#c5b0d5',
    '#c49c94',
    '#f7b6d2',
    '#dbdb8d',
    mycolor.1,
    mycolor.2,
    mycolor.3,
    colorseq
  )
}  




# see scripts on server
if(F){
if(F){
  
  
  sc <- qs::qread("../05-ad_merge_Major-Minor-filter.h5ad/sc.qs")
  sc$to_Tumor <- str_split(sc$Tissue_Site, pattern = "-", simplify = T)[,2]
  
  sc$Minor_cell %>% unique()
  sc$Tissue_Site %>% unique()
  
  sc <- sc %>% subset(Minor_cell %in% c("Epi_Tumor-Cycling", "Epi_Tumor"))
  
  cnt.cnv <- sc@assays$RNA@counts %>% as.matrix()
  
  anno.file <- sc@meta.data %>% 
    # rownames_to_column('cell_id') %>% 
    mutate(Normal = ifelse(Tissue_Site == "Om-dis", "Normal", "Malignent")) %>% 
    select(Cell_id, Normal, Minor_cell)
  write.table(anno.file, file = "anno.file.01.tsv", quote = F,sep = '\t', col.names = F, row.names = F)
  
  
  infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix=cnt.cnv,
    annotations_file='anno.file.01.tsv',
    delim="\t",
    gene_order_file="h://Rdata/genome_index/Ki_D/gencode.V37.MCpyV.gene.name.position.txt",
    # ref_group_names=c("Immune", "Endothelial", "Fibroblast")
    ref_group_names = "Normal"
  )
  
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff= 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir="infercnv.01/",  # dir is auto-created for storing outputs
                               
                               cluster_by_groups= T,   # cluster
                               denoise=T,
                               HMM=T,
                               num_threads = 16
                               
  )
  
}
#  save(infercnv_obj, file = "infercnv_obj.01.rdata")
if(F){
  
  setwd("H:/Work_Home/Jiwei_ICMC_data/OvarianCancer/AdiposeTissue/plotting.in.r/inferCNV")
  
  sc <- qs::qread("../05-ad_merge_Major-Minor-filter.h5ad/sc.qs")
  sc$to_Tumor <- str_split(sc$Tissue_Site, pattern = "-", simplify = T)[,2]
  
  sc$Minor_cell %>% unique()
  sc$Tissue_Site %>% unique()
  
  sc <- sc %>% subset(Minor_cell %in% c("Epi_Tumor-Cycling", "Epi_Tumor"))
  
  cnt.cnv <- sc@assays$RNA@counts %>% as.matrix()
  
  anno.file <- sc@meta.data %>% 
    # rownames_to_column('cell_id') %>% 
    mutate(Normal = ifelse(Tissue_Site == "Om-dis", "Normal", "Malignent")) %>% 
    select(Cell_id, Normal, Minor_cell)
  write.table(anno.file, file = "anno.file.01.tsv", quote = F,sep = '\t', col.names = F, row.names = F)
  
  
  infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix=cnt.cnv,
    annotations_file='anno.file.01.tsv',
    delim="\t",
    gene_order_file="h://Rdata/genome_index/Ki_D/gencode.V37.MCpyV.gene.name.position.txt",
    # ref_group_names=c("Immune", "Endothelial", "Fibroblast")
    ref_group_names = NULL
  
  )
  
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff= 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir="infercnv.02/",  # dir is auto-created for storing outputs
                               
                               cluster_by_groups= T,   # cluster
                               denoise=T,
                               HMM=T,
                               num_threads = 16
                               
  )
  
}
# save(infercnv_obj, file = "infercnv_obj.02_Null.rdata")
}

setwd("H:/Work_Home/Jiwei_ICMC_data/OvarianCancer/AdiposeTissue/plotting.in.r/inferCNV")
infer.obj <- qs::qread("infercnv_obj.03.qs")
# infer.obj <- qs::qread("infercnv_obj.04Null.qs")



eMat <- infer.obj@expr.data

eMat <- eMat-1


eMat <- apply(eMat, 2, function(x){
  lapply(x, function(i)i*i) %>% unlist %>% mean
  
})
df.meta <- data.table::fread("../05-ad_merge_Major-Minor-filter.h5ad/meta.data..new.csv")
df.meta.epi <- data.table::fread("../../pydata/04-ad-Epi_Tumor-Minor_cell-filter.new.csv")
df.infercluster <- data.table::fread("infercnv.04Null/infercnv.observation_groupings.txt", sep = ' ') %>% rename(Cell_id = V1)

df <- as.data.frame(eMat) %>% 
  `colnames<-`('CNV_score') %>% 
  rownames_to_column('Cell_id') %>% 
  left_join(df.meta %>% select( - contains("leiden"))) %>% 
  mutate(Normal = ifelse(Tissue_Site == "Dis", "Normal", "Malignent")) %>% 
  left_join(df.meta.epi %>% select(Cell_id, contains("leiden"))) %>% 
  left_join(df.infercluster)
  


ggplot(df)+
  theme_bw()+
  scale_fill_manual(values = colorseq)+
  scale_color_manual(values = colorseq)+
  
  geom_violin(aes(x= factor(`leiden_0.05`), y=CNV_score, fill=factor(`leiden_0.05`)))+
  
  ylim(c(0,0.01))+


ggplot(df)+
  theme_bw()+
  scale_fill_manual(values = colorseq)+
  scale_color_manual(values = colorseq)+
  
  geom_violin(aes(x= factor(`leiden_0.1`), y=CNV_score, fill=factor(`leiden_0.1`)))+
  
  ylim(c(0,0.01))+


ggplot(df)+
  theme_bw()+
  scale_fill_manual(values = colorseq)+
  scale_color_manual(values = colorseq)+
  
  geom_violin(aes(x= factor(`leiden_0.2`), y=CNV_score, fill=factor(`leiden_0.2`)))+
  
  ylim(c(0,0.01))+
  


ggplot(df)+
  theme_bw()+
  scale_fill_manual(values = colorseq)+
  scale_color_manual(values = colorseq)+
  
  geom_violin(aes(x= Minor_cell_2, y=CNV_score, fill=Minor_cell_2))+
  ylim(c(0,0.01))+
  
 ggplot(df)+
  theme_bw()+
  scale_fill_manual(values = colorseq)+
  scale_color_manual(values = colorseq)+
  
  geom_violin(aes(x= factor(`Annotation Group`), y=CNV_score, fill=factor(`Annotation Group`)))+
  ylim(c(0,0.01))+
  
  ggplot(df)+
  theme_bw()+
  scale_fill_manual(values = colorseq)+
  scale_color_manual(values = colorseq)+
  
  geom_violin(aes(x= Normal, y=CNV_score, fill=Normal))+
  ylim(c(0,0.01))+
  
  plot_layout(guides = 'collect',design = "1112222
              3333344
              5566777")
  







sc.raw <- qs::qread("05-ad_merge_Major-Minor-filter.qs")
sc <- subset(sc.raw, Cell_id %in% df$Cell_id)

df.meta.epi <- data.table::fread("../../pydata/04-ad_sub_epiTumor-before_reUMAP-meta_data.csv")

df.meta.epi <- data.frame(Cell_id = sc$Cell_id) %>% 
  left_join(df.meta.epi)


all(df.meta.epi$Cell_id == colnames(sc))

sc[['UMAP']] <- CreateDimReducObject(embeddings = as.matrix(data.frame(UMAP_1 = df.meta.epi$UMAP1, UMAP_2 = df.meta.epi$UMAP2, 
                                                                       row.names = df.meta.epi$Cell_id)),
                                     key = "UMAP_", global = T, assay = "RNA")

meta.data <- sc@meta.data %>% 
  rownames_to_column("x") %>% 
  left_join(df %>% select(Cell_id, `Annotation Group`, Normal)) %>% 
  as.data.frame %>% 
  column_to_rownames("x")
sc@meta.data <- meta.data

colnames(meta.data)

DimPlot(sc, reduction = "UMAP", group.by = c("Annotation Group", "Normal"))&
  scale_color_manual(values = colorseq)











# eMat <- eMat %>% t %>% scale() %>% t

cnv <- apply(eMat,2,mean) %>% 
  as.data.frame() %>% 
  `colnames<-`("CNV.Score") %>% 
  rownames_to_column('Cell_id')

cnv <- eMat %>% 
  t %>% 
  as.data.frame() %>% 
  rownames_to_column('Cell_id') %>% 
  gather(CNV.Score, val, (2:nrow(eMat)))


infer.obj@tumor_subclusters %>% unlist %>% table



df.meta <- data.table::fread("../05-ad_merge_Major-Minor-filter.h5ad/meta.data..new.csv")
df.meta.epi <- data.table::fread("../../pydata/04-ad-Epi_Tumor-Minor_cell-filter.new.csv")


df.cnv <- cnv %>% 
  left_join(df.meta %>% select( - contains("leiden"))) %>% 
  mutate(Normal = ifelse(Tissue_Site == "Dis", "Normal", "Malignent")) %>% 
  left_join(df.meta.epi %>% select(Cell_id, contains("leiden")))
  



ggplot(df.cnv)+
  theme_classic()+
  
  geom_violin(aes(Normal, CNV.Score, fill=Normal))




ggplot(df.cnv )+
  theme_classic()+
  
  geom_violin(aes(`leiden_0.05`, CNV.Score-1, fill=factor(`leiden_0.05`)))


ggplot(df.cnv )+
  theme_classic()+
  
  geom_violin(aes(`leiden_0.1`, CNV.Score-1, fill=factor(`leiden_0.1`)))


ggplot(df.cnv )+
  theme_classic()+
  
  geom_violin(aes(`leiden_0.2`, CNV.Score-1, fill=factor(`leiden_0.2`)))



ggplot(df.cnv )+
  theme_classic()+
  
  geom_violin(aes(Minor_cell_2, (CNV.Score-1), fill=Minor_cell_2))






### Use add_to_seurat
cnv.meta <- data.table::fread("x.meta.data.csv")
df.meta.epi <- data.table::fread("../../pydata/04-ad-Epi_Tumor-Minor_cell-filter.new.csv")

colnames(cnv.meta)


df <- cnv.meta %>% 
  select(Cell_id, contains("infercnv_")) %>% 
  left_join(df.meta.epi) 

ggdf <- df %>% 
  # gather(infer, Val, contains("infercnv_proportion_scaled_cnv_chr"))
  # gather(infer, Val, contains("infercnv_top"))
  # gather(infer, Val, contains("infercnv_has_cnv_chr"))
  # gather(infer, Val, contains("infercnv_has_loss_chr"))
gather(infer, Val, contains("infercnv_proportion_cnv_chr"))


ggplot(ggdf )+
  theme_classic()+
  geom_violin(aes(Minor_cell_2, Val, fill=Minor_cell_2))+
  
  ggplot(ggdf )+
  theme_classic()+
  geom_violin(aes(`leiden_0.05`, Val, fill=factor(`leiden_0.05`)))+
  
  ggplot(ggdf )+
  theme_classic()+
  geom_violin(aes(`leiden_0.1`, Val, fill=factor(`leiden_0.1`)))+
  
  ggplot(ggdf )+
  theme_classic()+
  geom_violin(aes(`leiden_0.2`, Val, fill=factor(`leiden_0.2`)))+
  
  plot_layout(ncol = 2)










#'@
#'@
#'@
#'#'@
#'#'@
#'#'@
#'#'@
#'#'@
#'#'@
#'#'@
#'#'@ add annotation cluster to pydata

meta.data <- data.table::fread("../../pydata/04-ad_sub_epiTumor-before_reUMAP-meta_data.csv", sep = ",")
df.infercluster <- data.table::fread("infercnv.04Null/infercnv.observation_groupings.txt", sep = ' ') %>% 
  rename(Cell_id = V1) %>% 
  mutate(infercnv_cluster = `Annotation Group`)



meta.data <- meta.data %>% 
  left_join(df.infercluster %>% select(Cell_id, infercnv_cluster)) %>% 
  as.data.frame %>% 
  column_to_rownames("V1")

data.table::fwrite(meta.data, file = "../../pydata/04-ad_sub_epiTumor-before_reUMAP-meta_data.infercnv_cluster.csv",
                   sep = ",", row.names = T, col.names = T, quote = F
                   )















