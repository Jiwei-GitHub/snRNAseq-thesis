options(string.as.factor = F)
# library(Seurat)
library(tidyverse)
library(patchwork)

library(future)
plan("multicore", workers = 4, future.seed=TRUE)
options(future.globals.maxSize= 20*1024*1024*1024, future.seed = TRUE)


library(ggsci)

## color
if(T){
  library(ggsci)
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


library(monocle)
setwd("H:/Work_Home/Jiwei_ICMC_data/OvarianCancer/AdiposeTissue/plotting.in.r")
# library(monocle3)




sc <- qs::qread("./05-ad_merge_Major-Minor-filter.h5ad/05-ad_merge_Major-Minor-filter.qs")
load("h://Rdata/19_geneIDtransform/gene.gtf.gencode.V37.MCPyV.rdata")

if(T){
sc$Major_cell_1 %>% unique()

# sc.sub <- sc %>% subset(Major_cell == "Adipocyte" & Tissue_Site != "Ovary")
sc.sub <- sc %>% subset(Major_cell_1 %in% c("Stromal") & Tissue_Site != "Ovary")

expCnt <- sc.sub@assays$RNA@counts
expressed_genes <- apply(expCnt, 1, function(x) sum(x!=0)) 
expCnt <- expCnt[expressed_genes>=50,]
cellMeta <- AnnotatedDataFrame(sc.sub@meta.data)
geneMeta <- data.frame(gene_short_name = rownames(expCnt), row.names = rownames(expCnt)) 
geneMeta <- AnnotatedDataFrame(data = geneMeta)

cds <- newCellDataSet(cellData = expCnt, 
                      phenoData = cellMeta, 
                      featureData = geneMeta,
                      expressionFamily=negbinomial.size())
rm(expCnt, geneMeta, cellMeta)
gc()



cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


diff_test_res <- differentialGeneTest(cds, cores = 10,
                                      fullModelFormulaStr = "~ Tissue_Site")
openxlsx::write.xlsx(diff_test_res, file = "monocle/diff_test_res-APSC-Om-Tissue_Site.xlsx")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

cds_ordering <- setOrderingFilter(cds = cds, ordering_genes = ordering_genes)
cds_ordering <- reduceDimension(cds_ordering, max_components = 2,
                                method = 'DDRTree')
qs::qsave(cds_ordering, 'monocle/cds_ordering.qs')

cds_ordering <- orderCells(cds_ordering)


plot_cell_trajectory(cds_ordering, color_by = "Tissue_Site")+
  plot_cell_trajectory(cds_ordering, color_by = "Minor_cell")
  # plot_cell_trajectory(cds_ordering, color_by = "State")+
  plot_cell_trajectory(cds_ordering, color_by = "Pseudotime")



diff_test_Pseudotime <- differentialGeneTest(cds_ordering, cores = 10,
                                             fullModelFormulaStr = "~sm.ns(Pseudotime)")
openxlsx::write.xlsx(diff_test_Pseudotime, file = "monocle/diff_test_Pseudotime.xlsx")

my_genes <- row.names(subset(fData(cds_ordering), gene_short_name %in% c("FABP4","MXRA8","ADIPOQ"))) 
plot_genes_in_pseudotime(cds_subset = cds_ordering[my_genes,], min_expr = 0.1, color_by="Minor_cell")

}






cds <- qs::qread("monocle/APSE-Om and Adiopocyte/cds_ordering.qs")


diff_tissue <- readxl::read_excel("monocle/APSE-Om and Adiopocyte/diff_test_res-APSC-Om-Tissue_Site.xlsx") %>% 
  left_join(gene.gtf %>% mutate(gene_short_name=gene_name))
diff_time <- readxl::read_excel("monocle/APSE-Om and Adiopocyte/diff_test_Pseudotime.xlsx") %>% 
  left_join(gene.gtf %>% mutate(gene_short_name=gene_name) %>% filter(!duplicated(gene_name)))


df.mono <- plot_cell_trajectory(cds, color_by = "Pseudotime", show_backbone = T, cell_size = 0.25)$data
data.table::fwrite(df.mono, file = "monocle/APSE-Om and Adiopocyte/df.mono.csv")
##


df.mono <- data.table::fread("monocle/APSE-Om and Adiopocyte/df.mono.csv")
df.meta <- data.table::fread("../pydata/04-ad_sub_adi-meta_data.csv")


df.mono <- df.mono %>% 
  dplyr::select(Cell_id, data_dim_1, data_dim_2, Pseudotime,State) %>% 
  filter(Cell_id %in% df.meta$Cell_id) %>% 
  left_join(df.meta) %>% 
  mutate(cluster = factor(cluster)) %>% 
  mutate(`leiden_0.4` = factor(`leiden_0.4`))
  
  
if(T){
p.1 <- ggplot(df.mono, aes(data_dim_1, data_dim_2))+
  xlab("Component1")+
  ylab("Component2")+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.background = element_blank())+
  scale_fill_manual(values = colorseq)+
  scale_color_manual(values = colorseq)+
  geom_point(aes(fill=Tissue_Site), shape = 21, size = 1, color=alpha("grey40", 0.3))+
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(n.breaks = 3)+
  
  theme(legend.position = c(0.8,0.8),
        legend.title = element_blank())+
  guides(fill = guide_legend(override.aes = list(size=3)))


p.2 <- ggplot(df.mono %>% mutate(cluster=factor(cluster)), aes(data_dim_1, data_dim_2))+
  xlab("Component1")+
  ylab("Component2")+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.background = element_blank())+
  scale_fill_manual(values = colorseq)+
  scale_color_manual(values = colorseq)+
  geom_point(aes(fill=`leiden_0.4`), shape = 21, size = 1, color=alpha("grey40", 0.3))+
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(n.breaks = 3)+
  
  theme(legend.position = c(0.8,0.8),
        legend.title = element_blank())+
  guides(fill = guide_legend(override.aes = list(size=3)))

p.3 <- ggplot(df.mono %>% mutate(cluster=factor(cluster)), aes(data_dim_1, data_dim_2))+
  xlab("Component1")+
  ylab("Component2")+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.background = element_blank())+
  # scale_fill_manual(values = colorseq)+
  scale_fill_gradient2()+
  scale_color_manual(values = colorseq)+
  geom_point(aes(fill=Pseudotime), shape = 21, size = 1, color=alpha("grey40", 0.3))+
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(n.breaks = 3)+
  
  theme(legend.position = c(0.8,0.7),
        legend.title = element_blank())


}
p <- p.1 + p.2 + p.3 + plot_layout(ncol = 3)
p
ggsave(plot = p, filename = "monocle/APSE-Om and Adiopocyte/plots/trajectory.leiden0.4.tiff",
       height = 8, width = 24,
       dpi = 600, units = "cm")


ggplot(df.mono # %>% 
         # mutate(`leiden_0.4` = as.integer(`leiden_0.4`))
       )+
  scale_fill_manual(values = colorseq)+
  theme_classic(base_line_size = 0.25)+
  # geom_violin(aes(`leiden_0.4`, Pseudotime, fill = cluster), alpha=0.6)+
  geom_boxplot(aes(factor(`leiden_0.4`), Pseudotime, fill = factor(`leiden_0.4`)), alpha=0.6, outlier.colour = NA)+
  # geom_jitter(aes(`leiden_0.4`, Pseudotime, fill = cluster), alpha=0.6, shape=21, height = 0)+
  stat_smooth(aes(as.integer(`leiden_0.4`), Pseudotime), geom = 'line', method = "glm", color='grey', linewidth=1)










ordering_genes <- row.names(subset(diff_tissue, qval < 0.01))

cds_ordering <- setOrderingFilter(cds = cds, ordering_genes = ordering_genes)
cds_ordering <- reduceDimension(cds_ordering, max_components = 2,
                                method = 'DDRTree')

cds_ordering <- orderCells(cds_ordering)




p <- plot_cell_trajectory(cds, color_by = "Tissue_Site", show_backbone = T, cell_size = 0.25)+
  plot_cell_trajectory(cds, color_by = "leiden_0.6", show_backbone = T, cell_size = 0.25)+
  plot_cell_trajectory(cds, color_by = "State", show_backbone = T, cell_size = 0.25)+
  plot_cell_trajectory(cds, color_by = "Pseudotime", show_backbone = T, cell_size = 0.25)+
  scale_color_gradient(low = "yellow", high = "red")
p
ggsave(plot = p, 
       filename = "monocle/APSE-Om and Adiopocyte/plots/cds_ordering-subStroma-Minor_cell.tiff", 
       width = 18, height=18, dpi = 300,
       units = 'cm')

plot_cell_trajectory(cds, color_by = "orig.ident")


genes <- diff_time %>% 
  filter(qval<=1e-6) %>% 
  filter(gene_type=='protein_coding') %>% 
  filter(use_for_ordering=="TRUE") %>% 
  
  arrange(qval)

# my_genes <- genes$gene_short_name %>% head(n = 1000)
# plot_genes_in_pseudotime(cds_subset = cds[my_genes,], min_expr = 0.1, color_by="Minor_cell")

n_top_gene = c(1000,2000,3000)
num_clusters = c(3,4,5,6,7)

n_top_gene <- c(4000,5000,6000)
num_clusters <- 3
for(i in n_top_gene){
  my_genes <- genes$gene_short_name %>% head(n = i)

  lapply(num_clusters, function(x){
    
    p <- plot_pseudotime_heatmap(cds[my_genes,], 
                                 cores = 16, show_rownames = F, 
                                 num_clusters = x,
                                 return_heatmap = T)
    ggsave(plot = p, 
           filename = paste0( "monocle/APSE-Om and Adiopocyte/plots/time-ht-top.",i,"genes_c",x,".tiff" ), 
           # "monocle/APSE-Om and Adiopocyte/plots/time-ht-c6.tiff",
           dpi = 300, width = 12, height = 20, units = 'cm')
    
    
    clusters <- cutree(p$tree_row, k = x)
    df_cluster <- data.frame(clusters)
    # df_cluster[,1] <- as.character(df_cluster[,1])
    colnames(df_cluster) <- "Gene_Clusters"
    df_cluster <- df_cluster %>% 
      rownames_to_column("gene_name") %>% 
      mutate(nCluster = paste0("nCluster-",x)) %>% 
      mutate(nTopGenes = i) %>% 
      as.data.frame()
    
    openxlsx::write.xlsx(df_cluster, overwrite = T,
                         file = paste0( "monocle/APSE-Om and Adiopocyte/plots/time-ht-top.",i,"genes_clusters-",x,".xlsx" ))
    
    
  }) 
  
  
}
rm(i, n_top_gene,num_clusters)

gc()


#


### clusterprofiler

library(clusterProfiler)
library(org.Hs.eg.db)

if(F){
df_gene <- readxl::read_excel("monocle/APSE-Om and Adiopocyte/plots/time-ht-top.5000genes_clusters-3.xlsx")

# kegg and go
if(F){
enrich.go.list <- lapply(1:3, function(x) {
  gene <-  df_gene$gene_name[df_gene$Gene_Clusters == x]
  ego <-   clusterProfiler::enrichGO(
    gene = gene,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05
  ) %>% as.data.frame() %>% mutate(geneCluster=x, Enrich="enrichGO")
}) %>% data.table::rbindlist()


enrich.kegg.list <- lapply(1:3, function(x) {
  gene <-  df_gene$gene_name[df_gene$Gene_Clusters == x]
  entrez_genes <- mapIds(org.Hs.eg.db, gene,'ENTREZID', 'SYMBOL')
  ekegg <- clusterProfiler::enrichKEGG(
    gene = entrez_genes,
    organism = 'hsa', 
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05
  ) %>% as.data.frame() %>% mutate(geneCluster=x, Enrich="enrichGO")
}) %>% data.table::rbindlist()

openxlsx::write.xlsx(list(enrich.go.list, enrich.kegg.list), 
                    file = "monocle/APSE-Om and Adiopocyte/time-ht-top.5000genes_clusters-3-enrichment.xlsx")
}

enrich.go.list <- readxl::read_excel("monocle/APSE-Om and Adiopocyte/time-ht-top.6000genes_clusters-3-enrichment.xlsx", sheet = "Sheet 1")
enrich.kegg.list <- readxl::read_excel("monocle/APSE-Om and Adiopocyte/time-ht-top.6000genes_clusters-3-enrichment.xlsx", sheet = "Sheet 2")

ggdf <- enrich.go.list %>%
# ggdf <- enrich.kegg.list %>%
  dplyr::group_by(geneCluster) %>% 
  arrange(desc(Count), .by_group = T) %>% 
  slice_head(n = 10) %>% 
  mutate(Description=factor(Description, levels=Description),
         ID=factor(ID, levels=ID)
         )

ggdf$Description %>% duplicated %>% table
ggdf$geneCluster %>% table
ggdf$Count %>% range()
xmax <- 205
if(T){
p.1 <- ggplot(ggdf %>% filter(geneCluster == 1), aes(Count, ID))+
  theme_classic()+
  scale_x_continuous(expand = c(0,0), limits = c(0,xmax))+
  geom_bar(stat = 'identity', width = 0.8, fill=alpha("darkred", 0.5))+
  ylab("")+
  geom_text(aes(label=Description, y=ID),x=2, hjust=0)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
  # scale_fill_gradient(low = alpha("pink", 0.5), high = alpha("darkred", 0.5))
  

p.2 <- ggplot(ggdf %>% filter(geneCluster == 2), aes(Count, ID))+
  theme_classic()+
  scale_x_continuous(expand = c(0,0), limits = c(0,xmax))+
  geom_bar(stat = 'identity', width = 0.8, fill=alpha("darkred", 0.5))+
  ylab("")+
  geom_text(aes(label=Description, y=ID),x=2, hjust=0)
  # scale_fill_gradient(low = alpha("pink", 0.5), high = alpha("darkred", 0.5))


p.3 <- ggplot(ggdf %>% filter(geneCluster == 3), aes(Count, ID))+
  theme_classic()+
  scale_x_continuous(expand = c(0,0), limits = c(0,xmax))+
  geom_bar(stat = 'identity', width = 0.8, fill=alpha("darkred", 0.5))+
  ylab("")+
  geom_text(aes(label=Description, y=ID),x=2, hjust=0)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
  # scale_fill_gradient(low = alpha("pink", 0.5), high = alpha("darkred", 0.5))
}

p <- p.1+p.3+p.2+plot_layout(ncol = 1, guides = 'collect', heights = c(10,9,10)) 
rm(xmax, p.1,p.3,p.2)
p 



markers = c("ERBB4", "CDH6","C1QA","C1R", "CD36","FATP4","FATP1")
markers <- diff_time %>%
  filter(gene_type == 'protein_coding') %>% 
  # filter(use_for_ordering=="TRUE") %>% 
  arrange(qval) %>% 
  slice_head(n=12) 
plot_cell_trajectory(cds, markers = markers$gene_name,
                     use_color_gradient = TRUE,
                     show_backbone = T, cell_size = 0.75)

ggsave(plot = last_plot(), 
       filename = "monocle/APSE-Om and Adiopocyte/plots/cds_ordering-Pseudotime-Markers.tiff", 
       width = 25, height=22, dpi = 300,
       units = 'cm')


markers <- diff_tissue %>%
  filter(gene_type == 'protein_coding') %>% 
  # filter(use_for_ordering=="TRUE") %>% 
  arrange(qval) %>% 
  slice_head(n=20) 
plot_cell_trajectory(cds, markers = markers$gene_name,
                     use_color_gradient = TRUE,
                     show_backbone = T, cell_size = 0.75)

ggsave(plot = last_plot(), 
       filename = "monocle/APSE-Om and Adiopocyte/plots/cds_ordering-Tissue-Markers.tiff", 
       width = 30, height=24, dpi = 300,
       units = 'cm')

}









### with cluster infor



cds <- qs::qread("monocle/ad_sub_ad/cds_ordering.qs")


p <- plot_cell_trajectory(cds, color_by = "Pseudotime", show_backbone = T, cell_size = 0.25)

























sc <- qs::qread("./05-ad_merge_Major-Minor-filter.h5ad/05-ad_merge_Major-Minor-filter.qs")
load("h://Rdata/19_geneIDtransform/gene.gtf.gencode.V37.MCPyV.rdata")

meta.data <- data.table::fread("../pydata/04-ad_sub_adi-meta_data.csv")

# cds <- cds[,cds$Cell_id %in% meta.data$Cell_id]
# 
# meta.data <- meta.data %>% 
#   filter(Cell_id %in% cds$Cell_id) %>% 
#   mutate(Cell_id = factor(Cell_id, levels = cds$Cell_id)) %>% 
#   dplyr::arrange(Cell_id)
# 
# cds$cluster = meta.data$cluster

if(T){
  
  sc.sub <- sc[rownames(sc) %in% gene.gtf$gene_name[gene.gtf$gene_type =='protein_coding'],
               sc$Cell_id %in% meta.data$Cell_id]
  
  all(sc.sub$Cell_id == meta.data$Cell_id)
  sc.sub$cluster = meta.data$cluster
  
  
  
  expCnt <- sc.sub@assays$RNA@counts
  expressed_genes <- apply(expCnt, 1, function(x) sum(x!=0)) 
  expCnt <- expCnt[expressed_genes>=50,]
  cellMeta <- AnnotatedDataFrame(sc.sub@meta.data)
  geneMeta <- data.frame(gene_short_name = rownames(expCnt), row.names = rownames(expCnt)) 
  geneMeta <- AnnotatedDataFrame(data = geneMeta)
  
  cds <- newCellDataSet(cellData = expCnt, 
                        phenoData = cellMeta, 
                        featureData = geneMeta,
                        expressionFamily=negbinomial.size())
  rm(expCnt, geneMeta, cellMeta)
  gc()
  
  
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  
  diff_test_res <- differentialGeneTest(cds, cores = 10,
                                        fullModelFormulaStr = "~ Tissue_Site")
  openxlsx::write.xlsx(diff_test_res, file = "monocle/ad_sub_ad/diff_test_res-ad_sub_ad-Tissue_Site.xlsx")
  ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
  
  cds_ordering <- setOrderingFilter(cds = cds, ordering_genes = ordering_genes)
  cds_ordering <- reduceDimension(cds_ordering, max_components = 2,
                                  method = 'DDRTree')
  qs::qsave(cds_ordering, 'monocle/ad_sub_ad/cds_ordering.qs')
  
  cds_ordering <- orderCells(cds_ordering)
  
  cds_ordering$cluster <- as.character(cds_ordering$cluster )
  plot_cell_trajectory(cds_ordering, color_by = "Tissue_Site")+
    plot_cell_trajectory(cds_ordering, color_by = "cluster")+
  plot_cell_trajectory(cds_ordering, color_by = "State", show_state_number = F)+
  plot_cell_trajectory(cds_ordering, color_by = "Pseudotime")
  
  
  
  diff_test_Pseudotime <- differentialGeneTest(cds_ordering, cores = 10,
                                               fullModelFormulaStr = "~sm.ns(Pseudotime)")
  openxlsx::write.xlsx(diff_test_Pseudotime, file = "monocle/ad_sub_ad/diff_test_Pseudotime.xlsx")
  
  my_genes <- row.names(subset(fData(cds_ordering), gene_short_name %in% c("FABP4","MXRA8","ADIPOQ"))) 
  plot_genes_in_pseudotime(cds_subset = cds_ordering[my_genes,], min_expr = 0.1, color_by="Minor_cell")
  
}





































##### old codes




if(F){
sc.sub <- sc %>% subset(Minor_cell  == "APSC-Om" )


expCnt <- sc.sub@assays$RNA@counts
cellMeta <- AnnotatedDataFrame(sc.sub@meta.data)
geneMeta <- AnnotatedDataFrame(data = data.frame(gene_short_name = rownames(expCnt), row.names = rownames(expCnt)))

cds <- newCellDataSet(cellData = expCnt, phenoData = cellMeta, featureData = geneMeta,expressionFamily=negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

diff_test_res <- differentialGeneTest(cds, cores = 10,
                                      fullModelFormulaStr = "~ Tissue_Site")
openxlsx::write.xlsx(diff_test_res, file = "monocle/diff_test_res-APSC-Om-Tissue_Site.xlsx")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))


cds_ordering <- setOrderingFilter(cds = cds, ordering_genes = ordering_genes)
plot_ordering_genes(cds_ordering)

cds_ordering <- reduceDimension(cds_ordering, max_components = 2,
                                method = 'DDRTree')

cds_ordering <- orderCells(cds_ordering)

plot_cell_trajectory(cds_ordering, color_by = "Tissue_Site")

plot_cell_trajectory(cds_ordering, color_by = "State")

qs::qsave(cds_ordering, 'monocle/cds_ordering.qs')

blast_genes <- c("FABP4", 'PDGFA')
plot_genes_jitter(cds_ordering[blast_genes,],
                  grouping = "Tissue_Site",
                  min_expr = 0.1)



subset(cds_ordering, Tissue_Site != "Ovary")


cds_ordering$Tissue_Site


cds_sub <- cds_ordering[,cds_ordering$Tissue_Site != "Ovary"]
cds_sub <- setOrderingFilter(cds = cds_sub, ordering_genes = ordering_genes)
cds_sub <- reduceDimension(cds_sub, max_components = 2,
                           method = 'DDRTree')
cds_sub <- orderCells(cds_sub, reverse = T)
plot_cell_trajectory(cds_sub, color_by = "Tissue_Site")+
  plot_cell_trajectory(cds_sub, color_by = "State")+
  plot_cell_trajectory(cds_sub, color_by = "Pseudotime")


diff_test_Pseudotime <- differentialGeneTest(cds_sub, cores = 10,
                                             fullModelFormulaStr = "~sm.ns(Pseudotime)")
openxlsx::write.xlsx(diff_test_Pseudotime, file = "monocle/diff_test_Pseudotime.xlsx")
qs::qsave(cds_sub, 'monocle/cds_sub.qs')

my_genes <- row.names(subset(fData(cds_sub), gene_short_name %in% c("FABP4","MXRA8","ADIPOQ"))) 
plot_genes_in_pseudotime(cds_subset = cds_sub[my_genes,], min_expr = 0.1, color_by="Tissue_Site")






cds_sub$


load("h://Rdata/19_geneIDtransform/gene.gtf.gencode.V37.MCPyV.rdata")
df_res <- diff_test_res %>% 
  mutate(gene_name = gene_short_name) %>% 
  left_join(gene.gtf) %>% 
  filter(gene_type == "protein_coding")


gg <- df_res$gene_short_name[df_res$qval <= 1e-6] %>% 
  head(n = 100)

}







## monocle3
if(F){
  
  
  expCnt <- sc.sub@assays$RNA@counts
  cellMeta <- sc.sub@meta.data
  geneMeta <- data.frame(gene_short_name = rownames(expCnt), row.names = rownames(expCnt))
  cds <- new_cell_data_set(expression_data = expCnt,
                           cell_metadata = cellMeta,
                           gene_metadata = geneMeta )
  rm(expCnt,cellMeta,geneMeta)
  
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- align_cds(cds, alignment_group = "orig.ident") 
  # residual_model_formula_str = "~ bg.300.loading + 
  # bg.400.loading + bg.500.1.loading + bg.500.2.loading + 
  # bg.r17.loading + bg.b01.loading + bg.b02.loading")
  cds <- reduce_dimension(cds)
  plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "Cell_type_minor")
  
  cds <- cluster_cells(cds, cluster_method = "louvain")
  cds <- learn_graph(cds)
  plot_cells(cds,
             color_cells_by = "Cell_type_minor",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=T)
  
  
  
  
  diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                        fullModelFormulaStr = "~Media")
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  
  
  
  
  
  BiocManager::install("monocle", update = F, force = T)
  library(monocle)
  
  
  
  
}



## monocle2 
sc.sub <- subset(x = sc.reUMAP, subset = Cell_type_minor.1 == "ASPC")
DimPlot(sc.sub, group.by = c('Cell_type_minor', 'Cell_type_minor.1'),raster = T, cols = colorseq)
if(F){
  expCnt <- sc.sub@assays$RNA@counts
  cellMeta <- AnnotatedDataFrame(sc.sub@meta.data)
  geneMeta <- AnnotatedDataFrame(data = data.frame(gene_short_name = rownames(expCnt), row.names = rownames(expCnt)))
  
  cds <- newCellDataSet(cellData = expCnt, phenoData = cellMeta, featureData = geneMeta,expressionFamily=negbinomial.size())
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  diff_test_res <- differentialGeneTest(cds,
                                        fullModelFormulaStr = "~ Cell_type_minor")
  openxlsx::write.xlsx(diff_test_res, file = "monocle/diff_test_res.xlsx")
  ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
  
  
  cds_ordering <- setOrderingFilter(cds = cds, ordering_genes = ordering_genes)
  plot_ordering_genes(cds_ordering)
  
  cds_ordering <- reduceDimension(cds_ordering, max_components = 2,
                                  method = 'DDRTree')
  
  cds_ordering <- orderCells(cds_ordering)
  
  plot_cell_trajectory(cds_ordering, color_by = "Cell_type_minor")
  
  plot_cell_trajectory(cds_ordering, color_by = "State")
  
  qs::qsave(cds_ordering, 'monocle/cds_ordering.qs')
  
  blast_genes <- c("FABP4", 'PDGFA')
  plot_genes_jitter(cds_ordering[blast_genes,],
                    grouping = "State",
                    min_expr = 0.1)
  
  
}



##     only omentum+adipose tissue  monocle2
sc.sub <- subset(x = sc.reUMAP, subset = Cell_type_minor.1 == "ASPC" & Sample_Type %in% c("Metastasis", "Adipose"))
# sc.reUMAP$Sample_Type %>% unique()

DimPlot(sc.sub, group.by = c('Cell_type_minor', 'Cell_type_minor.1'),raster = T, cols = colorseq)
if(F){
  expCnt <- sc.sub@assays$RNA@counts
  cellMeta <- AnnotatedDataFrame(sc.sub@meta.data)
  geneMeta <- AnnotatedDataFrame(data = data.frame(gene_short_name = rownames(expCnt), row.names = rownames(expCnt)))
  cds <- newCellDataSet(cellData = expCnt, phenoData = cellMeta, featureData = geneMeta,expressionFamily=negbinomial.size())
  expCnt <- as.matrix(expCnt)
  
  
  kept.gene <- apply(expCnt, MARGIN = 1, function(x){ sum(x>0) > 10})
  kept.gene <- rownames(expCnt)[kept.gene]
  
  
  rm(expCnt, cellMeta, geneMeta)
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  
  
  diff_test_res <- differentialGeneTest(cds[kept.gene,],
                                        fullModelFormulaStr = "~ Cell_type_minor")
  openxlsx::write.xlsx(diff_test_res, file = "monocle/diff_test_res.aspc-meta-adipose.xlsx")
  ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
  
  
  cds_ordering <- setOrderingFilter(cds = cds, ordering_genes = ordering_genes)
  plot_ordering_genes(cds_ordering)
  
  cds_ordering <- reduceDimension(cds_ordering, max_components = 2,
                                  method = 'DDRTree')
  
  cds_ordering <- orderCells(cds_ordering)
  
  plot_cell_trajectory(cds_ordering, color_by = "Cell_type_minor")
  
  plot_cell_trajectory(cds_ordering, color_by = "State")
  
  qs::qsave(cds_ordering, 'monocle/cds_ordering.aspc-meta-adipose.qs')
  
  blast_genes <- c("FABP4", 'PDGFA')
  plot_genes_jitter(cds_ordering[blast_genes,],
                    grouping = "State",
                    min_expr = 0.1)
  
  
}


##     only omentum tissue  monocle2
sc.sub <- subset(x = sc.reUMAP, subset = Cell_type_minor.1 == "ASPC" & Sample_Type %in% c("Metastasis"))
# sc.reUMAP$Sample_Type %>% unique()

DimPlot(sc.sub, group.by = c('Cell_type_minor', 'Cell_type_minor.1'),raster = T, cols = colorseq)
if(F){
  expCnt <- sc.sub@assays$RNA@counts
  cellMeta <- AnnotatedDataFrame(sc.sub@meta.data)
  geneMeta <- AnnotatedDataFrame(data = data.frame(gene_short_name = rownames(expCnt), row.names = rownames(expCnt)))
  cds <- newCellDataSet(cellData = expCnt, phenoData = cellMeta, featureData = geneMeta,expressionFamily=negbinomial.size())
  expCnt <- as.matrix(expCnt)
  
  
  kept.gene <- apply(expCnt, MARGIN = 1, function(x){ sum(x>0) > 10})
  kept.gene <- rownames(expCnt)[kept.gene]
  
  
  rm(expCnt, cellMeta, geneMeta)
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  
  
  diff_test_res <- differentialGeneTest(cds[kept.gene,],
                                        fullModelFormulaStr = "~ Cell_type_minor")
  openxlsx::write.xlsx(diff_test_res, file = "monocle/diff_test_res.aspc-Metastasis.xlsx")
  ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
  
  
  cds_ordering <- setOrderingFilter(cds = cds, ordering_genes = ordering_genes)
  plot_ordering_genes(cds_ordering)
  
  cds_ordering <- reduceDimension(cds_ordering, max_components = 2,
                                  method = 'DDRTree')
  
  cds_ordering <- orderCells(cds_ordering)
  
  plot_cell_trajectory(cds_ordering, color_by = "Cell_type_minor")
  
  plot_cell_trajectory(cds_ordering, color_by = "State")
  
  qs::qsave(cds_ordering, 'monocle/cds_ordering.aspc-Metastasis.qs')
  
  blast_genes <- c("FABP4", 'PDGFRA')
  plot_genes_jitter(cds_ordering[blast_genes,],
                    grouping = "State",
                    min_expr = 0.1)
  
  diff_test_res.time <- differentialGeneTest(cds_ordering[kept.gene,],
                                             fullModelFormulaStr = "~sm.ns(Pseudotime)")
  
  
  
  sig_gene_names <- diff_test_res.time %>% 
    filter(qval <= 0.01) %>% 
    arrange(qval) %>% 
    slice_head(n=100)
  plot_pseudotime_heatmap(cds_ordering[sig_gene_names$gene_short_name,],
                          num_clusters = 4,
                          cores = 12,
                          show_rownames = T)
  
  
  
  plot_genes_in_pseudotime(cds_ordering[blast_genes,], color_by = "Cell_type_minor")
  
  
  
  
  
  
}


