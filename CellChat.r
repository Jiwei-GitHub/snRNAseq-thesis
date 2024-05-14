
options(string.as.factor = F)
# library(Seurat)
library(tidyverse)
library(patchwork)

library(future)
plan("multicore", workers = 4, future.seed=TRUE)
options(future.globals.maxSize= 20*1024*1024*1024, future.seed = TRUE)


# library(ggsci)

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



setwd("H:/Work_Home/Jiwei_ICMC_data/OvarianCancer/AdiposeTissue/plotting.in.r/CellChat")



if(F){

df.adipocyte <- data.table::fread("../../pydata/04-ad_sub_adi-meta_data.csv") %>% 
  filter(Tissue_Site != "Ovary") %>% 
  mutate(cluster = `leiden_0.4`) %>% 
  select(Cell_id, orig.ident, Pat_ID, Tissue_Site, Minor_cell_2, cluster) %>% 
  mutate(cell_type = paste0(Minor_cell_2,"-C",cluster),
         Pat_Site = paste0(Pat_ID, "-",Tissue_Site))


df.Fibro <- data.table::fread("../../pydata/04-ad_sub_Fibro-meta_data.csv") %>% 
  filter(Tissue_Site != "Ovary") %>% 
  select(Cell_id, orig.ident, Pat_ID, Tissue_Site, Minor_cell_2, cluster) %>% 
  mutate(cell_type = ifelse(Minor_cell_2 == "Meso", paste0(Minor_cell_2,"-C",cluster), Minor_cell_2),
         Pat_Site = paste0(Pat_ID, "-",Tissue_Site)) 


df.macro <- data.table::fread("../../final_plots_in_python/pydata/adfi_immune_reUAMP.csv") %>% 
  filter(Tissue_Site != "Ovary") %>% 
  select(Cell_id, orig.ident, Pat_ID, Tissue_Site, Minor_cell_2, cluster) %>% 
  mutate(cell_type = ifelse(Minor_cell_2 == "Macro", paste0(Minor_cell_2,"-C",cluster), Minor_cell_2),
         Pat_Site = paste0(Pat_ID, "-",Tissue_Site))


df.EpiTumor <- data.table::fread("../../pydata/04-ad-Epi_Tumor-Minor_cell-filter.new.csv") %>% 
  filter(Tissue_Site != "Ovary") %>% 
  mutate(cluster = `leiden_0.2`) %>% 
  select(Cell_id, orig.ident, Pat_ID, Tissue_Site, Minor_cell_2, cluster) %>% 
  mutate(cell_type =  Minor_cell_2,
         Pat_Site = paste0(Pat_ID, "-",Tissue_Site))

openxlsx::write.xlsx(data.frame(cell_type=df_meta.detail.cell_type$cell_type %>% unique %>% sort),file = 'x.xlsx')
df_meta.detail.cell_type <- list(df.adipocyte, df.Fibro, df.macro, df.EpiTumor) %>% 
  data.table::rbindlist() %>% 
  left_join(readxl::read_excel('x.xlsx'))
  





data.table::fwrite(df_meta.detail.cell_type, file = "df_meta.detail.cell_type_major-CellChat.csv")
}

df_meta.detail.cell_type <- data.table::fread("df_meta.detail.cell_type-CellChat.csv")

sc <- qs::qread("../05-ad_merge_Major-Minor-filter.h5ad/05-ad_merge_Major-Minor-filter.qs")
sc <- sc[, sc$Cell_id %in% df_meta.detail.cell_type$Cell_id]

meta.data <- sc@meta.data %>% 
  rownames_to_column("x") %>% 
  select(x, Cell_id) %>% 
  left_join(df_meta.detail.cell_type) %>% 
  as.data.frame() %>% 
  column_to_rownames("x")
sc@meta.data <- meta.data
  

Idents(sc) <- sc$cell_type
set.seed(1234)
sc.ch <- subset(sc, downsample = 1000)
sc.ch <- subset(sc, cell_type != "Mast")


sites <- c("Adj", "Dis")
library(CellChat)
cellchat.list <- list()
for(i in c(1,2)){
  ss <- sites[i]
  Idents(sc) <- sc$cell_type
  sc.down <- subset(sc, Tissue_Site == ss)
  set.seed(1234)
  sc.down <- subset(sc.down, downsample = 2000)
  sc.down <- subset(sc.down, cell_type != "Mast")
  
  
  cellchat <- createCellChat(object = sc.down@assays$RNA@counts, meta = sc.down@meta.data, 
                             group.by = 'cell_type')
  
  cellchat <- setIdent(cellchat, ident.use = "cell_type")
  
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(object = cellchat)
  future::plan("multisession", workers = 10)
  
  # cellchat@idents
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat.tmp <- aggregateNet(cellchat)
  
  if(i==1){cellchat.list <- list(cellchat.tmp)}else{cellchat.list <- append(cellchat.list, cellchat.tmp)}
}

names(cellchat.list) <- sites
qs::qsave(cellchat.list, file = './cellchat.list.by.Tissue_Site.rdata')
  
  
  
  
cellchat.list <- qs::qread('./cellchat-all.list.by.Tissue_Site.qs')

cellchat.list <- qs::qread('./cellchat-all.list.by.Tissue_Site-major.qs')
names(cellchat.list)
cellchat.list <- list(Dis = cellchat.list[['Dis']], 
                      Adj = cellchat.list[['Adj']])
cellchat.list <- lapply(cellchat.list, function(x){
  x <- netAnalysis_computeCentrality(x)
  x
})


# cellchat.list <- qs::qread("cellchat-all.list.by.Tissue_Site.qs")
cellchat.list <- list(Dis = cellchat.list[['Dis']], 
                      Adj = cellchat.list[['Adj']])
ch.dis <- cellchat.list$Dis
ch.adj <- cellchat.list$Adj

a <- ch.dis@idents %>% unique %>% sort
b <- ch.adj@idents %>% unique %>% sort

length(a)
length(b)

group.new = levels(ch.adj@idents)
ch.dis <- liftCellChat(ch.dis, group.new)

cellchat.list <- list(Dis = ch.dis, Adj = ch.adj)
cellchat <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))



gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
ggsave(plot = gg1 + gg2,
  filename = "plots/compareInteractions-number-strength.tiff", 
  height = 12, width = 12,
  units = 'cm', dpi = 600)

dev.off()
tiff(filename = "plots/netVisual_diffInteraction-number.tiff",
     height = 16, width = 16,
     res = 600,units = 'cm',)
# par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)


dev.off()
tiff(filename = "plots/netVisual_diffInteraction-strength.tiff",
     height = 16, width = 16,
     res = 600,units = 'cm',)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

dev.off()
tiff(filename = "plots/netVisual_heatmap.tiff", 
     width = 22, height = 11,
     units = "cm",
     res = 600)
gg1 + gg2
dev.off()




weight.max <- getMaxWeight(cellchat.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cellchat.list)) {
  netVisual_circle(cellchat.list[[i]]@net$count, weight.scale = T, label.edge= F,
                   edge.weight.max = weight.max[2], 
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cellchat.list)[i]))
}











num.link <- sapply(cellchat.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(cellchat.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cellchat.list[[i]], title = names(cellchat.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

dev.off()
tiff(filename = "plots/Compare the major sources and targets.tiff", 
     width = 30, height = 14,
     units = "cm",
     res = 600)
patchwork::wrap_plots(plots = gg)
dev.off()




gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Epi_Tumor", signaling.exclude = NULL)
gg1
gg1 <- gg1+xlim(c(-2e-5, 1e-5))+ylim(c(-1e-5, 4e-6))
gg1

gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Adipocyte", signaling.exclude = NULL)
gg2
gg2 <- gg2+xlim(c(-2e-5, 1e-5))+ylim(c(-1e-5, 2.5e-6))
gg2



dev.off()
tiff(filename = "plots/signalingChanges-Adipocyte-Epitumor.tiff", 
     width = 42, height = 21,
     units = "cm",
     res = 600)
patchwork::wrap_plots(plots = list(gg1, gg2), guides = 'collect')
dev.off()








pathways.show <- c("NCAM", "SPP1", "NRXN", "CADM", "COLLAGEN", "LAMININ",'IGF','FGF','APP',"NEGR",'COMPLEMENT')
LR.show.adj <- extractEnrichedLR(ch.adj, signaling = pathways.show, geneLR.return = FALSE)[,1] %>% as.character()
LR.show.dis <- extractEnrichedLR(ch.dis, signaling = pathways.show, geneLR.return = FALSE)[,1] %>% as.character()
LR.show <- intersect(LR.show.adj, LR.show.dis)
# LR.show <- LR.show[1,1]
LR.show
dev.off()

for(i in 1:length(LR.show)){
  print(i)
  LR <- LR.show[i]
  print(LR)
  
 
  
  tiff(filename = paste0("plots/LR.show/",LR,"-Dis.tiff"), 
       width = 16, height = 16,
       res = 300, units = 'cm')
  
  # par(mfrow = c(2,2), xpd=TRUE)
  nname <- paste0(LR,"-Dis")
  
  netVisual_individual(ch.dis, signaling = pathways.show, 
                       # sources.use = c("Adipocyte", "Epi_Tumor"),
                       # targets.use = c("Adipocyte", "Epi_Tumor"),
                       pairLR.use = LR,
                       signaling.name = nname,
                       layout = "circle")
  title(main = nname,outer = F)

  print(nname)  
  dev.off()
  
  
  tiff(filename = paste0("plots/LR.show/",LR,"-Adj.tiff"), 
       width = 16, height = 16,
       res = 300, units = 'cm')
  nname <- paste0(LR,"-Adj")
  netVisual_individual(ch.adj, signaling = pathways.show, 
                       # sources.use = c("Adipocyte", "Epi_Tumor"),
                       # targets.use = c("Adipocyte", "Epi_Tumor"),
                       pairLR.use = LR, 
                       signaling.name = nname,
                       layout = "circle")
  title(main = nname,outer = F)
  
  print(nname) 
  dev.off()
  
  
}





pathways.shows <- c("NCAM", "SPP1", "NRXN", "CADM", "COLLAGEN", "LAMININ",'IGF','FGF','APP',"NEGR",'COMPLEMENT')

pathways.shows <- intersect(ch.dis@netP$pathways, ch.adj@netP$pathways)
dev.off()

for(j in 1:length(pathways.shows)){
pathways.show <- pathways.shows[j]
print(paste0(j, "-",pathways.show))
weight.max <- getMaxWeight(cellchat.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
# vertex.receiver = seq(1,10) # Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells


tiff(filename = paste0("plots/pathways.show/",pathways.show,".tiff"), 
      width = 35, height = 16,
      res = 300, units = 'cm')
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cellchat.list)) {
  netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, 
                      vertex.receiver = vertex.receiver, 
                      # edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste0(pathways.show, "-",names(cellchat.list)[i]))
}



dev.off()


}












# 
# 
# gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = c("Adipocyte-C2"), signaling.exclude = NULL)+
#   xlim(c(-.51e-5, 0.5e-5))+ylim(c(-1.5e-5, 1e-5))
# gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = c("Adipocyte-C3"), signaling.exclude = NULL)+
#   xlim(c(-1e-5, 2e-5))+ylim(c(-1e-5, 2e-5))
# gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = c("Adipocyte-C6"), signaling.exclude = NULL)+
#   xlim(c(-1e-5, 0.5e-5))+ylim(c(-1e-5, 1.2e-5))
# 
# 
# 
# patchwork::wrap_plots(plots = list(gg2,gg3, gg5), guides = 'collect')
# 
# 
# dev.off()
# tiff(filename = "plots/signalingChanges-Adipocyte-3clsts.vs.Firboblast.tiff", 
#      width = 42, height = 12,
#      units = "cm",
#      res = 600)
# patchwork::wrap_plots(plots = list(gg2,gg3, gg5), guides = 'collect')
# dev.off()







library(ComplexHeatmap)
i = 1
pathway.union <- union(cellchat.list[[i]]@netP$pathways, cellchat.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cellchat.list[[i]], pattern = "outgoing", 
                                        signaling = pathway.union,
                                        title = names(cellchat.list)[i], 
                                        width = 12, height = 30)
ht2 = netAnalysis_signalingRole_heatmap(cellchat.list[[i+1]], pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(cellchat.list)[i+1],
                                        width = 12, height = 30)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


dev.off()
tiff(filename = "plots/signalingRole_heatmap-outgoing.tiff", 
     width = 30, height = 40,
     units = "cm",
     res = 600)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


ht1 = netAnalysis_signalingRole_heatmap(cellchat.list[[i]], pattern = "incoming", 
                                        signaling = pathway.union,
                                        title = names(cellchat.list)[i], 
                                        width = 12, height = 30)
ht2 = netAnalysis_signalingRole_heatmap(cellchat.list[[i+1]], pattern = "incoming", 
                                        signaling = pathway.union, 
                                        title = names(cellchat.list)[i+1],
                                        width = 12, height = 30)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


dev.off()
tiff(filename = "plots/signalingRole_heatmap-incoming.tiff", 
     width = 30, height = 40,
     units = "cm",
     res = 600)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()










ch.dis@netP$pathways %>% sort

ch.adj@netP$pathways %>% sort












# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Adj"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "Tissue_Site", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "Adj",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Dis",ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        sources.use = 4, 
                        sources.use = c('Firboblast', "Meso-C0","Meso-C1"), 
                        targets.use = c("Macro-C1"), 
                        angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(cellchat.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                        sources.use = c('Firboblast', "Meso-C0","Meso-C1"), 
                        targets.use = c("Macro-C1"), 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(cellchat.list)[2]))
#> Comparing communications on a merged object
 gg2





 
 
 dev.off()
 
 # Chord diagram
 pathways.show <- c("NEGR") 
 par(mfrow = c(1,2), xpd=TRUE)
 for (i in 1:length(cellchat.list)) {
   netVisual_aggregate(cellchat.list[[i]], 
                       signaling = pathways.show, 
                       layout = "chord", 
                       signaling.name = paste(pathways.show, names(cellchat.list)[i]))
 }
 
 
 
 
 
 cellchat.list$Dis@netP$pathways
 
 
 
 
 
 
 
 
 
 
 
 
 
 







##
ids <- names(cellchat.list);ids
names(ids) <- ids
# id <- ids[4]
p.list <- lapply(ids, function(id){
  ch <- cellchat.list[[id]]
  groupSize <- as.numeric(table(ch@idents))
  df.signal <- ch@LR %>% as.data.frame()
  
  netVisual_circle(ch@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = id)
  
})
  
par(mfrow = c(1,2))




dev.off()








#'@
#'#'@
#'#'@
#'#'@
#'#'@
#'#'@
#'#'@
#'#'@
#'#'@


cellchat.list <- qs::qread("cellchat-all.list.by.Tissue_Site.qs")
cellchat.list <- list(Dis = cellchat.list[['Dis']], 
                      Adj = cellchat.list[['Adj']])
ch.dis <- cellchat.list$Dis
ch.adj <- cellchat.list$Adj

a <- ch.dis@idents %>% unique %>% sort
b <- ch.adj@idents %>% unique %>% sort

length(a)
length(b)

group.new = levels(ch.adj@idents)
ch.dis <- liftCellChat(ch.dis, group.new)

cellchat.list <- list(Dis = ch.dis, Adj = ch.adj)
cellchat <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))


mat <- ch.dis@net$weight
groupSize <- as.numeric(table(ch.dis@idents))
par(mfrow = c(4,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()
tiff(filename = "plots/single-data-ch.dis.tiff",
     height = 4*8.5, width = 5*8.5,
     units = 'cm',res = 600)






mat <- ch.adj@net$weight
groupSize <- as.numeric(table(ch.adj@idents))
par(mfrow = c(4,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
tiff(filename = "plots/single-data-ch.adj.tiff",
     height = 4*8.5, width = 5*8.5,
     units = 'cm',res = 600)






# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
ch.adj@idents %>% unique()
p1 <- netVisual_bubble(ch.adj, sources.use = c("Adipocyte-C2", "Adipocyte-C3", "Adipocyte-C6"), 
                       targets.use = c("Epi_Tumor", "Epi_Tumor-Cycling"),
                       remove.isolate = FALSE)
p2 <- netVisual_bubble(ch.adj, 
                       targets.use = c("Adipocyte-C2", "Adipocyte-C3", "Adipocyte-C6"), 
                       sources.use = c("Epi_Tumor", "Epi_Tumor-Cycling"),
                       remove.isolate = FALSE)



p1+p2


ch.adj@net




dev.off()
tiff(filename = "plots/single-data-ch.adj-adi-to-epi.tiff",
     height = 20, width = 20,
     units = 'cm',res = 600)
netVisual_chord_gene(ch.adj, sources.use = c("Adipocyte-C2", "Adipocyte-C3", "Adipocyte-C6"), 
                     targets.use = c("Epi_Tumor", "Epi_Tumor-Cycling"), lab.cex = 0.3,
                     slot.name = "netP",
                     legend.pos.x = 15)
dev.off()





ch.adj <- computeNetSimilarity(ch.adj, type = "functional")
ch.adj <- netEmbedding(ch.adj, type = "functional")
ch.adj <- netClustering(ch.adj, type = "functional")
netVisual_embedding(ch.adj, type = "functional", label.size = 3.5)


ch.adj <- computeNetSimilarity(ch.adj, type = "structural")
ch.adj <- netEmbedding(ch.adj, type = "structural")
ch.adj <- netClustering(ch.adj, type = "structural")
netVisual_embedding(ch.adj, type = "structural", label.size = 3.5)


pathways.show <- c("FLRT")
LR.show <- extractEnrichedLR(ch.adj, signaling = pathways.show, geneLR.return = FALSE)
LR.show
LR.show <- LR.show[1,1]
LR.show
netVisual_individual(ch.dis, signaling = pathways.show, 
                     # sources.use = c("Adipocyte-C2", "Adipocyte-C3", "Adipocyte-C6"),
                     # targets.use = c("Epi_Tumor", "Epi_Tumor-Cycling"),
                     
                     # targets.use = c("Adipocyte-C2", "Adipocyte-C3", "Adipocyte-C6"),
                     # sources.use = c("Epi_Tumor", "Epi_Tumor-Cycling"),
                     
                     pairLR.use = LR.show,
                     layout = "circle")


plotGeneExpression(ch.adj, features = "IL2")&ylim(0,2)

ch.adj@netP
