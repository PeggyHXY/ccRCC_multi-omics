library(data.table)
library(tidyverse)
library(patchwork)
library(Matrix)
library(CellChat) #v2.2.0
library(Seurat)
library(reticulate)
library(dplyr)
library(ggplot2)

use_python('/mnt/data/home/haoxueyu/.conda/envs/scRNA/bin/python')
setwd('/mnt/data/home/haoxueyu/Project/Lixuwen-ccRCC_20250509/scRNAseq/results/P5/Rdata')

sce <- readRDS("ccRCC_epi_immune.rds")#标准化后的数据
sce$samples <- sce$sample_id

dim(sce)

table(is.na(sce@meta.data$celltype))
# FALSE  TRUE 
# 81810  2419 
sce <- sce[, !is.na(sce$celltype)]
table(is.na(sce@meta.data$celltype))

my_cells <- c("cDC1","cDC2","cDC3", "pDC","Mast cell", "Neutrophil", "Mono_CD14","Mono_CD16", 
              "TAM_BNC2","TAM_CCL3","TAM_FN1","TAM_TREM2","TAM_CD163","TAM_FOLR2", "TAM_IL32")
lym_cells <- c("B_BANK1","B_AREG","Plasma_IGLC2+","Plasma_IGLC2-","CD16+NK","CD16+NK_FCER1G","CD56+NK","NKT",
               "CD4+Treg","CD4+T_FOXP1","CD4+T_CCR6","CD4+T_RGCC","CD4+T_FOS","CD4+T_TSHZ2",
               "CD8+Tex_HAVCR2","CD8+Tex_CD74","CD8+Tex_LAG3","CD8+ProT","CD8+T_PARP8","CD8+T_IFNG","CD8+Tact")

epi_state <- c("PT","Stress","MHCII","Cell_cycle","Cell_death","EMT")

sce$celltype <- factor(sce$celltype, levels = c(my_cells, lym_cells, epi_state))

DimPlot(sce,reduction="umap", group.by="celltype", raster=F,label=T,repel=T)+NoLegend()

sce$sample_id <- factor(sce$sample_id, levels = c("P6N1", "P6N2", "P6N3","P7N1","P7N2","P7N3","P8N1","P8N2","P8N3","P9N1", "P9N2","P9N3","P11N1", "P11N2", "P11N3"))

unique(sce$sample_id)

cellchat <- createCellChat(object = sce@assays$RNA@data, meta = sce@meta.data, group.by = "celltype")
dim(cellchat@data)          # 原始数据维度
dim(cellchat@data.signaling) # 处理后维度(应不为0)

cellchat <- setIdent(cellchat, ident.use = "celltype") 
levels(cellchat@idents) 

groupSize <- as.numeric(table(cellchat@idents))
cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat) 
dim(cellchat@data.signaling)#1389 81810

options(future.seed = TRUE)  # 启用并行安全的随机数生成
cellchat <- identifyOverExpressedGenes(cellchat) #耗时
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
future::plan(strategy = "multisession", workers = 10)
options(future.globals.maxSize =  10 * 1024^3)  # 增加到10GB
# cellchat <- computeCommunProb(cellchat, type = "triMean") #耗时
cellchat <- computeCommunProb(cellchat, type = "triMean", population.size =TRUE ) #耗时

cellchat <- filterCommunication(cellchat, min.cells = 20)  #原为10
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
save(cellchat,file='cellchat_epi_immune_v2.rdata') #2.7G

sum(cellchat@net$count) #50076
sum(cellchat@net$weight) #576.0701

df.net <- subsetCommunication(cellchat)
write.csv(df.net,"cellchat_epi_immune.df.net.csv")

# 查看所有互作强度的统计摘要
summary(c(cellchat@net$prob))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000000 0.0000000 0.0000000 0.0001188 0.0000000 0.2138416 

# 使用单线程
future::plan("sequential")

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
future::plan("multisession", workers = 4) 
options(future.globals.maxSize =  5 * 1024^3)
print(getOption("future.globals.maxSize")) #5368709120
cellchat <- netClustering(cellchat, type = "functional", nCores = 8)
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)


netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

#---------------------------------------------------------------------------------------------------------

load('cellchat_epi_immune_v2.rdata')

#为保持一致，仍命名为fig6 看着两次结果一致
pdf("../plot/fig6_signalingRole_scatter.pdf", width = 5, height = 5)
netAnalysis_signalingRole_scatter(cellchat,dot.alpha = 0)
dev.off()

pdf("../plot/fig6_signalingRole_heatmap.pdf", width = 11, height = 12)
p1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 8, width = 10, height = 25)
p2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 8, width = 10, height = 25)
p1 + p2
dev.off()

#Top6
source('../../../codes/netAnalysis_signalingRole_heatmap_ed.R')

choose_pathway=c("MHC-I","MIF","MHC-II","SPP1","APP","ADGRE")
pdf("../plot/fig6_signalingRole_heatmap_top6.pdf", width = 11, height = 12)
p1 <- netAnalysis_signalingRole_heatmap_ed(cellchat, signaling=choose_pathway, pattern = "incoming", font.size = 8, width = 9, height = 5, color.palette = c("#51C3CC","white", "#DF0000"))
p2 <- netAnalysis_signalingRole_heatmap_ed(cellchat, signaling=choose_pathway, pattern = "outgoing", font.size = 8, width = 9, height = 5, color.palette = c("#51C3CC","white", "#DF0000"))
p1 + p2
dev.off()


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE,mar = c(0, 1, 1, 0.5))

pdf("../plot/fig6_netVisual_circle.pdf", width = 8, height = 10)
netVisual_circle(cellchat@net$weight, 
                 sources.use = c("PT","Cell_cycle","Cell_death","MHCII","EMT","Stress"),
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge = F, 
                 title.name = "Outgoing weights/strength") 
netVisual_circle(cellchat@net$weight, 
                 targets.use = c("PT","Cell_cycle","Cell_death","MHCII","EMT","Stress"),
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge = F, 
                 title.name = "Incoming weights/strength")

dev.off()
