library(data.table)
library(tidyverse)
library(patchwork)
library(Matrix)
library(CellChat) #v2.2.0
library(Seurat)
library(reticulate)
library(dplyr)
library(ggplot2)

packageVersion('CellChat')

use_python('/mnt/data/home/haoxueyu/.conda/envs/scRNA/bin/python')
setwd('/mnt/data/home/haoxueyu/Project/Lixuwen-ccRCC_20250509/scRNAseq/results/P5/Rdata')

sce <- readRDS("ccRCC_epi_immune.rds")#标准化后的数据
sce$samples <- sce$sample_id

dim(sce) #33343 84229

table(is.na(sce@meta.data$celltype))
# FALSE  TRUE 
# 81810  2419 
sce <- sce[, !is.na(sce$celltype)]
table(is.na(sce@meta.data$celltype))
colnames(sce@meta.data)

####先排序再构建cellchat对象，便于后续展示，需要添加samples、celltype列
my_cells <- c("cDC1","cDC2","cDC3", "pDC","Mast cell", "Neutrophil", "Mono_CD14","Mono_CD16",
              "TAM_BNC2","TAM_CCL3","TAM_FN1","TAM_TREM2","TAM_CD163","TAM_FOLR2", "TAM_IL32")
lym_cells <- c("B_BANK1","B_AREG","Plasma_IGLC2+","Plasma_IGLC2-","CD16+NK","CD16+NK_FCER1G","CD56+NK","NKT",
               "CD4+Treg","CD4+T_FOXP1","CD4+T_CCR6","CD4+T_RGCC","CD4+T_FOS","CD4+T_TSHZ2",
               "CD8+Tex_HAVCR2","CD8+Tex_CD74","CD8+Tex_LAG3","CD8+ProT","CD8+T_PARP8","CD8+T_IFNG","CD8+Tact")

epi_state <- c("PT","Stress","MHCII","Cell_cycle","Cell_death","EMT")

sce$celltype <- factor(sce$celltype, levels = c(my_cells, lym_cells, epi_state))

sce$sample_id <- factor(sce$sample_id, levels = c("P6N1", "P6N2", "P6N3","P7N1","P7N2","P7N3","P8N1","P8N2","P8N3","P9N1", "P9N2","P9N3","P11N1", "P11N2", "P11N3"))
unique(sce$sample_id)
# P7N1  P7N2  P7N3  P9N1  P9N2  P9N3  P6N1  P6N2  P6N3  P8N1  P8N2  P8N3  P11N1 P11N2 P11N3

levels(sce$celltype) #42
table(sce$celltype)
table(sce$Category)

create_cellchat <- function(seu) {
  library(CellChat)
  library(future)
  
  cc <- createCellChat(seu, group.by = "celltype")
  cc@DB <- CellChatDB.human
  
  cc <- subsetData(cc)
  
  options(future.seed = TRUE)
  options(future.globals.maxSize = 20 * 1024^3)
  future::plan("multisession", workers = 5)
  
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 20)
  
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  
  future::plan("sequential")
  cc
}


#################################
nodule.list <- SplitObject(sce, split.by = "sample_id")
names(nodule.list)
nodule.list


cellchat.list <- lapply(nodule.list, function(seu){
  cellchat <- createCellChat(object = seu, group.by = "celltype")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 20) #默认10
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  return(cellchat)
})


sample_order <- c("P6N1", "P6N2", "P6N3","P7N1","P7N2","P7N3","P8N1","P8N2","P8N3","P9N1", "P9N2","P9N3","P11N1", "P11N2", "P11N3")
cellchat.merge <- mergeCellChat( cellchat.list[sample_order], add.names = sample_order )
save(cellchat.merge, file='ccRCC.cellchat_order.rdata')

sample_palette=c("#1F77B4FF","#AEC7E8FF","#9EDAE5FF","#9467BDFF","#9C9EDEFF","#C5B0D5FF","#FF7F0EFF","#FFBB78FF","#E7CB94FF","#2CA02CFF","#98DF8AFF","#91D1C2FF","#D6616BFF","#FF9896FF", "#F7B6D2FF") #"#7F7F7FFF","#C7C7C7FF"

prefix='../plot/'
pdf("../plot/compareInteractions_order.pdf", width = 12, height = 5)
p1 <- compareInteractions(cellchat.merge, show.legend = FALSE, group = names(cellchat.list), color.use = sample_palette, title.name = "count",  angle.x = 45, x.lab.rot = TRUE)
p2 <- compareInteractions(cellchat.merge, show.legend = FALSE, group = names(cellchat.list), measure = "weight", color.use = sample_palette, title.name = "weight", angle.x = 45, x.lab.rot = TRUE)
p1|p2
dev.off()


