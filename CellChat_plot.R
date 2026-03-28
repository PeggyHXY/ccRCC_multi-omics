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
load('cellchat_epi_immune.rdata')

my_cells <- c("cDC1","cDC2","cDC3", "pDC","Mast cell", "Neutrophil", "Mono_CD14","Mono_CD16", 
              "TAM_BNC2","TAM_CCL3","TAM_FN1","TAM_TREM2","TAM_CD163","TAM_FOLR2", "TAM_IL32")
lym_cells <- c("B_BANK1","B_AREG","Plasma_IGLC2+","Plasma_IGLC2-","CD16+NK","CD16+NK_FCER1G","CD56+NK","NKT",
               "CD4+Treg","CD4+T_FOXP1","CD4+T_CCR6","CD4+T_RGCC","CD4+T_FOS","CD4+T_TSHZ2",
               "CD8+Tex_HAVCR2","CD8+Tex_CD74","CD8+Tex_LAG3","CD8+ProT","CD8+T_PARP8","CD8+T_IFNG","CD8+Tact")

epi_state <- c("PT","Stress","MHCII","Cell_cycle","Cell_death","EMT")

# 使用人类的数据库（如果你的数据是小鼠，换成 CellChatDB.mouse）
CellChatDB <- CellChatDB.human



pdf("../../plot/fig6_netVisual_bubble_lym.pdf", width = 15, height = 4)
netVisual_bubble(cellchat, 
                 sources.use = c("PT","Stress","MHCII","Cell_cycle","Cell_death","EMT"), # 所有发送细胞群
                 targets.use = lym_cells, # 所有接收细胞群
                 signaling = checkpoint_interactions$pathway_name, # 指定我们要看的互作
                 remove.isolate = TRUE) + # 移除没有互作的组
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()


pdf("../../plot/fig6_netVisual_bubble_my_epi.pdf", width = 15, height = 5)
netVisual_bubble(cellchat, 
                 sources.use = my_cells, # 所有发送细胞群
                 targets.use = c("PT","Stress","MHCII","Cell_cycle","Cell_death","EMT"), # 所有接收细胞群
                 signaling = checkpoint_interactions$pathway_name, # 指定我们要看的互作
                 remove.isolate = TRUE) + # 移除没有互作的组
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()

pdf("../../plot/fig6_netVisual_bubble_lym_epi.pdf", width = 15, height = 5)
netVisual_bubble(cellchat, 
                 sources.use = lym_cells, # 所有发送细胞群
                 targets.use = c("PT","Stress","MHCII","Cell_cycle","Cell_death","EMT"), # 所有接收细胞群
                 signaling = checkpoint_interactions$pathway_name, # 指定我们要看的互作
                 remove.isolate = TRUE) + # 移除没有互作的组
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()

netVisual_bubble(cellchat, 
                 sources.use = my_cells, # 所有发送细胞群
                 targets.use = lym_cells, # 所有接收细胞群
                 signaling = checkpoint_interactions$pathway_name, # 指定我们要看的互作
                 remove.isolate = TRUE) + # 移除没有互作的组
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()

netVisual_heatmap(cellchat, color.heatmap = "Reds", measure='weight')
                  
checkpoint_pathway <- list("Immune_Checkpoint" = checkpoint_interactions$pathway_name)
cellchat@netP$pathways <- c(checkpoint_pathway, cellchat@netP$pathways) # 添加到pathways中

netVisual_heatmap(cellchat, 
                  signaling = c("CD86","PD-L1","PDL2"), 
                  color.heatmap = "Reds",
                  title.name = "Immune Checkpoint Interaction Strength")

# cellchat@meta[['celltype']] <- factor(cellchat@meta[['celltype']], levels = c(my_cells, lym_cells, epi_state))

cellchat@idents <- as.factor(cellchat@meta[['celltype']])
levels(cellchat@idents)

df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "epi_immune_net_pathway.csv")



pdf("../../plot/fig6_signalingRole_scatter.pdf", width = 5, height = 5)
netAnalysis_signalingRole_scatter(cellchat, dot.alpha = 0)
dev.off()

netAnalysis_signalingRole_scatter(cellchat, signaling = c("MIF", "SPP1"))
pdf("../../plot/fig6_signalingRole_scatter_MHC.pdf", width = 5, height = 5)
netAnalysis_signalingRole_scatter(cellchat, signaling = c("MHC-I"))
netAnalysis_signalingRole_scatter(cellchat, signaling = c("MHC-II"))
dev.off()

netAnalysis_signalingRole_scatter(cellchat, signaling = c("CSF"))


pdf("../../plot/fig6_signalingRole_heatmap.pdf", width = 11, height = 12)
p1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 8, width = 10, height = 25)
p2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 8, width = 10, height = 25)
p1 + p2
dev.off()

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways

# 可视化“整体通讯活跃度”（all），即既发送也接收的综合能力
pdf("../../plot/fig6_signalingRole_heatmap_all.pdf", width = 11, height = 12)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "all", signaling = pathways.show.all,  width = 15, height = 25)
dev.off()

head(pathways.show.all,50)
pdf("../../plot/fig6_signalingRole_heatmap_head30.pdf", width = 11, height = 10)
p1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = head(pathways.show.all,30), pattern = "incoming", font.size = 8, width = 12, height = 15)
p2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = head(pathways.show.all,30),  pattern = "outgoing", font.size = 8, width = 12, height = 15)
p1 + p2
dev.off()


pdf("../../plot/fig6_signalingRole_heatmap_head50.pdf", width = 13, height = 12)
p1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = head(pathways.show.all,50), pattern = "incoming", font.size = 8, width = 12, height = 15)
p2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = head(pathways.show.all,50),  pattern = "outgoing", font.size = 8, width = 12, height = 15)
p1 + p2
dev.off()

pdf("../../plot/fig6_signalingRole_heatmap_MHC.pdf", width = 13, height = 4)
p1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MHC-I","MHC-II"), pattern = "incoming", font.size = 8, width = 12, height = 5)
p2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MHC-I","MHC-II"),  pattern = "outgoing", font.size = 8, width = 12, height = 5)
p1 + p2
dev.off()


head_path <- c("MHC-I","MIF","MHC-II","SPP1","APP","ADGRE")
pdf("../../plot/fig6_signalingRole_network_head.pdf", width = 10, height = 3)
for (pathway.show in head_path) {
  netAnalysis_signalingRole_network(cellchat, signaling = pathway.show ,  width = 15, height = 2.5, font.size = 10) 
}
dev.off()


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE,mar = c(0.5, 1, 3, 0.5))
pdf("../../plot/fig6_netVisual_circle.pdf", width = 9, height = 9)
netVisual_circle(cellchat@net$weight, 
                 # sources.use=c("CD8+Tex_HAVCR2","CD8+Tex_CD74","CD8+Tex_LAG3","CD8+ProT","CD8+T_PARP8","CD8+T_IFNG","CD8+Tact"),
                 sources.use = c("PT","Cell_cycle","Cell_death","MHCII","EMT","Stress"),
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge = F, 
                 title.name = "Epithelial outgoing weights/strength")

netVisual_circle(cellchat@net$weight, 
                 targets.use = c("PT","Cell_cycle","Cell_death","MHCII", "EMT","Stress"),
                 # targets.use = c("CD8+Tex_HAVCR2","CD8+Tex_CD74","CD8+Tex_LAG3","CD8+ProT","CD8+T_PARP8","CD8+T_IFNG","CD8+Tact"),
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge = F, 
                 title.name = "Epithelial incoming weights/strength") 

dev.off()

# 定义一个包含常见免疫检查点的基因列表
immune_checkpoint_genes <- c('CD274', 'PDCD1', 'CTLA4', 'CD80', 'CD86', 'CD28', 
                             'TIGIT', 'CD226', 'CD96', 'PVR', 'PVRL2', 'NECTIN2',
                             'HAVCR2', 'LGALS9', 'LAG3', 'FGL1', 'LSECtin',
                             'SIGLEC10', 'CD24', 'CD47', 'SIRPA',
                             'VTCN1', 'HHLA2', 'TMIGD2')
# 注意：CD274 = PD-L1, PDCD1 = PD-1, HAVCR2 = TIM-3, VTCN1 = B7-H4

# 从CellChatDB中筛选出涉及我们感兴趣基因的互作对
# 这里我们在interaction_name中搜索包含这些基因的条目
checkpoint_interactions <- subsetDB(CellChatDB, 
                                    search = "Immune Checkpoint") 
# 更精确的方法：直接通过基因名匹配
# 找出所有配体或受体在我们的基因列表中的互作对
idx <- which((CellChatDB$interaction$ligand %in% immune_checkpoint_genes) | 
               (CellChatDB$interaction$receptor %in% immune_checkpoint_genes))
checkpoint_interactions <- CellChatDB$interaction[idx, ]

colnames(checkpoint_interactions)
# [1] "interaction_name"         "pathway_name"             "ligand"                   "receptor"                
# [5] "agonist"                  "antagonist"               "co_A_receptor"            "co_I_receptor"           
# [9] "annotation"               "interaction_name_2"       "evidence"                 "is_neurotransmitter"     
# [13] "ligand.symbol"            "ligand.family"            "ligand.location"          "ligand.keyword"          
# [17] "ligand.secreted_type"     "ligand.transmembrane"     "receptor.symbol"          "receptor.family"         
# [21] "receptor.location"        "receptor.keyword"         "receptor.surfaceome_main" "receptor.surfaceome_sub" 
# [25] "receptor.adhesome"        "receptor.secreted_type"   "receptor.transmembrane"   "version" 

# 查看我们找到了哪些免疫检查点互作对
print(checkpoint_interactions[, c('ligand', 'receptor', 'annotation')])

print(checkpoint_interactions[,'pathway_name'])
# [1] "GALECTIN" "GALECTIN" "GALECTIN" "THBS"     "THBS"     "THBS"     "THBS"     "THBS"     "CD80"     "CD80"     "CD80"    
# [12] "CD86"     "CD86"     "CD96"     "CD96"     "ICOS"     "ICOS"     "NECTIN"   "NECTIN"   "NECTIN"   "NECTIN"   "NECTIN"  
# [23] "PD-L1"    "PDL2"     "PVR"      "PVR"      "GALECTIN" "GALECTIN" "SIRP" 

# 绘制所有细胞群间免疫检查点互作的气泡图
pdf("../../plot/fig6_netVisual_bubble_my.pdf", width = 15, height = 4)
netVisual_bubble(cellchat, 
                 sources.use = c("PT","Stress","MHCII","Cell_cycle","Cell_death","EMT"), # 所有发送细胞群
                 targets.use = my_cells, # 所有接收细胞群
                 signaling = checkpoint_interactions$pathway_name, # 指定我们要看的互作
                 remove.isolate = TRUE) + # 移除没有互作的组
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()


pdf("../../plot/fig6_netVisual_circle_single.pdf", width = 6, height = 6)
for (s in c("PT","Cell_cycle","Cell_death","MHCII","EMT","Stress")){
  
  netVisual_circle(cellchat@net$weight, 
                   sources.use = s, 
                   vertex.weight = groupSize,
                   weight.scale = T, 
                   label.edge = F, 
                   title.name = s )
}
dev.off()

pheatmap::pheatmap(cellchat@net$count, border_color = "black", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T,number_color="black",number_format = "%.0f")



### unique receive pathway in States
sig = cellchat@LR$LRsig
centr_path = names(cellchat@netP$centr)
centr_path

df = data.frame(
  pathway = centr_path,
  PT_out = unlist(lapply(cellchat@netP$centr,function(x){x[['outdeg']]['PT']})) ,
  PT_income = unlist(lapply(cellchat@netP$centr,function(x){x[['indeg']]['PT']})) ,
  Cell_death_out = unlist(lapply(cellchat@netP$centr,function(x){x[['outdeg']]['Cell_death']})) ,
  Cell_death_income = unlist(lapply(cellchat@netP$centr,function(x){x[['indeg']]['Cell_death']})) ,
  Cell_cycle_out = unlist(lapply(cellchat@netP$centr,function(x){x[['outdeg']]['Cell_cycle']})) ,
  Cell_cycle_income = unlist(lapply(cellchat@netP$centr,function(x){x[['indeg']]['Cell_cycle']})) ,
  EMT_out = unlist(lapply(cellchat@netP$centr,function(x){x[['outdeg']]['EMT']})) ,
  EMT_income = unlist(lapply(cellchat@netP$centr,function(x){x[['indeg']]['EMT']})) ,
  Stress_out = unlist(lapply(cellchat@netP$centr,function(x){x[['outdeg']]['Stress']})) ,
  Stress_income = unlist(lapply(cellchat@netP$centr,function(x){x[['indeg']]['Stress']})) ,
  MHCII_out = unlist(lapply(cellchat@netP$centr,function(x){x[['outdeg']]['MHCII']})) ,
  MHCII_income = unlist(lapply(cellchat@netP$centr,function(x){x[['indeg']]['MHCII']})) ,
  row.names = centr_path)

head(df)

dim(df)
df = filter(df,rowSums(df[,-1])!=0)
df[df$PT_out!=0 & df$Cell_death_out == 0 & df$Cell_cycle_out==0 & df$EMT_out==0 & df$Stress_out == 0 & df$MHCII_out==0,]$pathway
#
df[df$PT_out==0 & df$Cell_death_out == 0 & df$Cell_cycle_out==0 & df$EMT_out==0 & df$Stress_out != 0 & df$MHCII_out==0,]$pathway
#
df[df$PT_out==0 & df$Cell_death_out != 0 & df$Cell_cycle_out==0 & df$EMT_out==0 & df$Stress_out == 0 & df$MHCII_out==0,]$pathway
# "IL2"
df[df$PT_out==0 & df$Cell_death_out == 0 & df$Cell_cycle_out!=0 & df$EMT_out==0 & df$Stress_out == 0 & df$MHCII_out==0,]$pathway
# "DHT"   "ADGRB"
df[df$PT_out==0 & df$Cell_death_out == 0 & df$Cell_cycle_out==0 & df$EMT_out!=0 & df$Stress_out == 0 & df$MHCII_out==0,]$pathway
# "FN1"       "NEGR"      "SEMA3"     "Adenosine" 
df[df$PT_out==0 & df$Cell_death_out == 0 & df$Cell_cycle_out==0 & df$EMT_out==0 & df$Stress_out == 0 & df$MHCII_out!=0,]$pathway
# "CSF" "BMP"


pdf("../../plot/fig6_netVisual_bubble_outgoing.pdf", width = 16, height = 5)
netVisual_bubble(cellchat,signaling = c("IL2","DHT","ADGRB","FN1","NEGR","SEMA3", "Adenosine","CSF","BMP"), #pathways,
                 remove.isolate = T,
                 sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"),thresh=0.01)
dev.off()

pdf("../../plot/fig6_netVisual_bubble_EMToutgoing.pdf", width = 10, height = 5)
netVisual_bubble(cellchat,signaling = c("FN1","NEGR","SEMA3", "Adenosine"), #pathways,
                 remove.isolate = T, font.size = 10,font.size.title = 12,
                 sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"),thresh=0.01)
dev.off()

pdf("../../plot/fig6_netVisual_bubble_leftoutgoing.pdf", width = 10, height = 5)
netVisual_bubble(cellchat,signaling = c("IL2","DHT","ADGRB","CSF","BMP"), #pathways,
                 remove.isolate = T, font.size = 10,font.size.title = 12,
                 sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"),thresh=0.01)
dev.off()


netVisual_bubble(cellchat,signaling = c("CSF","MHC-I"), #pathways,
                 remove.isolate = T, font.size = 10,font.size.title = 12,
                 sources.use = c("MHCII"),thresh=0.01)

pdf("../../plot/fig6_signalingRole_network_Cell_death_outgoing.pdf", width = 10, height = 3)
pathways=c("IL2")
for (pathway.show in pathways) {
  netAnalysis_signalingRole_network(cellchat, signaling = pathway.show ,  width = 15, height = 2.5, font.size = 10) 
}
dev.off()

pdf("../../plot/fig6_signalingRole_network_Cell_cycle_outgoing.pdf", width = 10, height = 3)
pathways=c("DHT","ADGRB")
for (pathway.show in pathways) {
  netAnalysis_signalingRole_network(cellchat, signaling = pathway.show ,  width = 15, height = 2.5, font.size = 10) 
}
dev.off()

pdf("../../plot/fig6_signalingRole_network_EMT_outgoing.pdf", width = 10, height = 3)
pathways=c('FN1',"NEGR","SEMA3", "Adenosine")
for (pathway.show in pathways) {
netAnalysis_signalingRole_network(cellchat, signaling = pathway.show ,  width = 15, height = 2.5, font.size = 10) 
}
dev.off()

pdf("../../plot/fig6_signalingRole_network_MHCII_outgoing.pdf", width = 10, height = 3)
pathways=c("CSF","BMP")
for (pathway.show in pathways) {
  netAnalysis_signalingRole_network(cellchat, signaling = pathway.show ,  width = 15, height = 2.5, font.size = 10) 
}
dev.off()

netVisual_heatmap(cellchat, signaling = "FN1", color.heatmap = "Reds", title.name = "FN1 Signaling Strength")
netVisual_heatmap(cellchat, signaling = "IL6", color.heatmap = "Reds", title.name = "IL6 Signaling Strength")


pdf("../../plot/fig6_netVisual_aggregate_EMT_outgoing.pdf", width = 6, height = 6)
netVisual_aggregate(cellchat, signaling = 'FN1',layout = "chord",remove.isolate = F,
                    # color.use = bar,
                    sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"),
                    # targets.use = lym_cells
                    )

netVisual_aggregate(cellchat, signaling = c('NEGR'),layout = "chord",remove.isolate = F, # color.use = bar,
                    sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"),
                    # targets.use = lym_cells
  )
netVisual_aggregate(cellchat, signaling = c('SEMA3'),layout = "chord",remove.isolate = F, # color.use = bar,
                    sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"),
                    # targets.use = lym_cells
)
netVisual_aggregate(cellchat, signaling = c('Adenosine'),layout = "chord",remove.isolate = F, # color.use = bar,
                    sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"),
                    # targets.use = lym_cells
)

dev.off()

pdf("../../plot/fig6_netVisual_aggregate_Cell_death_outgoing.pdf", width = 6, height = 6)
netVisual_aggregate(cellchat, signaling = 'IL2',layout = "chord",remove.isolate = F,
                    # color.use = bar,
                    sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"),
                      )
dev.off()
pdf("../../plot/fig6_netVisual_aggregate_Cell_cycle_outgoing.pdf", width = 6, height = 6)
netVisual_aggregate(cellchat, signaling = 'DHT',layout = "chord",remove.isolate = F,
                    # color.use = bar,
                    sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"))

netVisual_aggregate(cellchat, signaling = 'ADGRB',layout = "chord",remove.isolate = F,
                    # color.use = bar,
                    sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"))


dev.off()


pdf("../../plot/fig6_netVisual_aggregate_MHCII_outgoing.pdf", width = 6, height = 6)
netVisual_aggregate(cellchat, signaling = 'CSF',layout = "chord",remove.isolate = F,
                    # color.use = bar,
                    sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"))

netVisual_aggregate(cellchat, signaling = 'BMP',layout = "chord",remove.isolate = F,
                    # color.use = bar,
                    sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"))
dev.off()


netVisual_aggregate(cellchat, signaling = 'ICAM',layout = "chord",remove.isolate = F,
                    # color.use = bar,
                    sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"))

netVisual_aggregate(cellchat, signaling = 'NECTIN',layout = "chord",remove.isolate = F,
                    # color.use = bar,
                    sources.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"))


df = filter(df,rowSums(df[,-1])!=0)

df[df$PT_income!=0 & df$Cell_death_income == 0 & df$Cell_cycle_income==0 & df$EMT_income==0 & df$Stress_income == 0 & df$MHCII_income==0,]$pathway
df[df$PT_income==0 & df$Cell_death_income == 0 & df$Cell_cycle_income==0 & df$EMT_income==0 & df$Stress_income != 0 & df$MHCII_income==0,]$pathway

df[df$PT_income==0 & df$Cell_death_income != 0 & df$Cell_cycle_income==0 & df$EMT_income==0 & df$Stress_income == 0 & df$MHCII_income==0,]$pathway
df[df$PT_income==0 & df$Cell_death_income == 0 & df$Cell_cycle_income!=0 & df$EMT_income==0 & df$Stress_income == 0 & df$MHCII_income==0,]$pathway
#  "ADGRE" "TGFb"  "GRN"   
df[df$PT_income==0 & df$Cell_death_income == 0 & df$Cell_cycle_income==0 & df$EMT_income!=0 & df$Stress_income == 0 & df$MHCII_income==0,]$pathway
#"VISFATIN"  "NEGR"
df[df$PT_income==0 & df$Cell_death_income == 0 & df$Cell_cycle_income==0 & df$EMT_income==0 & df$Stress_income == 0 & df$MHCII_income!=0,]$pathway

##接收信号
pdf("../../plot/fig6_netVisual_bubble_incoming.pdf", width = 6, height = 4)
netVisual_bubble(cellchat,signaling = c("ADGRE","TGFb","GRN","VISFATIN","NEGR"),
                 # sources.use = c(), 
                 remove.isolate = T,
                 targets.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"),thresh=0.01)

dev.off()

pdf("../../plot/fig6_netVisual_aggregate_Cell_cycle_incoming.pdf", width = 7, height = 7)
netVisual_aggregate(cellchat, signaling = 'ADGRE',layout = "chord",remove.isolate = F,
                    # color.use = bar,
                    targets.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"))

netVisual_aggregate(cellchat, signaling = 'TGFb',layout = "chord",remove.isolate = F, # color.use = bar,
                    targets.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"))
netVisual_aggregate(cellchat, signaling = 'GRN',layout = "chord",remove.isolate = F,
                    # color.use = bar,
                    targets.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"))
dev.off()


pdf("../../plot/fig6_netVisual_aggregate_EMT_incoming.pdf", width = 7, height = 7)
netVisual_aggregate(cellchat, signaling = 'VISFATIN',layout = "chord",remove.isolate = F,
                    # color.use = bar,
                    targets.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"))
netVisual_aggregate(cellchat, signaling = 'NEGR',layout = "chord",remove.isolate = F, # color.use = bar,
                    targets.use = c("PT","Cell_death","Stress","MHCII","EMT","Cell_cycle"))
dev.off()



