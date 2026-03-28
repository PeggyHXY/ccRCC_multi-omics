library(Seurat)
library(pheatmap)
library(ggplot2)
library(tidyr)
library(tidydr)

setwd('/mnt/data/home/haoxueyu/Project/Lixuwen-ccRCC_20250509/scRNAseq/results/P5/Rdata')
lym <- readRDS("lym.rds")
dim(lym) #33343 50793

unique(lym$celltype_refine)
levels(Idents(lym))
# [1] CD16+NK        CD16+NK_FCER1G CD8+T_IFNG     NKT            CD4+T_FOXP1    CD8+T_PARP8    CD4+Treg       CD8+Tex_HAVCR2 CD4+T_RGCC    
# [10] CD4+T_CCR6     CD56+NK        B_BANK1        CD4+T_FOS      CD8+ProT       CD8+Tex_CD74   CD8+Tact       Plasma_IGLC2+  B_AREG        
# [19] Plasma_IGLC2-  CD4+T_TSHZ2    CD8+Tex_LAG3  


colors=c("B_BANK1"="#7BAFDE","B_AREG"="#17BECFFF", "Plasma_IGLC2+"="#AEC7E8FF","Plasma_IGLC2-"= "#9EDAE5FF", 
           "CD16+NK_FCER1G"= "#D99BBD","CD16+NK"="#FF9896FF", "CD56+NK"="#F7B6D2FF","NKT"="#C49C94",
          "CD4+Treg" = "#E7BA52FF","CD4+T_FOXP1"="#FFBB78FF","CD4+T_CCR6"="#F79B30FF","CD4+T_RGCC"="#E7CB94FF","CD4+T_FOS"="#CCCC59", 
          "CD4+T_TSHZ2"="#F6EF8E",
          "CD8+Tex_HAVCR2"="#91D1C2FF","CD8+Tex_CD74"="#2CA02CFF","CD8+Tex_LAG3"="#7E6148FF",
         "CD8+ProT"="#CEDB9CFF","CD8+T_PARP8"="#299364","CD8+T_IFNG"= "#B5CF6BFF","CD8+Tact"="#8CA252FF"
          )

sample_palette=c("#1F77B4FF","#AEC7E8FF","#9EDAE5FF","#9467BDFF","#9C9EDEFF","#C5B0D5FF","#FF7F0EFF","#FFBB78FF","#E7CB94FF","#2CA02CFF","#98DF8AFF","#91D1C2FF","#D6616BFF","#FF9896FF", "#F7B6D2FF") #"#7F7F7FFF","#C7C7C7FF"

pid_palette=c("#51C3CC80","#C5B0D5FF","#FFBB78FF", "#98DF8AFF","#FF9896FF")
# 检查颜色数量
ct_order=c("B_BANK1","B_AREG","Plasma_IGLC2+","Plasma_IGLC2-","CD16+NK","CD16+NK_FCER1G","CD56+NK","NKT",
           "CD4+Treg","CD4+T_FOXP1","CD4+T_CCR6","CD4+T_RGCC","CD4+T_FOS","CD4+T_TSHZ2",
           "CD8+Tex_HAVCR2","CD8+Tex_CD74","CD8+Tex_LAG3","CD8+ProT","CD8+T_PARP8","CD8+T_IFNG","CD8+Tact")

lym@meta.data$celltype_refine <- factor(lym@meta.data$celltype_refine,levels = ct_order) 


prefix='../../plot/fig5_lym_'
l=0.03
p <- DimPlot(lym, reduction = "umap", group.by = 'celltype_refine',label = TRUE, label.size = 4,raster=FALSE, repel=TRUE, cols = colors) + 
  ggtitle("Cell Type")+
  theme_dr(xlength = l*4, ylength = l*4.8,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  # NoLegend()+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 1,hjust = 0,size = 8),
        plot.title = element_text(hjust = 0.5,face = 'bold')) 
p
ggsave(p, filename=paste0(prefix,"umap_celltype.pdf"),width=8, height=5)


lym$pid <- factor(lym$pid, levels = c('P6', 'P7', 'P8', 'P9', 'P11'))
p <- DimPlot(lym, reduction = "umap", group.by = 'celltype_refine',split.by = 'pid',label = FALSE, label.size = 5,raster=FALSE, repel=TRUE, cols = colors) + 
  ggtitle("Cell Type")+
  theme_dr(xlength = l*4, ylength = l*4.8,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  NoLegend()+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 1,hjust = 0,size = 8),
        plot.title = element_text(hjust = 0.5,face = 'bold')) 
p
ggsave(paste0(prefix,"umap_pid.pdf"),p,width=8,height=4)


lym$sample_id <- factor(lym$sample_id, levels = c("P6N1", "P6N2", "P6N3","P7N1","P7N2","P7N3","P8N1","P8N2","P8N3","P9N1", "P9N2","P9N3","P11N1", "P11N2", "P11N3"))
levels(lym$sample_id)
p <- DimPlot(lym, reduction = "umap", group.by = 'celltype_refine',split.by = 'sample_id',label = FALSE, label.size = 5,raster=FALSE, repel=TRUE, cols = colors, ncol=6) + 
  ggtitle("Cell Type")+
  theme_dr(xlength = l*4, ylength = l*4.8,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  NoLegend()+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 1,hjust = 0,size = 8),
        plot.title = element_text(hjust = 0.5,face = 'bold'))
p
ggsave(paste0(prefix,"umap_sample_id.pdf"),p,width=8,height=10)

#as_tibble()
unique(lym@meta.data[,c("leiden_res2.0", "celltype_refine")])
table(lym@meta.data[,c( "sample_id","celltype_refine")])

table(lym@meta.data[,c("celltype_refine")])

lym_marker_dict=list(
  'B cell'=c('CD79A','MS4A1','BANK1','AREG'), 
  'Plasma cell'=c('IGHA1','TNFRSF17','IGLC2','IGHG1'), #,'CD38'
  # 'NK cell'=c('GNLY', 'KLRD1', 'FGFBP2'),
  'CD16 NK'=c('GNLY','KLRD1','FCGR3A','FGFBP2','FCER1G'),#'EOMES'
  'CD56 NK'=c('NCAM1'),
  'NKT'=c('CD3D','CD3E','GNLY'),
  
  # "CD4+ T cell"=c('CD4','IL7R','PTPRC'),
  'CD4+ T'=c('FOXP1','RGCC'),#'FOS',
  'CD4+ Tn'=c('CD4','CCR7','SELL'),
  'CD4+ Tm'=c('GPR183','AREG','IL7R'),
  'CD4+ Th'=c('PTPN13','CCR6','CXCR3'),
  'CD4+ Tex'=c('NR4A1','PDCD1'),#'GADD45B'
  'CD4+ Treg'=c('FOXP3','IL2RA'),
  
  'CD8+ T cell'=c('CD8A', 'CD8B','TNFRSF9',"CD74"),
  # 'CD8+ cycling'=c('PCNA','MCM5','MKI67','TOP2A'),
  'CD8+ Tact'=c('GZMM','ICOS','CD69','CD38'), #activated 
  'CD8+ MAIT'=c('SLC4A10','KLRB1','ZBTB16'),#'TRAV1-2'
  'CD8+ Pro-T'=c('MKI67','TK1','STMN1'),
  'CD8+ Teff'=c('ZEB2','AOAH','IFNG','TOX'),
  'CD8+ Tex'=c('LAG3', 'HAVCR2','TCF7','TIGIT','TOX2','EOMES'), #exhausted
  'CD8+ Tm'=c('VCAM1','NAB1')
)
#TNFRSF9(CD137/4-1BB) TNFSF9(4-1BBL/CD137L)
# HAVCR2(TIM3) PDCD1(PD1) CD274(PDL1) TCF7(TCF1)  ITGAE(CD103)  CTLA4 
if(TRUE){

lym_marker_dict=list(
  'B cell'=c('CD79A','MS4A1','BANK1','AREG'), 
  'Plasma cell'=c('IGHA1','TNFRSF17','IGLC2','IGHG1'),#,'CD38'
  'CD16 NK'=c('GNLY','KLRD1','FCGR3A','FGFBP2','FCER1G','MYOM2'),
  'CD56 NK'=c('NCAM1'),
  'NKT'=c('CD3D','CD3E','GNLY'),
  
  'CD4+ Treg'=c('CD4','FOXP3','IL2RA'),
  'CD4+ Tn'=c('SELL','CCR7'),
  "CD4+ T cell"=c('IL7R','FOXP1','CCR6','RGCC','FOS'),#'PTPRC'
  'CD4+ Tex'=c('NR4A1','PDCD1','TSHZ2'),#'CTLA4'
  
  # 'CD8+ MAIT'=c('SLC4A10','KLRB1','ZBTB16'),#,'TRAV1-2' 
  'CD8+ T cell'=c('CD8A', 'CD8B','CD74','GZMK'),
  'CD8+ pTex'=c('TCF7','LEF1','STUB1'),#'CD28','EOMES','CD27'
  # 'CD8+Trm'=c('ITGAE'), #组织驻留
  'CD8+Tex_HAVCR2'=c('HAVCR2','TOX','TOX2','PDCD1','TNFRSF9','HLA-DRA'),#'VCAM1',
  'CD8+Tex_LAG3'=c('LAG3','ENTPD1','GZMB'),#exhausted 'TIGIT', ENTPD1(CD39)
  'CD8+Pro-T'=c('MKI67','TK1','STMN1'),
  'CD8+T_PARP8'=c('PARP8','ANKRD28','S100A4','BTD'),#'ITGA1',
  'CD8+T_IFNG'=c('IFNG','TNFSF9'),
  'CD8+Tact'=c('GZMM','ICOS','CD69') #activated 'CD38'
)
core_gene <-unlist(lym_marker_dict, use.names = FALSE)
lym$celltype_refine <- factor(lym$celltype_refine, levels = rev(ct_order))
levels(lym$celltype_refine)
options(repr.plot.width=12, repr.plot.height=5,repr.plot.res=300)
p <- DotPlot(lym, features = unique(core_gene), assay = 'RNA', group.by = "celltype_refine",col.min = -3,
             col.max = 3)+ 
  scale_color_gradient2( high = "#DF0000", mid = "white", low = "#51C3CC", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        legend.position = "bottom",  # 图例在下方
        legend.direction = "horizontal",  # 水平排列
        legend.box = "horizontal"         # 颜色和大小图例并排
  )
ggsave(paste0(prefix,'dotplot.pdf'),p, width = 12, height = 6)
p}

levels(lym@meta.data$celltype_refine)

p <- VlnPlot(lym, features= unique(core_gene), stack = TRUE, raster=FALSE, flip = TRUE, sort = TRUE,) + 
  theme(legend.position = "none")
ggsave(p, filename=paste0(prefix,"VlnPlot_core.pdf"), width = 10, height = 10)

checkpoint <- c('CD274','PDCD1', 'CTLA4', 'LAG3', 'TIGIT', 'HAVCR2', 'CD28', 'ICOS', 'CD40LG', 'TNFRSF4')
p <- VlnPlot(lym, features= unique(checkpoint), stack = FALSE, raster=FALSE, flip = TRUE, sort=TRUE) + theme(legend.position = "none")
p
ggsave(p, filename=paste0(prefix,"VlnPlot_checkpoint.pdf"), width = 8, height = 8)


#-----------------------------------------------------------------------------------------------------------------------
lym@meta.data$celltype_refine <- factor(
  lym@meta.data$celltype_refine,
  levels = rev(ct_order)
)
levels(lym$sample_id) 
p <- dittoBarPlot(lym, "celltype_refine", group.by = "sample_id", main = '', color.panel = colors)
p #+ guides(fill = guide_legend(reverse = TRUE))
ggsave(paste0(prefix,'dittoBarPlot.pdf'),p, width = 6, height = 5)


cell_counts <- table(lym@meta.data[,c('sample_id', 'celltype_refine')])
cell_prop <- prop.table(cell_counts, margin = 1) * 100  # margin=1表示按行标准化
cell_prop <- round(cell_prop, 2)
print(cell_prop)

write.table(cell_counts, file = "lym_celltype_counts.csv", row.names = TRUE, sep='\t',quote = F)
write.table(cell_prop, file = "lym_celltype_percent.csv", row.names = TRUE, sep='\t',quote = F)


# 计算距离矩阵
cell_prop_dist <- dist(cell_prop, method = "euclidean")

# 层次聚类
hc <- hclust(cell_prop_dist, method = "ward.D2")

# 绘制聚类树
plot(hc, hang = -1, main = "Sample Clustering by Cell Proportion")

p <- pheatmap::pheatmap(
  t(cell_prop),          # 转置使细胞类型为行，样本为列
  cluster_rows = TRUE,     # 对细胞类型聚类
  cluster_cols = hc,       # 使用层次聚类结果
  clustering_method = "ward.D2",
  scale = "row",           # 按行标准化（突出细胞类型比例差异）
  colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100),
  main = "Cell Proportion Heatmap",angle_col = "45"
)
ggsave(paste0(prefix,'CellProportionHeatmap_hc.pdf'),p, width = 6, height = 5)



p <- pheatmap::pheatmap(
  t(cell_prop),          # 转置使细胞类型为行，样本为列
  # cluster_rows = TRUE,     # 对细胞类型聚类
  # cluster_cols = hc,       # 使用层次聚类结果
  # clustering_method = "ward.D2",
  angle_col = "45",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "correlation",
  scale = "row",           # 按行标准化（突出细胞类型比例差异）
  colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100),
  main = "Cell Proportion Heatmap"
)
ggsave(paste0(prefix,'CellProportionHeatmap.pdf'),p, width = 6, height = 5)


p<- pheatmap::pheatmap(
  t(cell_prop),
  cellwidth = 15,     # 单元格宽度
  cellheight = 10,    # 单元格高度
  treeheight_row = 20, # 行聚类树高度
  treeheight_col = 20,  # 列聚类树高度
  scale = 'row',
  angle_col = "45",
  clustering_method = "ward.D2", #default: complete This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "correlation",
  color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100), #viridis::viridis(100),
  main = "Cell Proportion Heatmap"
)
ggsave(paste0(prefix,"CellProportionHeatmap_cor.pdf"),p,width=6,height=5)



library(ggplot2)

if(TRUE){
# 计算细胞比例（按样本和细胞类型分组）
prop_data <- as.data.frame(lym@meta.data) %>%
  group_by(sample_id, celltype_refine) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

desired_sample_order=c("P6N1", "P6N2", "P6N3","P7N1","P7N2","P7N3","P8N1","P8N2","P8N3","P9N1", "P9N2","P9N3","P11N1", "P11N2", "P11N3")
# 设置样本和细胞类型的顺序
prop_data$sample_id <- factor(prop_data$sample_id, levels = desired_sample_order)
prop_data$celltype_refine <- factor(prop_data$celltype_refine, levels = ct_order)

# 绘制条形图
p <- ggplot(prop_data, aes(x = sample_id, y = prop, fill = celltype_refine)) +
  
  geom_col(position = "stack", width = 0.7) +  # width控制柱宽
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    labels = scales::percent,  # y轴显示百分比
    expand = c(0, 0),         # 移除轴两端空白
    limits = c(0, 1.05),      # 扩展上限留白
    breaks = seq(0, 1, 0.2)   # 自定义刻度位置
  ) +
  labs(x = "Sample", y = "Proportion", fill = "Cell Type") +
  theme_minimal(base_size = 12) +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),  # 轴线
    axis.ticks = element_line(color = "black"),                 # 刻度线
    axis.ticks.length = unit(0.2, "cm"),                        # 刻度长度
    
    # x轴标签旋转
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      margin = margin(t = 5)  # 标签与轴的距离
    ),
    
    # 图例和背景
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = unit(c(1, 1, 2, 1), "cm")  # 下边距更大
  )
ggsave(paste0(prefix,'BarPlot.pdf'),p, width = 5, height = 5)
p}

#------------------------------------------------------------------------------------------------------------------------------------
# Integrated single-cell transcriptome and T cell receptor profiling reveals defects of T cell exhaustion in pulmonary tuberculosis
t_cell_marker = list("Naive"= c("CCR7", "SELL", "IL7R", "TCF7", "LEF1", "LTB", "FOXO1"), # "GZMK"
                     # 活化状态 (Activation)
                     "Activation"= c("SELL", "CD40LG", "ANXA1", "IL2RA", "CD69"),
                     # 记忆状态 (Memory)
                     "Memory"= c("CCR7", "TCF7", "CD69", "NR4A1", "MYADM", "GATA3", "TBX21"),
                     # 效应状态 (Effector)
                     "Effector"= c("BATF3", "IL2RA", "IFNG", "IRF4", "MYC", "SLC7A5", "SLC7A1", "XCL1", "GZMB", "CCL3", "CCL4", "IL2", "PRF1", "NKG7", "GNLY"),
                     # 耐受状态 (Tolerant) IZUMO1R (FOLR4)
                     "Tolerant"= c("EGR2", "IKZF2", "FOLR4", "CD200", "DGKZ", "BTLA", "TOX", "CTLA4"),
                     # 耗竭状态 (Exhausted) PDCD1 (PD-1) HAVCR2 (TIM3) CD244 (284) "IL10R",
                     "Exhausted"= c("PDCD1", "LAG3", "HAVCR2", "CD244", "CD160", "EOMES", "NR4A2", "PTGER4", "TOX", "TOX2", "TIGIT", "CTLA4", "ENTPD1","IL10RA"),
                     # 增殖状态 (Proliferation)
                     "Proliferation"= c("MKI67", "TK1", "STMN1"),
                     # 细胞毒性 (Cytotoxicity)
                     "Cytotoxicity"= c("GZMK", "GZMA", "GZMB", "NKG7", "PRF1", "IFNG", "GNLY"),
                     # 调节性功能 (Regulatory)
                     "Regulatory"= c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "TNFRSF18", "IL10", "IKZF2", "CCR8"),
                     # 早期记忆状态 (Early Memory)
                     "Early Memory"= c("CCR2", "CX3CR1", "IL18RAP", "ZEB2"))


t_sce <- subset(lym,subset=grepl("NKT|CD8|CD4", celltype_refine))
unique(t_sce$celltype_refine)

t_sce <- AddModuleScore_UCell(
  t_sce,assay='RNA',slot = 'data',
  features = t_cell_marker,
  ncores = 30, maxRank=2000
)

saveRDS(t_sce,'Tcell.rds')

t_sce <- readRDS('Tcell.rds')

tcell <- subset(t_sce, subset=grepl("CD8|CD4", celltype_refine))

# sample_scores <- t_sce@meta.data %>%
#   group_by(sample_id) %>%
#   summarise(across(
#     ends_with("_UCell"),
#     list(
#       mean = ~mean(., na.rm = TRUE)    # 均值
#       # median = ~median(., na.rm = TRUE) # 中位数
#       # total = ~sum(., na.rm = TRUE)      # 总和
#     ),
#     .names = "{.col}_{.fn}"  # 新列名格式：原列名_统计量
#   )) %>%
#   tibble::column_to_rownames("sample_id")  # 将sample_id设为行名

sample_scores <- tcell@meta.data %>%
  group_by(sample_id) %>%
  summarise(across(
    ends_with("_UCell"),
    list(
      mean = ~mean(., na.rm = TRUE)    # 计算均值
    ),
    .names = "{.col}" 
  )) %>% rename_with(~gsub("_UCell$", "", .x), .cols = ends_with("_UCell")) %>% tibble::column_to_rownames("sample_id")
 
write.table(sample_scores, file = "lym_sample_Ucell.csv", row.names = TRUE, quote = F, sep='\t')

prefix='../../plot/fig5_lym_' 
p <- pheatmap::pheatmap(
  t(sample_scores),
  cellwidth = 15,     # 单元格宽度
  cellheight = 15,    # 单元格高度
  treeheight_row = 20, # 行聚类树高度
  treeheight_col = 20,  # 列聚类树高度
  scale = 'row',
  angle_col = "45",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100), #viridis::viridis(100),
  main = "UCell Scores of T Cell State"
)
ggsave(paste0(prefix,"UCell_pheatmap_hc.pdf"),p,width=6,height=5)

p <- pheatmap::pheatmap(t(sample_scores), show_colnames = T, 
                        scale = "row",angle_col = "45",
                        color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100))
ggsave(paste0(prefix,"UCell_pheatmap.pdf"),p,width=5,height=5)


p <- pheatmap::pheatmap(
  t(sample_scores),
  cellwidth = 15,     # 单元格宽度
  cellheight = 15,    # 单元格高度
  treeheight_row = 20, # 行聚类树高度
  treeheight_col = 20,  # 列聚类树高度
  scale = 'row',
  angle_col = "45",
  clustering_method = "ward.D2", #default: complete This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "correlation",
  color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100), #viridis::viridis(100),
  main = "UCell Scores of T Cell State"
)
ggsave(paste0(prefix,"UCell_pheatmap_cor.pdf"),p,width=6,height=5)


sample_scores_dist <- dist(sample_scores, method = "euclidean")
# 层次聚类
hc <- hclust(sample_scores_dist, method = "ward.D2")

# 绘制聚类树
plot(hc, hang = -1, main = "Sample Clustering by UCell Score")


# p <- pheatmap::pheatmap(cor(t(sample_scores),method="pearson"), show_colnames = T, 
#                         # scale = "column",
#                         angle_col = "45",
#                         color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100))
# ggsave(paste0(prefix,"UCell_cor_pheatmap.pdf"),p,width=7,height=5)

library(reticulate)
use_condaenv('/mnt/data/home/haoxueyu/.conda/envs/scRNA')
library(SCP)

colnames(tcell@meta.data)

prefix='../../plot/fig5_Tcell_'

states=paste0(names(t_cell_marker),"_UCell")
states
comparisons=list(c("P6N1", "P6N2"), c("P6N2", "P6N3"),c("P6N1", "P6N3"),
                 c("P7N1", "P7N2"), c("P7N2", "P7N3"), c("P7N1", "P7N3"),
                 c("P8N1", "P8N2"), c("P8N2", "P8N3"), c("P8N1", "P8N3"),
                 c("P9N1", "P9N2"), c("P9N2", "P9N3"), c("P9N1", "P9N3"),
                 c("P11N1", "P11N2"), c("P11N2", "P11N3"),c("P11N1", "P11N3"))
for (ucell_col in states) {
  # 生成图形
  p <- FeatureStatPlot(
    tcell,
    stat.by = ucell_col,  # 动态指定当前 UCell Score 列
    group.by = "sample_id",
    pairwise_method = "wilcox.test",
    comparisons = comparisons,
    stat_color = sample_palette,
    legend.position = "none",
    title = ucell_col
  )
  
  # 保存图形（文件名基于列名）
  output_file <- paste0(prefix, ucell_col, ".pdf")
  ggsave(output_file, p, width = 6, height = 5)
  
  message("Saved: ", output_file)
}

tcell@meta.data$sample_id <- as.character(tcell@meta.data$sample_id)

p <- FeatureStatPlot(tcell, stat.by = 'Proliferation_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "sample_id", 
                     comparisons = list(c("P6N1", "P6N2"), c("P6N2", "P6N3"), c("P6N1", "P6N3"),
                                        c("P7N1", "P7N2"), c("P7N2", "P7N3"), c("P7N1", "P7N3"),
                                        c("P8N1", "P8N2"), c("P8N2", "P8N3"), c("P8N1", "P8N3"),
                                        c("P9N1", "P9N2"), c("P9N2", "P9N3"), c("P9N1", "P9N3"),
                                        c("P11N1", "P11N2"), c("P11N2", "P11N3"), c("P11N1", "P11N3")), #选择比较的细胞簇
                     stat_color=sample_palette, legend.position="none",sort="decreasing")
p
ggsave(paste0(prefix,"UCell_Proliferation_sort.pdf"),p,width=6,height=5)

p <- FeatureStatPlot(tcell, stat.by = 'Cytotoxicity_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "sample_id",   
                     comparisons = list(c("P6N1", "P6N2"), c("P6N2", "P6N3"),c("P6N1", "P6N3"),
                                        c("P7N1", "P7N2"), c("P7N2", "P7N3"), c("P7N1", "P7N3"),
                                        c("P8N1", "P8N2"), c("P8N2", "P8N3"), c("P8N1", "P8N3"),
                                        c("P9N1", "P9N2"), c("P9N2", "P9N3"), c("P9N1", "P9N3"),
                                        c("P11N1", "P11N2"), c("P11N2", "P11N3"),c("P11N1", "P11N3")), #选择比较的细胞簇
                     stat_color=sample_palette, legend.position="none",sort="decreasing" )

p
ggsave(paste0(prefix,"UCell_Cytotoxicity_sort.pdf"),p,width=6,height=5)

p <- FeatureStatPlot(tcell, stat.by = 'Exhausted_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "sample_id",   
                     comparisons = list(c("P6N1", "P6N2"), c("P6N2", "P6N3"),c("P6N1", "P6N3"),
                                        c("P7N1", "P7N2"), c("P7N2", "P7N3"), c("P7N1", "P7N3"),
                                        c("P8N1", "P8N2"), c("P8N2", "P8N3"), c("P8N1", "P8N3"),
                                        c("P9N1", "P9N2"), c("P9N2", "P9N3"), c("P9N1", "P9N3"),
                                        c("P11N1", "P11N2"), c("P11N2", "P11N3"),c("P11N1", "P11N3")), #选择比较的细胞簇
                     stat_color=sample_palette, legend.position="none",sort="decreasing" )
                     
p
ggsave(paste0(prefix,"UCell_Exhausted_sort.pdf"),p,width=6,height=5)

p <- FeatureStatPlot(tcell, stat.by = 'Regulatory_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "sample_id",   
                     comparisons = list(c("P6N1", "P6N2"), c("P6N2", "P6N3"),c("P6N1", "P6N3"),
                                        c("P7N1", "P7N2"), c("P7N2", "P7N3"), c("P7N1", "P7N3"),
                                        c("P8N1", "P8N2"), c("P8N2", "P8N3"), c("P8N1", "P8N3"),
                                        c("P9N1", "P9N2"), c("P9N2", "P9N3"), c("P9N1", "P9N3"),
                                        c("P11N1", "P11N2"), c("P11N2", "P11N3"),c("P11N1", "P11N3")), #选择比较的细胞簇
                     stat_color=sample_palette, legend.position="none",sort="decreasing" )
p
ggsave(paste0(prefix,"UCell_Regulatory_sort.pdf"),p,width=6,height=5)


p1 <- FeaturePlot(tcell, 
                  features = "Exhausted_UCell", 
                  cols = c("#51C3CC", "#DF0000"), # 低分灰色，高分红色
                  order = TRUE) & # order=TRUE让高分细胞点显示在最上层
  theme(plot.title = element_text(face = "bold", size = 16),
        legend.position = "right")
p1


#-------------------------------------------------------------------------------------
ct_scores <- tcell@meta.data %>%
  group_by(celltype_refine) %>%
  summarise(across(
    ends_with("_UCell"),
    list(
      mean = ~mean(., na.rm = TRUE)    # 均值
      # median = ~median(., na.rm = TRUE) # 中位数
      # total = ~sum(., na.rm = TRUE)      # 总和
    ),
    .names = "{.col}" 
  )) %>% rename_with(~gsub("_UCell$", "", .x), .cols = ends_with("_UCell")) %>% 
  tibble::column_to_rownames("celltype_refine")  # 将sample_id设为行名


prefix="../../plot/fig5_Tcellct_"
p <- pheatmap::pheatmap(t(ct_scores), show_colnames = T, cluster_rows = TRUE, 
                        cellwidth = 15,     # 单元格宽度
                        cellheight = 15,    # 单元格高度
                        treeheight_row = 20, # 行聚类树高度
                        treeheight_col = 20,  # 列聚类树高度
                        scale = "row",angle_col = "45",
                        color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100))
ggsave(paste0(prefix,"UCell_pheatmap.pdf"),p,width=5,height=3)


p <- pheatmap::pheatmap(
  t(ct_scores),
  scale = 'row',
  cluster_rows = TRUE, 
  cellwidth = 15,     # 单元格宽度
  cellheight = 15,    # 单元格高度
  treeheight_row = 20, # 行聚类树高度
  treeheight_col = 20,  # 列聚类树高度
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100), #viridis::viridis(100),
  main = "UCell Scores by Cell Type"
)
ggsave(paste0(prefix,"UCell_pheatmap_hc.pdf"),p,width=6,height=5)

p<- pheatmap::pheatmap(
  t(ct_scores),
  cellwidth = 15,     # 单元格宽度
  cellheight = 15,    # 单元格高度
  treeheight_row = 20, # 行聚类树高度
  treeheight_col = 20,  # 列聚类树高度
  scale = 'row',
  angle_col = "45",
  clustering_method = "ward.D2", #default: complete This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "correlation",
  color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100), #viridis::viridis(100),
  main = "UCell Scores of T Cell State"
)
ggsave(paste0(prefix,"UCell_pheatmap_cor.pdf"),p,width=6,height=5)


unique(tcell$celltype_refine)
table(lym@meta.data[,c('celltype_refine')])

levels(tcell@meta.data$celltype_refine) #查看是否为因子，可能会影响排序，看结果celltype未影响，sample_id有影响

tcell@meta.data$celltype_refine <- as.character(tcell@meta.data$celltype_refine)

library(scplotter)

p <- FeatureStatPlot(tcell, stat.by = 'Exhausted_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "celltype_refine",   
                     comparisons = list(c("CD8+Tex_LAG3", "CD8+Tex_HAVCR2"), c("CD8+Tex_HAVCR2", "CD8+Tex_CD74"),c("CD8+Tex_CD74", "CD4+T_TSHZ2"),
                                        c("CD4+T_TSHZ2", "CD8+ProT"),c("CD8+ProT", "CD8+Tact"), c("CD8+Tact", "CD4+Treg"), c("CD4+Treg","CD8+T_IFNG"),
                                        c("CD8+T_IFNG","CD8+T_PARP8"),c("CD8+T_PARP8","CD4+T_FOXP1"),c("CD4+T_FOXP1","CD4+T_RGCC"),
                                        c("CD4+T_RGCC","CD4+T_FOS"),c("CD4+T_FOS","CD4+T_CCR6")
                                        ),stat_color=colors,sort="decreasing", legend.position="none" )
                     
p
ggsave(paste0(prefix,"UCell_Exhausted.pdf"),p,width=6,height=5)

p <- FeatureStatPlot(tcell, stat.by = 'Cytotoxicity_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "celltype_refine",   
                     comparisons = list(c("CD8+Tex_LAG3", "CD8+Tex_HAVCR2"),c("CD8+Tex_HAVCR2", "CD8+Tex_CD74"),
                                        c("CD8+Tex_CD74", "CD8+ProT"),c("CD8+ProT", "CD8+Tact"),
                                        c("CD8+Tact", "CD8+T_IFNG"), c("CD8+T_IFNG", "CD8+T_PARP8"),
                                        c("CD8+T_PARP8","CD4+T_FOS"), c("CD4+T_FOS","CD4+T_CCR6"),c("CD4+T_CCR6","CD4+T_FOXP1"),c("CD4+T_FOXP1","CD4+T_RGCC"),c("CD4+T_RGCC","CD4+Treg"),
                                        c("CD4+Treg","CD4+T_TSHZ2"))
                     ,stat_color=colors,sort="decreasing", legend.position="none" )
p
ggsave(paste0(prefix,"UCell_Cytotoxicity.pdf"),p,width=6,height=5)


p <- FeatureStatPlot(tcell, stat.by = 'Proliferation_UCell',#目的基因   c("","")
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "celltype_refine",   
                     comparisons = list(c("CD8+ProT", "CD8+T_IFNG"), c("CD8+T_IFNG", "CD4+T_FOXP1"),c("CD4+T_FOXP1","CD8+T_PARP8"),
                                        c("CD8+T_PARP8","CD4+Treg"),c("CD4+Treg","CD8+Tex_HAVCR2"),c("CD8+Tex_HAVCR2","CD4+T_RGCC"),c("CD4+T_RGCC","CD4+T_CCR6"),
                                        c("CD4+T_CCR6","CD4+T_FOS"),c("CD4+T_FOS","CD8+Tex_CD74"),c("CD8+Tex_CD74","CD8+Tact"),c("CD8+Tact","CD4+T_TSHZ2"),c("CD4+T_TSHZ2","CD8+Tex_LAG3")
                     ), #选择比较的细胞簇 Treg
                     stat_color=colors, sort="decreasing", legend.position="none" )
p
ggsave(paste0(prefix,"UCell_Proliferation.pdf"),p,width=6,height=5)
