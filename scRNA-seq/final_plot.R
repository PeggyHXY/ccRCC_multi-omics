library(Seurat)
library(pheatmap)
library(ggplot2)
library(tidyr)
library(tidydr)

setwd('/mnt/data/home/haoxueyu/Project/Lixuwen-ccRCC_20250509/scRNAseq/results/P5/Rdata')

sce <- readRDS('ccRCC_annotated.rds')
dim(sce)
# 33343 110004

sample_palette=c("#1F77B4FF","#AEC7E8FF","#9EDAE5FF","#9467BDFF","#9C9EDEFF","#C5B0D5FF","#FF7F0EFF","#FFBB78FF","#E7CB94FF","#2CA02CFF","#98DF8AFF","#91D1C2FF","#D6616BFF","#FF9896FF", "#F7B6D2FF") #"#7F7F7FFF","#C7C7C7FF"

pid_palette=c("#51C3CC80","#C5B0D5FF","#FFBB78FF", "#98DF8AFF","#FF9896FF")

colors=c( "Epithelial cell" ="#D62728FF", "Endothelial cell"="#E377C2FF", "Fibroblast"="#9467BDFF", 
  "B cell"="#1F77B4FF","Plasma cell"="#AEC7E8FF","NK cell"="#BCBD22FF", "NKT cell"="#DBDB8DFF",  "CD4 T cell" = "#98DF8AFF", "CD8 T cell"="#2CA02CFF", "Dendritic cell" ="#8C564BFF", 
  "pDC"="#C49C94FF","Macrophage" ="#FF7F0EFF",  "Monocyte" ="#FFBB78FF", "Neutrophil cell"="#EEC362","Mast cell" ="#FF9896FF"
   )

ct_order=c('Epithelial cell','Endothelial cell','Fibroblast', 'B cell','Plasma cell','NK cell','NKT cell','CD4 T cell', 'CD8 T cell', 'Dendritic cell','pDC','Macrophage','Monocyte','Neutrophil cell','Mast cell')

Idents(sce) <- 'celltype_manual'
levels(Idents(sce))
# [1) "B cell"           "CD4 T cell"       "CD8 T cell"       "Dendritic cell"   "Endothelial cell" "Epithelial cell" 
# [7) "Fibroblast"       "Macrophage"       "Mast cell"        "Monocyte"         "NK cell"          "NKT cell"        
# [13) "Neutrophil cell"  "Plasma cell"      "pDC" 

# sce$celltype_manual <- factor(sce$celltype_manual, levels = names(Colclus))
marker_list <- list("Epithelial cell"=c('ALDH1A1', 'CD24','CA9','POU5F1','PDZK1IP1','KRT8'), #'NDUFA4L2' Fibroblast也表达 ,'VEGFA','CRYAB','STMN1'
                    'Endothelial cell'=c( 'PTPRB', 'KDR', 'PLVAP','PECAM1','VWF'),#'ESM1','CD34',
                    'Fibroblast'=c('ACTA2', 'COL1A1', 'COL1A2', 'TAGLN'),
                    'B cell'=c('MS4A1','CD79A'), 
                    "Plasma cell"=c('IGHA1','TNFRSF17'), 
                    'NK cell'=c('GNLY', 'KLRD1', 'FGFBP2'),
                    'CD4 T cell'=c('CD3D', 'CD3E','CD4'),
                    'CD8 T cell'=c('CD8A','CD8B','TNFRSF9'),
                    # 'γδ T cell'=c('TRGV9','TRDV2'),
                    'Dendritic cell'=c('CD1C','CLEC10A'),
                    'pDC'=c('SERPINF1',"CLEC4C"),
                    'Macrophage'=c('MSR1','MS4A7','CD14','MRC1'),#'CD68',
                    'Monocyte'=c('FCGR3A','FCN1','CSTA'), # 'CTSS'偏髓系的marker
                    'Neutrophil cell'=c('FCGR3B', 'CMTM2','CXCR2'), #'S100A8','S100A9'
                    'Mast cell'=c('MS4A2', 'TPSAB1','TPSB2')
                    )

all_marker_list <- list("Epithelial cell"=c('ALDH1A1', 'CD24','CA9','POU5F1','PDZK1IP1','KRT8'), #'NDUFA4L2' Fibroblast也表达 ,'VEGFA','CRYAB','STMN1'
                    'Endothelial cell'=c( 'PTPRB', 'KDR', 'PLVAP','PECAM1','VWF'),#'ESM1','CD34',
                    'Fibroblast'=c('ACTA2', 'COL1A1', 'COL1A2', 'TAGLN'),
                    'B cell'=c('MS4A1','CD79A'), 
                    "Plasma cell"=c('IGHA1','TNFRSF17'), 
                    'NK cell'=c('GNLY', 'KLRD1', 'FGFBP2'),
                    'CD4 T cell'=c('CD3D', 'CD3E','CD4'),
                    'CD8 T cell'=c('CD8A','CD8B','TNFRSF9'),
                    # 'γδ T cell'=c('TRGV9','TRDV2'),
                    'Dendritic cell'=c('CD1C','CLEC10A'),
                    'pDC'=c('LILRA4',"CLEC4C",'IRF7'),
                    'Macrophage'=c('MSR1','MS4A7','CD14'),
                    'Monocyte'=c('FCN1','CSTA'),
                    'Neutrophil cell'=c('FCGR3B', 'CMTM2','CXCR2','S100A8','S100A9'), 
                    'Mast cell'=c('MS4A2', 'TPSAB1','TPSB2','KIT','CPA3')
)


# core_gene <- c( marker_list$`Epithelial cell`, marker_list$Fibroblast)
core_gene <-unlist(marker_list, use.names = FALSE)

prefix='../../plot/fig3_'
# p=DimPlot(sce, group.by = c("celltype_manual"), label=T, label.size=4, reduction = "umap", raster=FALSE, repel=TRUE,
#           alpha=0.8) + scale_color_manual(values = colors)

l=0.03
sce$celltype_manual <- factor(sce$celltype_manual, levels = ct_order)
p <- DimPlot(sce, reduction = "umap", group.by = 'celltype_manual',label = TRUE, label.size = 5,raster=FALSE, repel=TRUE, cols = colors) + 
  ggtitle("Cell Type")+
  theme_dr(xlength = l*4, ylength = l*4.8,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  # NoLegend()+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 1,hjust = 0,size = 8),
        plot.title = element_text(hjust = 0.5,face = 'bold')) 
p
ggsave(p, filename=paste0(prefix,"umap_celltype.pdf"),width=6, height=5, dpi =250)

sce$pid <- factor(sce$pid, levels = c('P6', 'P7', 'P8', 'P9', 'P11'))
p <- DimPlot(sce, reduction = "umap", group.by = 'pid',label = FALSE, label.size = 5,raster=FALSE, repel=TRUE, cols = pid_palette) + 
  ggtitle("Patient")+
  theme_dr(xlength = l*4, ylength = l*4.8,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  # NoLegend()+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 1,hjust = 0,size = 8),
        plot.title = element_text(hjust = 0.5,face = 'bold')) 
p
ggsave(paste0(prefix,"broad_umap_bypid.pdf"),p,width=5,height=5)


sce$sample_id <- factor(sce$sample_id, levels = c("P6N1", "P6N2", "P6N3","P7N1","P7N2","P7N3","P8N1","P8N2","P8N3","P9N1", "P9N2","P9N3","P11N1", "P11N2", "P11N3"))
levels(sce$sample_id)
p <- DimPlot(sce, reduction = "umap", group.by = 'sample_id',label = FALSE, label.size = 5,raster=FALSE, repel=TRUE, cols = sample_palette) + 
  ggtitle("Sample")+
  theme_dr(xlength = l*4, ylength = l*4.8,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  # NoLegend()+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 1,hjust = 0,size = 8),
        plot.title = element_text(hjust = 0.5,face = 'bold'))
p
ggsave(paste0(prefix,"broad_umap_bysample_id.pdf"),p,width=5,height=5)


# p <- DimPlot(sce, reduction = "umap", group.by = 'celltype_manual',split.by = 'pid',label = FALSE, label.size = 5,raster=FALSE, repel=TRUE, cols = colors) + 
#   ggtitle("Cell Type")+
#   theme_dr(xlength = l*4, ylength = l*4.8,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
#   NoLegend()+
#   theme(panel.grid = element_blank(),axis.title = element_text(face = 1,hjust = 0,size = 8),
#         plot.title = element_text(hjust = 0.5,face = 'bold')) 
# p
# ggsave(paste0(prefix,"broad_umap_pid.pdf"),p,width=15,height=5)

# p <- UMAPPlot(sce, group.by = 'celltype_manual', split.by = 'sample_id', raster=FALSE) + scale_color_manual(values = colors)
# p <- DimPlot(sce, reduction = "umap", group.by = 'celltype_manual',split.by = 'sample_id',label = FALSE, label.size = 5,raster=FALSE, repel=TRUE, cols = colors) + 
#   ggtitle("Cell Type")+
#   theme_dr(xlength = l*4, ylength = l*4.8,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
#   NoLegend()+
#   theme(panel.grid = element_blank(),axis.title = element_text(face = 1,hjust = 0,size = 8),
#         plot.title = element_text(hjust = 0.5,face = 'bold'))
# p
# ggsave(paste0(prefix,"broad_umap_sample_id1.pdf"),p,width=15,height=5)

p <- VlnPlot(sce, features= unique(core_gene), stack = TRUE, raster=FALSE, flip = TRUE, sort = TRUE) + theme(legend.position = "none")
ggsave(p, filename=paste0(prefix,"VlnPlot_celltype.pdf"), width = 10, height = 10)


sce$celltype_manual <- factor(sce$celltype_manual, levels = rev(ct_order))
levels(sce$celltype_manual)
options(repr.plot.width=12, repr.plot.height=5,repr.plot.res=300)
p <- DotPlot(sce, features = core_gene, assay = 'RNA', group.by = "celltype_manual" )+ 
  scale_color_gradient2( high = "#DF0000", mid = "white", low = "#51C3CC", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        # legend.position = "bottom",  # 图例在下方
        # legend.direction = "horizontal",  # 水平排列
        # legend.box = "horizontal"         # 颜色和大小图例并排
        )
p
ggsave(paste0(prefix,'marker_dotplot.pdf'),p, width = 12, height = 5)


sce@meta.data$celltype_manual <- factor(
  sce@meta.data$celltype_manual,
  levels = rev(ct_order)
)
levels(sce$sample_id) 
p <- dittoBarPlot(sce, "celltype_manual", group.by = "sample_id", main = '', color.panel = colors)
p #+ guides(fill = guide_legend(reverse = TRUE))
ggsave(paste0(prefix,'dittoBarPlot.pdf'),p, width = 6, height = 5)


cell_counts <- table(sce@meta.data[,c('sample_id', 'celltype_manual')])
write.table(cell_counts, file = "Total_celltype_counts.csv", row.names = TRUE, sep='\t', quote = F)

cell_prop <- prop.table(cell_counts, margin = 1) * 100 
cell_prop <- round(cell_prop, 2)
print(cell_prop)
write.table(cell_prop, file = "Total_celltype_prop.csv", row.names = TRUE, sep='\t', quote = F)


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
  clustering_method = "ward.D2", angle_col = "45",
  scale = "row",           # 按行标准化（突出细胞类型比例差异）
  colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100),
  main = "Cell Proportion Heatmap"
)
ggsave(paste0(prefix,'CellProportionHeatmap.pdf'),p, width = 5, height = 4)



library(ggplot2)

# 计算细胞比例（按样本和细胞类型分组）
prop_data <- as.data.frame(sce@meta.data) %>%
  group_by(sample_id, celltype_manual) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

desired_sample_order=c("P6N1", "P6N2", "P6N3","P7N1","P7N2","P7N3","P8N1","P8N2","P8N3","P9N1", "P9N2","P9N3","P11N1", "P11N2", "P11N3")
# 设置样本和细胞类型的顺序
prop_data$sample_id <- factor(prop_data$sample_id, levels = desired_sample_order)
prop_data$celltype_manual <- factor(prop_data$celltype_manual, levels = ct_order)

# 绘制条形图
p <- ggplot(prop_data, aes(x = sample_id, y = prop, fill = celltype_manual)) +
  
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
  ) +
  # 添加水平参考线（可选）
  geom_hline(yintercept = seq(0, 1, 0.2), 
             color = "gray90", 
             linewidth = 0.2, 
             linetype = "dashed")

  # geom_bar(stat = "identity", position = "stack") +
  # scale_fill_manual(values = colors) +
  # labs(x = "Sample", y = "Proportion", fill = "Cell Type") + theme_minimal()+
  # theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  # 
p
ggsave(paste0(prefix,'broad_BarPlot.pdf'),p, width = 5, height = 5)


#############################################
library(UCell)
h.gene <- list('HALLMARK_HYPOXIA' = h_genesets$HALLMARK_HYPOXIA)
length(h.gene)

sce <- AddModuleScore_UCell(
  sce, assay='RNA',slot = 'data',
  features = h.gene ,
  ncores = 30, maxRank=2000
)

colnames(sce@meta.data)

library(reticulate)
use_condaenv('/mnt/data/home/haoxueyu/.conda/envs/scRNA')
library(SCP)
p <- FeatureStatPlot(sce, stat.by = 'HALLMARK_HYPOXIA_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "sample_id",   
                     comparisons = list(c("P6N1", "P6N2"), c("P6N2", "P6N3"),c("P6N1", "P6N3"),
                                        c("P7N1", "P7N2"), c("P7N2", "P7N3"), c("P7N1", "P7N3"),
                                        c("P8N1", "P8N2"), c("P8N2", "P8N3"), c("P8N1", "P8N3"),
                                        c("P9N1", "P9N2"), c("P9N2", "P9N3"), c("P9N1", "P9N3"),
                                        c("P11N1", "P11N2"), c("P11N2", "P11N3"),c("P11N1", "P11N3")), #选择比较的细胞簇
                     stat_color=sample_palette, legend.position="none") #sort="decreasing"
p
prefix='../../plot/fig3_'
ggsave(paste0(prefix,"UCell_HALLMARK_HYPOXIA.pdf"),p,width=6,height=5)


saveRDS(sce,file='ccRCC_annotated.rds')

sample_scores <- sce@meta.data %>%
  group_by(sample_id) %>%
  summarise(across(
    ends_with("_UCell"),
    list(
      mean = ~mean(., na.rm = TRUE),    # 均值
      median = ~median(., na.rm = TRUE) # 中位数
      # total = ~sum(., na.rm = TRUE)      # 总和
    ),
    .names = "{.col}_{.fn}"  # 新列名格式：原列名_统计量
  )) %>%
  tibble::column_to_rownames("sample_id") 

sample_scores

write.table(sample_scores, file = "sample_Ucell.csv", row.names = TRUE,quote = F, sep='\t')



