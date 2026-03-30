library(Seurat)
library(pheatmap)
library(ggplot2)
library(harmony)
library(tidyr)
library(tidydr)
library(SCP)
library(dplyr)

setwd('/mnt/data/home/haoxueyu/Project/Lixuwen-ccRCC_20250509/scRNAseq/results/P5/Rdata')
my <- readRDS("my.rds")
dim(my) #33343 19385

unique(my@meta.data[, c("pid", "donor","sample","sample_id")])
unique(my@meta.data$celltype_refine)

levels(Idents(my))
# [1] "TAM_VEGFA"  "cDC3"       "TAM_IL32"   "TAM_FOLR2"  "TAM_CD9"    "TAM_CD163"  "TAM_CCL3"   "TAM_BNC2"   "Mast cell"  "Mono_CD16"  "Mono_CD14" 
# [12] "Neutrophil" "pDC"        "cDC2"       "cDC1" 

colors=c("cDC1"="#AEC7E8FF","cDC2"="#17BECFFF", "cDC3"="#1F77B4FF","pDC"= "#9EDAE5FF", 
          "Mast cell"="#FFBB78FF", "Neutrophil"="#E7CB94FF", 
          "Mono_CD14"="#B5CF6BFF" ,"Mono_CD16"="#98DF8AFF", 
          "TAM_BNC2"="#DE9ED6FF","TAM_CCL3"="#E7969CFF","TAM_CD163"="#C5B0D5FF","TAM_TREM2"="#9467BDFF", 
          "TAM_FOLR2"="#F7B6D2FF","TAM_IL32"="#9C9EDEFF","TAM_FN1"="#C49C94FF"
          )

sample_palette=c("#1F77B4FF","#AEC7E8FF","#9EDAE5FF","#9467BDFF","#9C9EDEFF","#C5B0D5FF","#FF7F0EFF","#FFBB78FF","#E7CB94FF","#2CA02CFF","#98DF8AFF","#91D1C2FF","#8C564BFF","#C49C94FF",'#daa08f') #"#7F7F7FFF","#C7C7C7FF"

pid_palette=c("#AEC7E8FF","#C5B0D5FF","#FFBB78FF","#91D1C2FF","#C49C9490")
# 检查颜色数量

ct_order=c("cDC1","cDC2","cDC3", "pDC","Mast cell", "Neutrophil", "Mono_CD14","Mono_CD16", "TAM_BNC2","TAM_CCL3","TAM_CD163","TAM_TREM2","TAM_FOLR2", 
           "TAM_IL32","TAM_FN1")
my$celltype_refine <- factor(my$celltype_refine, levels = ct_order)

prefix='../../plot/fig5_my_'
l=0.03
p <- DimPlot(my, reduction = "umap", group.by = 'celltype_refine',label = TRUE, label.size = 5,raster=FALSE, repel=TRUE, cols = colors) + 
  ggtitle("Cell Type")+
  theme_dr(xlength = l*4, ylength = l*4.8,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  # NoLegend()+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 1,hjust = 0,size = 8),
        plot.title = element_text(hjust = 0.5,face = 'bold')) 
p
ggsave(p, filename=paste0(prefix,"umap_celltype.pdf"),width=6, height=4)


my$pid <- factor(my$pid, levels = c('P6', 'P7', 'P8', 'P9', 'P11'))
p <- DimPlot(my, reduction = "umap", group.by = 'pid', label = FALSE, label.size = 5,raster=FALSE, repel=TRUE, cols = pid_palette) + 
  ggtitle("Patient")+
  theme_dr(xlength = l*4, ylength = l*4.8,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  # NoLegend()+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 1,hjust = 0,size = 8),
        plot.title = element_text(hjust = 0.5,face = 'bold')) 
p
ggsave(paste0(prefix,"umap_pid.pdf"),p,width=4,height=4)


my$sample_id <- factor(my$sample_id, levels = c("P6N1", "P6N2", "P6N3","P7N1","P7N2","P7N3","P8N1","P8N2","P8N3","P9N1", "P9N2","P9N3","P11N1", "P11N2", "P11N3"))
levels(my$sample_id)
p <- DimPlot(my, reduction = "umap", group.by = 'sample_id',label = FALSE, label.size = 5,raster=FALSE, repel=TRUE, cols = sample_palette) + 
  ggtitle("Sample")+
  theme_dr(xlength = l*4, ylength = l*4.8,arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  # NoLegend()+
  theme(panel.grid = element_blank(),axis.title = element_text(face = 1,hjust = 0,size = 8),
        plot.title = element_text(hjust = 0.5,face = 'bold'))
p
ggsave(paste0(prefix,"umap_sample_id.pdf"),p,width=4,height=4)

colnames(my@meta.data)
df_plot = aggregate(my@meta.data[,c('nCounts_RNA','pct_counts_mt')], by = list(my$leiden_res1.0), FUN = mean)
df_plot = as.data.frame(df_plot)
colnames(df_plot)[1] <- 'Cluster'
pdf('../../plot/my_umap_nCount_percent_mt.pdf', width = 2.5, height = 2.6)
df_plot %>% 
  ggplot(aes(x = pct_counts_mt, y = nCounts_RNA,size=pct_counts_mt,color=pct_counts_mt)) +
  geom_point(alpha=0.8) +
  ggrepel::geom_text_repel(aes(label=Cluster), size=5, hjust=0.5, vjust=-1,color='black') +
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_gradient2(high = "#3A3A3A", mid = "white", low = "#A41125", midpoint = 10)+
  scale_x_continuous(expand = c(0.1, 0.1)) +
  scale_y_continuous(expand = c(0.1, 0.1)) 
dev.off()

# 'CHI3L1','CHI3L2','CHID1','TGM2'在巨噬细胞中几乎不表达
myel_marker_dict=list(
  'Myeloid'=c('CD68','CD163','LYZ'),
  'Mast cell'=c('MS4A2','TPSAB1','CTSG'),#
  'Neutrophil'=c('FCGR3B', 'CMTM2','CXCR2'),#
  'pDC'=c('SERPINF1',"CLEC4C",'STMN1'),#
  'CD1C+ DC'=c('CD1C','CLEC10A','CADM1','XCR1'),#
  'CLEC9A+ DC'=c('CLEC9A','SERPINF1','LSP1','COTL1'),#
  'Macrophage'=c('MSR1','MS4A7','CD14'),#,'AREG','IL1B'
  # 'M1'=c('FCGR1A','FCGR2A','FCGR3A','ITGAX','SOCS1','CXCL10','IFNG','TNF'),
  # 'M2'=c('MRC1','CD163','CD204','CD206','F13A1','STAB1','CHI3L1','CHI3L2','CHID1','TGM2','FN1','MSR1','TGFBI','OLR1'),
  # 'BAG3+ Macro'=c('BAG3','G0S2'),
  # 'MSR1+ TAM'=c('MSR1','STAT1'),
  'CD163+ TAM'=c('CD163','MAF','RNASE1'),#SEPP1
  'FOLR2+ TAM'=c('FOLR2','LYVE1','MARCO'),
  'GPNMB+ TAM'=c('GPNMB','CTSD','FTL'),#
  'VSIR+ TAM'=c('VSIR','STAT3','STAT6'),#
  'FN1+ TAM'=c('FN1','MARCO'),#
  'PLXDC2+ TAM'=c('PLXDC2','BMP2K'),#
  'SPP1+ TAM'=c('SPP1','TREM2','MRC1'),#
  'MHC II TAM'=c('APOE','APOC1'),#
  'Pro-Infla TAM'=c('MSR1','STAT1','CXCL10', 'CXCL9'),#
  'Monocyte'=c('FCN1','CSTA','S100A8','S100A9'),
  'CD16- Mono'=c('CD14'),#
  'CD16+ Mono'=c('FCGR3A'),#
  'cDC3'=c('CCL17', 'CLEC10A','CD74')
)
table(my$celltype_refine)
table(my$leiden_res1.0) #cluster17 表达Mastcell+macrophage marker，仅13个，归为mastcell

if(TRUE){
ct_order=c("cDC1","cDC2","cDC3","pDC", "Mast cell", "Neutrophil", "Mono_CD14","Mono_CD16", "TAM_BNC2","TAM_CCL3","TAM_TREM2","TAM_CD163","TAM_FOLR2", 
           "TAM_IL32","TAM_FN1")
  
my_marker_dict=list(
  # 'Myeloid'=c('CD68','CD163','LYZ'),
  'cDC1'=c('CLEC9A','SERPINF1','LSP1','CADM1','XCR1'),#,'COTL1'
  'cDC2'=c('CD1C','CLEC10A','FCER1A'),#
  'cDC3'=c('LAMP3','CCL19','CCR7'),#,'CD74','CCL17','NRP2',
  'pDC'=c('LILRA4',"CLEC4C",'IRF7'),#
  'Mast cell'=c('MS4A2','TPSAB1','KIT','CPA3'),#6,16 'CTSG',
  'Neutrophil'=c('FCGR3B', 'CMTM2','CXCR2'),#
  'Monocyte' =c('FCN1','CSTA'),#'S100A8','S100A9','CSF3R'
  'Mono_CD14'=c('CD14'),#
  'Mono_CD16'=c('FCGR3A'),#
  # 'Macrophage'=c('MSR1','MS4A7','PLXDC2','GPNMB'),#'IL1B','SEPP1'
  # 'TAM'=c('SPP1','LYVE1','CCL2'),#'CD9','CCR2','PDGFB','PDCD1LG2','TLR7',
  'M'=c('MSR1','LGMN','APOE','CD81'),#'MS4A7','C1QA', 'SPP1','LYVE1','CCL2'
  #'M1'=c('CD80','TNF')
  # 'M2'=c('CD200R1','ARG1'),
  'TAM_BNC2'=c('BNC2'),#14
  'TAM_CCL3'=c('CCL3'),#2
  'TAM_TREM2'=c('TREM2','CD9'),#0 'MARCO'
  'TAM_CD163'=c('CD163','MAF','RNASE1','IGF1','MRC1'),#1
  'TAM_FOLR2'=c('FOLR2'),#3,12 ,'CTSD','FTL'
  'TAM_IL32'=c('IL32','BCL11B'),#11 'FYN'
  'TAM_FNA'=c('FN1','VEGFA')#8
)

core_gene <-unlist(my_marker_dict, use.names = FALSE)
my$celltype_refine <- factor(my$celltype_refine, levels = rev(ct_order))
levels(my$celltype_refine)
# options(repr.plot.width=12, repr.plot.height=5,repr.plot.res=300)
p <- DotPlot(my, features = unique(core_gene), assay = 'RNA', group.by = "celltype_refine")+ 
  scale_color_gradient2( high = "#DF0000", mid = "white", low = "#51C3CC", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        legend.position = "bottom",  # 图例在下方
        legend.direction = "horizontal",  # 水平排列
        legend.box = "horizontal"         # 颜色和大小图例并排
  )
ggsave(paste0(prefix,'dotplot.pdf'),p, width = 10, height = 5)
p}
colnames(my@meta.data)

p <- DotPlot(my, features = unique(core_gene), assay = 'RNA', group.by = "leiden_res1.0")+ 
  scale_color_gradient2( high = "#DF0000", mid = "white", low = "#51C3CC", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        legend.position = "bottom",  # 图例在下方
        legend.direction = "horizontal",  # 水平排列
        legend.box = "horizontal"         # 颜色和大小图例并排
  )
p
ggsave(paste0(prefix,'dotplot_cluster.pdf'),p, width = 10, height = 5)

p <- VlnPlot(my, features= unique(core_gene), stack = TRUE, raster=FALSE, flip = TRUE, sort = TRUE,) + 
  theme(legend.position = "none")
ggsave(p, filename=paste0(prefix,"VlnPlot_core.pdf"), width = 10, height = 10)


table(my@meta.data[c('sample_id','celltype_refine')])

table(my@meta.data[c('celltype_refine')])

#-----------------------------------------------------------------------------------------------------------------------
my@meta.data$celltype_refine <- factor(
  my@meta.data$celltype_refine,
  levels = rev(ct_order)
)
levels(my$sample_id) 
p <- dittoBarPlot(my, "celltype_refine", group.by = "sample_id", main = '', color.panel = colors)
p #+ guides(fill = guide_legend(reverse = TRUE))
ggsave(paste0(prefix,'dittoBarPlot.pdf'),p, width = 6, height = 5)


cell_counts <- table(my@meta.data[,c('sample_id', 'celltype_refine')])
write.table(cell_counts, file = "myeloid_celltype_counts.csv", row.names = TRUE, sep='\t', quote = F)

cell_prop <- prop.table(cell_counts, margin = 1) * 100  # margin=1表示按行标准化
cell_prop <- round(cell_prop, 2)
print(cell_prop)
write.table(cell_prop, file = "myeloid_celltype_prop.csv", row.names = TRUE, sep='\t', quote = F)


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
  angle_col = "45",
  scale = "row",           # 按行标准化（突出细胞类型比例差异）
  colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100),
  main = "Cell Proportion Heatmap"
)
p
ggsave(paste0(prefix,'CellProportionHeatmap_hc.pdf'),p, width = 5, height = 4)

write.table(cell_counts, file = "my_celltype_counts.csv", row.names = TRUE, sep='\t')


p <- pheatmap::pheatmap(
  t(cell_prop),          # 转置使细胞类型为行，样本为列
  cellwidth = 15,     # 单元格宽度
  cellheight = 10,    # 单元格高度
  treeheight_row = 20, # 行聚类树高度
  treeheight_col = 20,  # 列聚类树高度
  angle_col = "45",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",#"correlation",
  scale = "row",           # 按行标准化（突出细胞类型比例差异）
  colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100),
  main = "Cell Proportion Heatmap"
)
ggsave(paste0(prefix,'CellProportionHeatmap.pdf'),p, width = 5, height = 4)


p <- pheatmap::pheatmap(
  t(cell_prop),          # 转置使细胞类型为行，样本为列
  cellwidth = 15,     # 单元格宽度
  cellheight = 10,    # 单元格高度
  treeheight_row = 20, # 行聚类树高度
  treeheight_col = 20,  # 列聚类树高度
  angle_col = "45",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "correlation",
  scale = "row",           # 按行标准化（突出细胞类型比例差异）
  colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100),
  main = "Cell Proportion Heatmap"
)
ggsave(paste0(prefix,'CellProportionHeatmap_cor.pdf'),p, width = 5, height = 4)

library(ggplot2)
my_cn <- table(my$sample_id, my$celltype_refine)

write.table(my_cn, file = "my_cellcount.csv", row.names = TRUE,quote = F, sep='\t')
# 计算细胞比例（按样本和细胞类型分组）
prop_data <- as.data.frame(my@meta.data) %>%
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
p
ggsave(paste0(prefix,'BarPlot.pdf'),p, width = 5, height = 5)

#-------------------------------------------------------------------------------------------------------------------------------------
library(UCell)
gene_set<- read.csv('mmc3_my.csv',sep='\t',header = T)##读取已经下载好的免疫细胞和对应基因列表，来源见文献附件
gene_set<-gene_set[, 1:2]#选取特异基因和对应的免疫细胞两行
score_list<- split(as.matrix(gene_set)[,1], gene_set[,2])

my <- AddModuleScore_UCell(
  my, assay='RNA',slot = 'data',
  features = score_list,
  ncores = 30, maxRank=2000
)

# saveRDS(my,'my.rds')

prefix='../../plot/fig5_my_'
sample_scores <- my@meta.data %>%
  group_by(sample_id) %>%
  summarise(across(
    ends_with("_UCell"),
    list(
      mean = ~mean(., na.rm = TRUE)    # 均值
      # median = ~median(., na.rm = TRUE) # 中位数
      # total = ~sum(., na.rm = TRUE)      # 总和
    ),
    .names = "{.col}"  # 新列名格式：原列名_统计量
  )) %>% rename_with(~gsub("_UCell$", "", .x), .cols = ends_with("_UCell")) %>% 
  tibble::column_to_rownames("sample_id")  # 将sample_id设为行名

write.table(sample_scores, file = "my_sample_Ucell.csv", row.names = TRUE,quote = F, sep='\t')

p <- pheatmap::pheatmap(t(sample_scores), show_colnames = T, 
                        cellwidth = 15,     # 单元格宽度
                        cellheight = 10,    # 单元格高度
                        treeheight_row = 20, # 行聚类树高度
                        treeheight_col = 20,  # 列聚类树高度
                        scale = "row",angle_col = "45",
                        color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100),
                        main = "UCell Scores of Myeloid Cells")
ggsave(paste0(prefix,"UCell_pheatmap.pdf"),p,width=7,height=4)

sample_scores_dist <- dist(sample_scores, method = "euclidean")
# 层次聚类
hc <- hclust(sample_scores_dist, method = "ward.D2")

# 绘制聚类树
plot(hc, hang = -1, main = "Sample Clustering by UCell Score")

p <- pheatmap::pheatmap(
  t(sample_scores),
  scale = 'row',
  cellwidth = 15,     # 单元格宽度
  cellheight = 10,    # 单元格高度
  treeheight_row = 20, # 行聚类树高度
  treeheight_col = 20,  # 列聚类树高度
  clustering_method = "ward.D2",
  angle_col = "45",
  color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100), #viridis::viridis(100),
  main = "UCell Scores of Myeloid Cells"
)
ggsave(paste0(prefix,"UCell_pheatmap_hc.pdf"),p,width=7,height=4)

p <- pheatmap::pheatmap(t(sample_scores), show_colnames = T,
                        cellwidth = 15,     # 单元格宽度
                        cellheight = 10,    # 单元格高度
                        treeheight_row = 20, # 行聚类树高度
                        treeheight_col = 20,  # 列聚类树高度
                        scale = "row",
                        angle_col = "45",
                        clustering_distance_rows = "euclidean",
                        clustering_distance_cols = "correlation",
                        color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100))
ggsave(paste0(prefix,"UCell_pheatmap_cor.pdf"),p,width=6,height=5)



#------------------------------------------------------------------------------------------------------------------------------
Mmarker1=list('M1'=c('FCGR1A','FCGR2A','FCGR3A','ITGAX','SOCS1','CXCL10'), 'M2'=c('MRC1','CD163','CD204','CD206','F13A1','STAB1','CHI3L1','CHI3L2','CHID1','TGM2','FN1','MSR1','TGFBI','OLR1'))

macro <- subset(my,subset=grepl("TAM", celltype_refine))

macro_marker<- read.csv('M1M2.marker.xls',sep='\t',header = T)
macro_marker<- as.list(macro_marker)
macro_marker$M2 <- macro_marker$M2[macro_marker$M2 != ""]
macro_marker$TAM <- macro_marker$TAM[macro_marker$TAM!= ""]
macro_marker

macro <- AddModuleScore_UCell(
  macro, assay='RNA',slot = 'data',
  features = macro_marker,
  ncores = 30, maxRank=2000
)

# saveRDS(macro,'macro.rds')
macro <- readRDS('macro.rds')

p <- VlnPlot(macro, features = c('M1_UCell',"M2_UCell","TAM_UCell"), group.by = "sample_id", cols=sample_palette,sort=TRUE)
ggsave(paste0(prefix,"UCell_Vlnplot_sample.pdf"),p,width=10,height=3)

p <- VlnPlot(macro, features = c('M1_UCell',"M2_UCell","TAM_UCell"), group.by = "celltype_refine", cols=colors,sort=TRUE)
ggsave(paste0(prefix,"UCell_Vlnplot_ct.pdf"),p,width=10,height=3)


prefix='../../plot/fig5_TAMct_'
ct_scores <- macro@meta.data %>%
  group_by(celltype_refine) %>%
  summarise(across(
    c('M1_UCell','M2_UCell',"TAM_UCell"),
    list(
      mean = ~mean(., na.rm = TRUE)    # 均值
      # median = ~median(., na.rm = TRUE) # 中位数
      # total = ~sum(., na.rm = TRUE)      # 总和
    ),
    .names = "{.col}" 
  )) %>% rename_with(~gsub("_UCell$", "", .x), .cols = ends_with("_UCell")) %>% 
  tibble::column_to_rownames("celltype_refine")  # 将sample_id设为行名

write.table(ct_scores, file = "TAM_sample_Ucell.csv", row.names = TRUE,quote = F, sep='\t')

p <- pheatmap::pheatmap(t(ct_scores), show_colnames = T, cluster_rows = TRUE, 
                        cellwidth = 15,     # 单元格宽度
                        cellheight = 15,    # 单元格高度
                        treeheight_row = 20, # 行聚类树高度
                        treeheight_col = 20,  # 列聚类树高度
                        scale = "row",angle_col = "45",
                        color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100))
ggsave(paste0(prefix,"UCell_pheatmap.pdf"),p,width=5,height=3)

sample_scores_dist <- dist(ct_scores, method = "euclidean")
# 层次聚类
hc <- hclust(sample_scores_dist, method = "ward.D2")
plot(hc, hang = -1, main = "Sample Clustering by UCell Score")

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
ggsave(paste0(prefix,"UCell_pheatmap_hc.pdf"),p,width=5,height=3)

p <- pheatmap::pheatmap(t(ct_scores), show_colnames = T,
                        cellwidth = 15,     # 单元格宽度
                        cellheight = 15,    # 单元格高度
                        treeheight_row = 20, # 行聚类树高度
                        treeheight_col = 20,  # 列聚类树高度
                        scale = "row",
                        angle_col = "45",
                        clustering_distance_rows = "euclidean",
                        clustering_distance_cols = "correlation",
                        color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100))
ggsave(paste0(prefix,"UCell_pheatmap_cor.pdf"),p,width=6,height=5)



colnames(macro@meta.data)
unique(macro@meta.data$celltype_refine)

p <- FeatureStatPlot(macro, stat.by = 'TAM_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "celltype_refine",   
                     comparisons = list(c("TAM_CD163", "TAM_IL32"), c("TAM_IL32", "TAM_FOLR2"),c("TAM_CD163", "TAM_FOLR2"),c("TAM_FOLR2", "TAM_CCL3"),
                                        c("TAM_CCL3", "TAM_TREM2"), c("TAM_TREM2", "TAM_FN1"), c("TAM_FN1","TAM_BNC2")
                                        ), #选择比较的细胞簇
                     stat_color=sample_palette,sort=TRUE, legend.position="none" )
p
ggsave(paste0(prefix,"UCell_TAM.pdf"),p,width=6,height=5)


p <- FeatureStatPlot(macro, stat.by = 'M1_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "celltype_refine",   
                     comparisons = list(c("TAM_CCL3", "TAM_FN1"), c("TAM_FN1", "TAM_CD163"),c("TAM_CCL3", "TAM_CD163"),
                                        c("TAM_CD163", "TAM_FOLR2"),
                                        c("TAM_FOLR2", "TAM_TREM2"), c("TAM_TREM2", "TAM_IL32"),
                                        c("TAM_IL32","TAM_BNC2")
                     ), #选择比较的细胞簇
                     stat_color=sample_palette,sort=TRUE, legend.position="none")
ggsave(paste0(prefix,"UCell_M1.pdf"),p,width=6,height=5)
p
p <- FeatureStatPlot(macro, stat.by = 'M2_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "celltype_refine",   
                     comparisons = list(c("TAM_CD163", "TAM_FOLR2"), c("TAM_FOLR2", "TAM_CCL3"),c("TAM_CD163", "TAM_CCL3"),
                                        c("TAM_CCL3","TAM_IL32"),c("TAM_IL32","TAM_FN1"),c("TAM_FN1", "TAM_TREM2"),c("TAM_TREM2","TAM_BNC2")
                     ), #选择比较的细胞簇
                     stat_color=sample_palette,sort=TRUE,legend.position="none" )
p
ggsave(paste0(prefix,"UCell_M2.pdf"),p,width=6,height=5)



#————————————————————————————————————————————————————————————————————————————————————————
prefix='../../plot/fig5_TAM_'
sample_scores <- macro@meta.data %>%
  group_by(sample_id) %>%
  summarise(across(
    c('M1_UCell','M2_UCell',"TAM_UCell"),
    list(
      mean = ~mean(., na.rm = TRUE)    # 均值
      # median = ~median(., na.rm = TRUE) # 中位数
      # total = ~sum(., na.rm = TRUE)      # 总和
    ),
    .names = "{.col}" 
  )) %>% rename_with(~gsub("_UCell$", "", .x), .cols = ends_with("_UCell")) %>% 
  tibble::column_to_rownames("sample_id")  # 将sample_id设为行名


p <- pheatmap::pheatmap(t(sample_scores), show_colnames = T, cluster_rows = TRUE, 
                        cellwidth = 15,     # 单元格宽度
                        cellheight = 15,    # 单元格高度
                        treeheight_row = 20, # 行聚类树高度
                        treeheight_col = 20,  # 列聚类树高度
                        scale = "row",
                        angle_col = "45",
                        color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100))
ggsave(paste0(prefix,"UCell_pheatmap.pdf"),p,width=4,height=3)

sample_scores_dist <- dist(sample_scores, method = "euclidean")
# 层次聚类
hc <- hclust(sample_scores_dist, method = "ward.D2")
plot(hc, hang = -1, main = "Sample Clustering by UCell Score")

p <- pheatmap::pheatmap(
  t(sample_scores),
  scale = 'row',
  cluster_rows = TRUE, 
  cellwidth = 15,     # 单元格宽度
  cellheight = 15,    # 单元格高度
  treeheight_row = 20, # 行聚类树高度
  treeheight_col = 20,  # 列聚类树高度
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100), #viridis::viridis(100),
  main = "UCell Scores by Sample"
)
ggsave(paste0(prefix,"UCell_pheatmap_hc.pdf"),p,width=5,height=3)


p <- pheatmap::pheatmap(t(sample_scores), show_colnames = T,
                        cellwidth = 15,     # 单元格宽度
                        cellheight = 15,    # 单元格高度
                        treeheight_row = 20, # 行聚类树高度
                        treeheight_col = 20,  # 列聚类树高度
                        scale = "row",
                        angle_col = "45",
                        clustering_distance_rows = "euclidean",
                        clustering_distance_cols = "correlation",
                        color = colorRampPalette(c("#51C3CC", "white", "#DF0000"))(100))
ggsave(paste0(prefix,"UCell_pheatmap_cor.pdf"),p,width=6,height=5)


mat <- AverageExpression(macro, features = unique(c(macro_marker$M1, macro_marker$M2, macro_marker$TAM)), group.by = "sample_id")$RNA

p <- pheatmap::pheatmap(as.matrix(mat), show_colnames = T, 
                        cluster_rows = FALSE, 
                        cluster_cols = TRUE, 
                        scale = "row",
                        angle_col = "45",
                        color = colorRampPalette(c("#1F77B4FF", "white", "#FF7F0EFF"))(150))
ggsave(paste0(prefix,"gene_pheatmap.pdf"),p,width=6,height=8)


library(reticulate)
use_condaenv('/mnt/data/home/haoxueyu/.conda/envs/scRNA')
library(SCP)

prefix='../../plot/fig5_TAM_'
p <- FeatureStatPlot(macro, stat.by = 'TAM_UCell',#目的基因  
                pairwise_method = "wilcox.test", #检验方法  
                group.by = "sample_id",   
                comparisons = list(c("P6N1", "P6N2"), c("P6N2", "P6N3"),c("P6N1", "P6N3"),
                                   c("P7N1", "P7N2"), c("P7N2", "P7N3"), c("P7N1", "P7N3"),
                                   c("P8N1", "P8N2"), c("P8N2", "P8N3"), c("P8N1", "P8N3"),
                                   c("P9N1", "P9N2"), c("P9N2", "P9N3"), c("P9N1", "P9N3"),
                                   c("P11N1", "P11N2"), c("P11N2", "P11N3"),c("P11N1", "P11N3")), #选择比较的细胞簇
                stat_color=sample_palette, legend.position="none" ,sort=TRUE)
ggsave(paste0(prefix,"UCell_TAM_sort.pdf"),p,width=6,height=5)

p <- FeatureStatPlot(macro, stat.by = 'M2_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "sample_id",   
                     comparisons = list(c("P6N1", "P6N2"), c("P6N2", "P6N3"),c("P6N1", "P6N3"),
                                        c("P7N1", "P7N2"), c("P7N2", "P7N3"), c("P7N1", "P7N3"),
                                        c("P8N1", "P8N2"), c("P8N2", "P8N3"), c("P8N1", "P8N3"),
                                        c("P9N1", "P9N2"), c("P9N2", "P9N3"), c("P9N1", "P9N3"),
                                        c("P11N1", "P11N2"), c("P11N2", "P11N3"),c("P11N1", "P11N3")), #选择比较的细胞簇
                     stat_color=sample_palette, legend.position="none",sort=TRUE )
ggsave(paste0(prefix,"UCell_M2_sort.pdf"),p,width=6,height=5)

p <- FeatureStatPlot(macro, stat.by = 'M1_UCell',#目的基因  
                     pairwise_method = "wilcox.test", #检验方法  
                     group.by = "sample_id",   
                     comparisons = list(c("P6N1", "P6N2"), c("P6N2", "P6N3"),c("P6N1", "P6N3"),
                                        c("P7N1", "P7N2"), c("P7N2", "P7N3"), c("P7N1", "P7N3"),
                                        c("P8N1", "P8N2"), c("P8N2", "P8N3"), c("P8N1", "P8N3"),
                                        c("P9N1", "P9N2"), c("P9N2", "P9N3"), c("P9N1", "P9N3"),
                                        c("P11N1", "P11N2"), c("P11N2", "P11N3"),c("P11N1", "P11N3")), #选择比较的细胞簇
                     stat_color=sample_palette, legend.position="none",sort=TRUE )
ggsave(paste0(prefix,"UCell_M1_sort.pdf"),p,width=6,height=5)



