require(maftools)
library(ggplot2)
library(tidyr)
library(dplyr)
options(stringsAsFactors = F)
library(ComplexHeatmap)
## annovar

setwd('/Project/Lixuwen-ccRCC_20250509/WES/Results')
laml.maf = 'ccRCC_somatic_maftools.maf'
laml.clin = 'clinical_with_samples.csv'

laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

getClinicalData(laml)

tmb_res=tmb(laml,logScale = FALSE, captureSize = 60.507855)

tmb_res=tmb_res[, .(Tumor_Sample_Barcode ,total_perMB )]

tmb_vector <- setNames(tmb_res$total_perMB, tmb_res$Tumor_Sample_Barcode)
tmb_vector


col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(col) = c('Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins','Frame_Shift_Del',
               'In_Frame_Ins', 'In_Frame_Del', 'Splice_Site')
fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
fabcolors
colors= c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598","#ABDDA4","#66C2A5", "#3288BD")


mat=table(laml@data$Hugo_Symbol,laml@data$Tumor_Sample_Barcode)



# 读取MAF文件
maf_df <- read.delim(laml.maf, comment.char = "#", stringsAsFactors = FALSE)

# 读取临床数据
clinical_data <- read.delim(laml.clin, stringsAsFactors = FALSE)
# 将Sample_ID设置为行名，便于后续匹配
rownames(clinical_data) <- clinical_data$Tumor_Sample_Barcode

# 选择最感兴趣的基因（例如，已知的癌症基因）
genes_of_interest <- c('VHL','MTOR','TTN','PKHD1L1','OBSCN','ARHGEF12', 'ARID1A','BIRC6','CACNA1D','FAT3','FLT4', 'KDM5C','LRP1B','NF1','PBRM1','PIK3CA','PTEN','SETD2','SMARCA4','TRRAP')


# 构建矩阵：行是基因，列是样本，值是突变类型
mut_mat <- reshape2::acast(
  data = maf_df[maf_df$Hugo_Symbol %in% genes_of_interest, ],
  formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
  value.var = "Variant_Classification",
  fun.aggregate = function(x) {
    ifelse(length(x) > 0, paste(unique(x), collapse = ";"), '')
  }
)
# 将NA替换为空
mut_mat[is.na(mut_mat)] <- ""
mut_mat['TTN', c('P6N1','P6N2','P6N3')] <- 'Multi_Hit'
mut_mat['CACNA1D','P10N1'] <- 'Multi_Hit'

library(readr)

# write_tsv(as.data.frame.matrix(mut_mat),"mut_mat.tsv")

# 定义每种突变类型对应的绘图函数
alter_fun <- function(x, y, w, h, v) {
  # 使用 grid.rect 的向量化版本
  # 背景：为所有单元格先画一个白色背景
  grid.rect(x, y, w, h, 
            gp = gpar(fill = "white", col = "grey70"))
  
  # 错义突变 - 实心方块
  idx <- which(v == "Missense_Mutation")
  if (length(idx)) {
    grid.rect(x[idx], y[idx], w[idx] * 0.9, h[idx] * 0.9,
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  }
  
  # 无义突变 - 实心黑色方块
  idx <- which(v == "Nonsense_Mutation")
  if (length(idx)) {
    grid.rect(x[idx], y[idx], w[idx] * 0.9, h[idx] * 0.9,
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  }
  
  # 移码缺失 - 实心长条
  idx <- which(v == "Frame_Shift_Del")
  if (length(idx)) {
    grid.rect(x[idx], y[idx], w[idx] * 0.9, h[idx] * 0.1,
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  }
  
  # 移码插入 - 实心长条（不同颜色）
  idx <- which(v == "Frame_Shift_Ins")
  if (length(idx)) {
    grid.rect(x[idx], y[idx], w[idx] * 0.9, h[idx] * 0.1,
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  }
  
  # 框内缺失 - 空心矩形
  idx <- which(v == "In_Frame_Del")
  if (length(idx)) {
    grid.rect(x[idx], y[idx], w[idx] * 0.7, h[idx] * 0.7,
              gp = gpar(fill = "white", col = col["In_Frame_Del"], lwd = 2))
  }
  
  # 剪切位点 - 实心三角形
  idx <- which(v == "Splice_Site")
  if (length(idx)) {
    grid.polygon(
      x = unit.c(x[idx] - w[idx]*0.5, x[idx] + w[idx]*0.5, x[idx]),
      y = unit.c(y[idx] - h[idx]*0.5, y[idx] - h[idx]*0.5, y[idx] + h[idx]*0.5),
    gp = gpar(fill = col["Splice_Site"], col = NA)
    )
  }
}


if(TRUE){
# 确保临床数据的样本顺序与突变矩阵的列（样本）顺序一致
# sample_order <- colnames(mut_mat)
sample_order <- c('P6N1','P6N2','P6N3','P1N1','P1N2','P1N3','P7N1','P7N2','P7N3','P2N1','P2N2','P2N3',
                  'P4N1','P4N2','P4N3','P9N1','P9N2','P9N3','P10N1','P10N2','P10N3','P3N1','P3N2','P3N3',
                  'P5N1','P5N2','P5N3', 'P8N1','P8N2','P8N3')

clinical_annot <- clinical_data[sample_order, ]
clinical_annot <- clinical_annot %>% left_join(tmb_res, by = "Tumor_Sample_Barcode")

clinical_annot <- clinical_annot[match(sample_order, clinical_annot$Tumor_Sample_Barcode), ]

mut_mat <- mut_mat[, sample_order]

col <- c(
  Missense_Mutation = "#1F77B4FF",
  Nonsense_Mutation = "#FF7F0EFF",
  Frame_Shift_Ins = "#D62728FF",
  Frame_Shift_Del = "#9467BDFF",
  In_Frame_Ins = "#2CA02CFF",
  In_Frame_Del = "#BCBD22FF",
  Splice_Site = "#E377C2FF",
  Multi_Hit = "#C49C94FF"
)


# 创建顶部的注释（TMB是连续变量）
  alter_fun = list(
    background = alter_graphic("rect", fill = "#CCCCCC60"),  # 背景色
    Missense_Mutation = alter_graphic("rect", fill = col["Missense_Mutation"]),
    Nonsense_Mutation = alter_graphic("rect", fill = col["Nonsense_Mutation"]),
    Multi_Hit = alter_graphic("rect", fill = col["Multi_Hit"]),
    Frame_Shift_Ins = alter_graphic("rect", fill = col["Frame_Shift_Ins"]),
    Frame_Shift_Del = alter_graphic("rect", fill = col["Frame_Shift_Del"]),
    In_Frame_Ins = alter_graphic("rect", fill = col["In_Frame_Ins"]),
    In_Frame_Del = alter_graphic("rect", fill = col["In_Frame_Del"]),
    Splice_Site = alter_graphic("rect", fill = col["Splice_Site"])  
  )
col_pid=c('P1'="#AEC7E8FF", 'P2'="#DBDB8DFF", 'P3'="#F7B6D2FF",'P4'= "#FF9896FF",'P5'= "#C7C7C7FF", 
          'P6'="#51C3CC",'P7'="#C5B0D5FF",'P8'="#FFBB78FF",'P9'="#98DF8AFF",'P10'="#91D1C2FF")

clinical_annot$Patient=as.factor(clinical_annot$Patient)
ha <- HeatmapAnnotation(Patient=clinical_annot$Patient, col = list(Patient = col_pid), show_annotation_name = TRUE, annotation_name_side = "left", 
                        annotation_name_gp = gpar(fontsize = 10))
  
top_annotation <- HeatmapAnnotation(
  # 绘制TMB柱状图
  TMB = anno_barplot(
    clinical_annot$total_perMB,
    # height = unit(3, "cm"), # 设置高度
    bar_width = 0.8,
    border = FALSE,
    gp = gpar(fill = "#66C2A5", col = NA), # 设置颜色
    axis_param = list(at = c(0, 3, 5),  # 坐标轴刻度
      labels = c("0", "3", "5\n(mut/Mb)"),  # 刻度标签
      gp = gpar(fontsize = 10)),
    annotation_name_side = "left"  # 注释名称位置
  ),
  # 添加癌症类型标注
  Gender = clinical_annot$Gender,
  Stage = clinical_annot$Stage,
  Grade = clinical_annot$Grade,
  Location = clinical_annot$Location,
  scRNA.seq = clinical_annot$scRNA.seq,
  # 注释图例参数
  annotation_name_side = "left", 
  annotation_name_gp = gpar(fontsize = 10),
  show_annotation_name = TRUE,
  annotation_height = c(
    TMB = unit(3, "cm"),
    Grade = unit(0.9, "cm"),
    Stage = unit(1, "cm"),
    Gender = unit(1, "cm"),  
    Location = unit(1, "cm"),
    scRNA.seq = unit(1, "cm")
    
  ), 
  # 为分类变量设置颜色
  col = list(
    Gender = c("M" = "#BCBD22FF", "F" = "#C5B0D5FF"),
    Stage = c("T1b"="#AEC7E8FF", "T2a"="#17BECFFF", "T3a"="#1F77B4FF" ),
    Grade = c("G1-2" = "#DBDB8DFF","G2"="#9EDAE5FF", "G2-1" = "#C7C7C7FF", "G2-3" = "#C49C94FF", "G3-4"="#8C564BFF"),
    Location = c("L" = "#FFBB78FF", "R" = "#F7B6D2FF"),
    scRNA.seq=c("Y"="#98DF8AFF","N"="#FF9896FF")
  ),
  gap = unit(0.5, "mm") # 注释之间的间隔
)

onco_plot <- oncoPrint(
  mut_mat,
  top_annotation = top_annotation,  # 添加TMB等注释
  alter_fun = alter_fun,      # 自定义变异图形函数
  col = col,                  # 变异类型颜色映射
  row_names_side = "left",
  row_names_rot = 15,
  pct_side = "right",
  show_column_names = TRUE,
  column_names_rot = 45,
  column_order = sample_order,
  remove_empty_columns = FALSE,
  column_names_gp = gpar(fontsize = 10),#default 10
  column_title = "WES oncoplot",
  # cluster_columns= hclust,
  alter_fun_is_vectorized = FALSE,
  row_gap = unit(0.1, "mm"),
  column_gap = unit(0.1, "mm"),
  heatmap_width = unit(1.5, "npc"),
  heatmap_height = unit(1, "npc"),
  bottom_annotation = ha #注释信息在底部
)

# 绘制图形
pdf("plot/WES_oncoprint_final.pdf", width = 12, height = 8)
# png("plot/WES_oncoprint_final.png", width = 11, height = 8)
draw(onco_plot,annotation_legend_side = "right", heatmap_legend_side = "right", padding = unit(c(2, 2, 2, 15), "mm"),
     align_heatmap_legend = "global_center")
dev.off()
}
