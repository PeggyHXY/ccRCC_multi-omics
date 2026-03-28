library(MesKit)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
options(stringsAsFactors = F)
setwd('/Project/Lixuwen-ccRCC_20250509/WES/Results')
laml.clin='ccRCC_somatic_MesKit_clinical.txt'
clinical_data <- read.delim(laml.clin, stringsAsFactors = FALSE)
clinical_data
match.sa <- list('HRM'='P1','LRF'='P10','LYM'='P11','CXS'='P12','ZJY'='P13','LYS'='P2','JCH'='P3','BHY'='P4','ZHL'='P5'
,'LLF'='P6','ZSQ'='P7','QZX'='P8','CRS'='P9')

clinical_data$Tumor_ID=clinical_data$Tumor_Sample_Barcode

write.table(clinical_data, "ccRCC_somatic_MesKit_clinical.txt", quote=FALSE,row.names = FALSE,sep = "\t")

laml.maf='ccRCC_somatic_MesKit.maf'
maf_df <- read.delim(laml.maf, comment.char = "#", stringsAsFactors = FALSE)
maf_df$Tumor_ID=maf_df$Tumor_Sample_Barcode
maf_df$Patient_ID=maf_df$Patient
maf_df$Tumor_Sample_Label=str_sub(maf_df$Tumor_Sample_Barcode, -2, -1)

write.table(maf_df, "ccRCC_somatic_MesKit.maf", quote=FALSE,row.names = FALSE,sep = "\t")
########################################
# Ref_allele_depth,Alt_allele_depth,Tumor_ID,Patient_ID,Tumor_Sample_Label
# RefDepth,AltDepth

maf <- readMaf(mafFile = "./ccRCC_somatic_MesKit.maf",
               clinicalFile  = "./ccRCC_somatic_MesKit_clinical.txt",
               refBuild = "hg38") 

plotMutProfile(maf, use.tumorSampleLabel = TRUE,topGenesCount =20, class='SP',
               geneList = c('VHL','MTOR','TTN','PKHD1L1','OBSCN','ARHGEF12', 'ARID1A','BIRC6','CACNA1D','FAT3','FLT4', 'KDM5C','LRP1B','NF1','PBRM1','PIK3CA','PTEN','SETD2','SMARCA4','TRRAP')
            )

# pdf("PhyloTree.pdf", width = 5, height = 4)
# pid=c("BHY", "CRS", "HRM", "JCH", "LLF", "LRF", "LYS", "QZX", "ZHL", "ZSQ")
pid=c('P1','P2','P3','P4','P5','P6','P7','P8','P9','P10')
for (p in pid) {
  
  phyloTree <- getPhyloTree(maf, patient.id = p, method = "NJ", min.vaf = 0.05)
  plotree <- plotPhyloTree(phyloTree, patient.id = p, use.tumorSampleLabel = TRUE)
  print(plotree)
  ggsave(paste0(p,"_PhyloTree.pdf"),plotree, width = 5, height = 4)
  dev.off()
}

# dev.off()


for (p in pid) {
  pdf(paste0(p,"_PhyloTree.pdf"), width = 5, height = 4)
  tree.NJ <- getPhyloTree(maf, patient.id = p, method = "NJ")
  tree.MP <- getPhyloTree(maf, patient.id = p, method = "MP")
  # compare phylogenetic trees constructed by two approaches
  compareTree(tree.NJ, tree.MP, plot = TRUE, use.tumorSampleLabel = TRUE)
  dev.off()
}

phyloTree <- getPhyloTree(maf, patient.id = p, method = "NJ", min.vaf = 0.05)
# pdf(paste0(p,"_PhyloTree.pdf"), width = 5, height = 4)
plotree <- plotPhyloTree(phyloTree,patient.id = p, use.tumorSampleLabel = TRUE)
ggsave(paste0(p,"_PhyloTree.pdf"),plotree, width = 5, height = 4)

dev.off()


pdf(paste0(p,"_PhyloTree.pdf"), width = 5, height = 4)
tree.NJ <- getPhyloTree(maf, patient.id = p, method = "NJ")
tree.MP <- getPhyloTree(maf, patient.id = p, method = "MP")
# compare phylogenetic trees constructed by two approaches
compareTree(tree.NJ, tree.MP, plot = TRUE, use.tumorSampleLabel = TRUE)
dev.off()





