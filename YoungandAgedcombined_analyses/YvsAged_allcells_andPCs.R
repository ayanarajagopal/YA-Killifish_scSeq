library(Seurat)
library(scater)
library(scuttle)
library(BiocParallel)
library(ggplot2)
install.packages("rliger")
library(EnhancedVolcano)
library(tidyverse)
library(patchwork)
library(dplyr)
library(scCustomize)
library(limma)
library(dplyr) 
install.packages("enrichR")
library(ggpubr)
library(MAST)
library(scclusteval)
library(SingleCellExperiment)
library(fgsea)
library(Matrix)
library(edgeR)
library(writexl)
BiocManager::install("EnhancedVolcano")
BiocManager::install("scclusteval")
BiocManager::install("MAST")
BiocManager::install("DESeq2")
library("scater")
library("scran")
library("DropletUtils")
library("BiocFileCache")
BiocManager::install("Nebulosa")
library(Nebulosa)
install.packages("scuttle")
BiocManager::install(c("goseq","GO.db"))
install.packages("devtools")
devtools::install_github("neurorestore/Libra")
devtools::install_github("rpolicastro/scProportionTest")

source("C:/Users/u0129074/Documents/ScRNAseq_paper/withsample2_young/custom_seurat_functions.R")
rm(killi_youngandaged_integrated_final)
raw.data_YandA_merged = list(youngisoseq,youngisoseq_S2,agedisoseq)
raw.data_YandA_new
DefaultAssay(youngisoseq) <- "SCT"
DefaultAssay(youngisoseq_S2) <- "SCT"
DefaultAssay(agedisoseq) <- "SCT"
dim(agedisoseq)
DefaultAssay(killi_youngandaged_integrated_final) <- "SCT"
agedisoseq@meta.data$condition <- "AgedK"
youngisoseq@meta.data$condition <- "YoungK"
youngisoseq_S2@meta.data$condition <- "YoungK"

VlnPlot(killi_youngandaged_integrated_final, features = c("percent.mt","percent.ribo", "nGene","nUMI"), cols=c("turquoise","#F8766D"), group.by="condition", ncol = 4)
VlnPlot(killi_youngandaged_integrated_final, c("nFeature_RNA","nCount_RNA"),stack=T, flip=, split.by="condition", cols=c("turquoise","#F8766D"), pt.size = 0.1, ncol = 1, group.by = "seurat_clusters")

combined.features_YandA_merged <- SelectIntegrationFeatures(object.list = raw.data_YandA_merged, nfeatures = 2400)
raw.data_YandA_merged <- PrepSCTIntegration(object.list = raw.data_YandA_merged, anchor.features = combined.features_YandA_merged)
combined.features_YandA_merged

int.anchors_YandA_merged <- FindIntegrationAnchors(object.list = raw.data_YandA_merged, reduction="rpca",
                                            anchor.features = combined.features_YandA_merged)

View(int.anchors_YandA@anchor.features)
killi_youngandaged_integrated_final <- IntegrateData(anchorset = int.anchors_YandA_merged, normalization.method = "SCT")
?IntegrateData
killi_youngandaged_integrated_final <- subset(killi_youngandaged_integrated_final,idents=clusters_to_recluster)
DefaultAssay(killi_youngandaged_integrated_final) <- "RNA"
DefaultAssay(killi_youngandaged_integrated_final)
dim(killi_youngandaged_integrated_final)
table(Idents(killi_youngandaged_integrated_final))
killi_youngandaged_integrated_final$condition <- ifelse(killi_youngandaged_integrated_final$orig.ident == "Young_KILLI_isoseq", "YoungK",
                                                 ifelse(killi_youngandaged_integrated_final$orig.ident == "Young_KILLI_isoseq_S2", "YoungK",
                                                        "AgedK"))

killi_youngandaged_integrated_final <- ScaleData(killi_youngandaged_integrated_final, verbose = FALSE)
killi_youngandaged_integrated_final <- RunPCA(object = killi_youngandaged_integrated_final, verbose = FALSE)
ElbowPlot(killi_youngandaged_integrated_final, ndims=40)

# Determine the K-nearest neighbor graph
killi_youngandaged_integrated_final <- FindNeighbors(object = killi_youngandaged_integrated_final, 
                                               dims = 1:40)

# Determine the clusters for various resolutions                                
killi_youngandaged_integrated_final <- FindClusters(object = killi_youngandaged_integrated_final,
                                              resolution = 0.6)

Idents(object = killi_youngandaged_integrated_final) <- "integrated_snn_res.0.55"
killi_youngandaged_integrated_final$
Idents(killi_youngandaged_integrated_final)
Idents(killi_youngandaged_integrated_final) = killi_youngandaged_integrated_final$seurat_clusters
cluster_ids_0.6 <- Idents(killi_youngandaged_integrated_final)
cluster_ids_0.6
Idents(killi_youngandaged_integrated_final)
#Calc. of TSNE and UMAP
killi_youngandaged_integrated_final <- RunUMAP(object = killi_youngandaged_integrated_final,  reduction= "pca", dims = 1:20)
killi_youngandaged_integrated_final <- RunTSNE(object = killi_youngandaged_integrated_final, reduction= "pca", dims = 1:20)

overall_aged_upgenes <- c("CCL-C5A","ISG15 (2 OF 2)","BX572630.2","OLA.23920","NFU-G-1-010204","NFU-G-1-017699","NFU-G-1-007454","RNF213","ORLA-UAA","NFU-G-1-014949","B2M (1 OF 2)","LDLRAP1B.1","PVRL2L.1","NFU-G-1-019352","CD68","KRT18","SFTPB","OLA.9601","NFU-G-1-016915","NFU-G-1-002581","RASGEF1BB","ZC3HDC1L","NFU-G-1-021973","OLA.15992","HELZ2.1","PVRL2.2","GABBR1.1","FXYD6","NFU-G-1-010317","HMGB2B.1","CD74","FKBP5.1","EEF1DA","LSR","NFU-G-1-021000.1","BX005256.5","SNORA74","EPD","KRTCAP2","KRT5","FABP7A","NEBL.1","SNCGB","NFU-G-1-003229","TFA","OLA.10208","CASKIN1 (1 OF 4)","KIAA0247","TSC22D3","ANXA5B","PLCB1 (2 OF 2)","APOA2","PABPC1A","PB.6299.1","OCS-05354","CR396586.1","PSAP","BX088712.2.1","NFU-G-1-019948","MIR-944","BAHCC1.2","STMN1A","ADCYAP1R1.1","PPP3R2 (2 OF 3).2","NFU-G-1-019981","6-SEP","RBP4","AHSG","BCL11A (1 OF 3)","LPHN2 (2 OF 2)","TNRC6B (2 OF 2)")
overall_aged_downgenes <- c("HSD17B12B.1","DBI","HSPBP1","SAFB.1","CADM2A","HSPB1","SERBP1 (3 OF 3)","NYAP2","DUSP4","HMGB2B","GGCTA","NFU-G-1-013518","SHC1","BASP1","RRAD","PCDH2G16","DUSP5","ATAD2B.1","HES5 (2 OF 2)","ITM2CB","5-8S-RRNA","FOS","PPP1R15B","BAG3","PLPPR4","NR4A1 (2 OF 2)","SLC25A25A","DNAJA4","HSPA1A.1","HSPA1A","MKNK2","SLC25A6","ZFAND2A","ATF3","HSPA8","OLA.5669","NFU-G-1-023540","09-MAR","NPAS4","DNAJB1B","03-SEP.1","NFU-G-1-026687","FOSL1","PB.4162.1","08-MAR","NFU-G-1-026237","NFU-G-1-006772","LOC107373011","FOSB","DNAJB2","HSP90AA1.2","06-SEP","LOC107396491","LOC107394792","LOC107379568","LOC107391080","LOC107373448","LOC107386688","LOC107374976","LOC107384710","LOC107388042","LOC107387014","LOC107388051","LOC107384722","LOC107381309","LOC107375033","LOC107372965","LOC107381735","LOC107383762","LOC107387769","LOC107394257","LOC107376678","NFIXA","LOC107391378","LOC107393582","LOC107378076","LOC107390365","LOC107393497","LOC107390924","LOC107386604","LOC107380342","LOC107397147","LOC107392070","LOC107374626","LOC107386157","CELF2","LOC107387181","LOC107372823","LOC107379492","LOC107374467","28S-RRNA","LOC107396109","SSU-RRNA-EUKARYA","LOC107377899","LOC107376386","LOC107394272","HSPA8.1","LOC107384901","LOC107379395","LOC107395803","LOC107373838","LOC107381616","LOC107381811","LOC107378566","LOC107383746","LOC107377705","LOC107396261","LOC107394573","LOC107387694","LOC107385042","LOC107394932")

#expression_table
killi_youngandaged_integrated_final <- NormalizeData(killi_youngandaged_integrated_final, verbose = FALSE)

killi_youngandaged_integrated_final_diffexpr_final2 <- FindAllMarkers(killi_youngandaged_integrated_final,min.pct = 0.25, only.pos = F, logfc.threshold = 0.25)
View(killi_youngandaged_integrated_final_diffexpr_final)
killi_youngandaged_integrated_top500_final2<- killi_youngandaged_integrated_final_diffexpr_final2 %>% group_by(cluster) %>% top_n(n=500, wt = avg_log2FC)
young_aged_top30_markers <- youngandaged_isoseq_combinedclusters %>% group_by(cluster) %>% top_n(n=30, wt = avg_log2FC)
View(killi_youngandaged_integrated_final_diffexpr_final2)
write_xlsx(killi_youngandaged_integrated_top500_final2,path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/y_a_markers_270224_res0.8_top500.xlsx")
write_xlsx(young_aged_top30_markers, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Aug2021/y_a_markers_top30_res1.txt", sep="\t", col.names=TRUE)
write.table(killi_youngandaged_integrated_top500_final2, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/y_a_markers_270224_res0.6_top502.xlsx", sep="\t", col.names=TRUE)

saveRDS(killi_youngandaged_integrated_final, file = "C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/aged/killi_young_aged_merged.rds")
saveRDS(onlyrevelantPCs, file = "C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/aged/killi_young_aged_PCs_merged.rds")

View(killi_youngandaged_integrated_final_diffexpr_final)
head(young_aged_isoseq_markers)
DefaultAssay(killi_youngandaged_integrated_final) <- "SCT"
####heatmap
top10_combined <- killi_youngandaged_integrated_final_diffexpr_final %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10

DoHeatmap(killi_youngandaged_integrated_final, features = killi_youngandaged_integrated_top10$gene,   draw.lines = TRUE,
          lines.width = NULL,hjust=0,angle=60, size=5,
          group.bar.height = 0.02) + scale_fill_gradientn(colors = c("skyblue","black","red"))
?DoHeatmap
##labels
cols_cols= c("#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#24B700","#00ACFC","#24B700","#F8766D","#24B700","#F8766D","#00ACFC","#00ACFC","#00ACFC","#F8766D","#24B700","#F8766D","#00ACFC","#F8766D","#24B700","#24B700","#00ACFC","#24B700")
killi_youngandaged_integrated_final$condition <- factor(killi_youngandaged_integrated_final$condition, levels = c("YoungK", "AgedK"))
head(killi_youngandaged_integrated_final$condition)
DimPlot(killi_youngandaged_integrated_final, reduction = "tsne", pt.size = 1,label.size=3,label=T, split.by="condition")
DimPlot(killi_youngandaged_integrated_final,label=F,col=cols_cols,reduction = "tsne", order=T,pt.size = 1.2) 
dev.off()
##Plot joint density of PCNA and glul genes; glial/nonglial signature##
plot_density(killi_youngandaged_integrated_final, joint=T,c("PCNA","SHH")) + theme_classic()
plot_density(pbmc, "CD4")
killi_youngandaged_integrated_final
genes.to.label.clusters <- c("WNT8B","ZIC2","PDGFRB","KLF2","SOX17","STAB2","PCNA","STMN1A","HMGB2A","CD9B","OLIG1","CX43","FABP7A","GLUL (2 OF 2)","SLC1A2B","APOEB","AIF1","LCP1","SLC17A6B","SV2A","ELAVL3","NEUROD2","MEX3A","TUBB5","EPD","APOA1")
genes.to.label.clusters.2 <- c("EPD","APOA1","WNT8B","ZIC2","PDGFRB","KLF2","SOX17","STAB2","PCNA","STMN1A","HMGB2A","CD9B","OLIG1","CX43","FABP7A","GLUL (2 OF 2)","SLC1A2B","APOEB","AIF1","LCP1","SLC17A6B","SV2A","ELAVL3","NEUROD2","MEX3A","TUBB5")

DotPlot(killi_youngandaged_integrated_final, split.by ="condition", features= genes.to.label.clusters.2, scale=T, cols=c("turquoise","#F8766D")) + RotatedAxis() + coord_flip() + theme(axis.text.y = element_text(size=12, color="black")) + theme(axis.text.x = element_text(size=12, color="black")) + theme(axis.text.y = element_text(size=12, color="black")) +  theme(axis.line.x = element_line(size=1, color="black")) + theme(axis.line.y = element_line(size=1, color="black"))
DefaultAssay(killi_youngandaged_integrated_final) <- "SCT"
DefaultAssay(killi_youngandaged_integrated_final)
Yvsaged_annotations <- c("Neu1","Neu2","Neu3","Neu4","Neu5","Prog1","Ex-mN1","MG","Prog2","Vas1","Vas2","Prog3","Neu6","OD","Neu7","Neu8","NE-RG3","Neu9","Neu10","Vas3")

names(youngPCs_annotations) <- levels(killi_youngandaged_integrated_final)
killi_youngandaged_integrated_final <- RenameIdents(killi_youngandaged_integrated_final, youngPCs_annotations)

killi_youngandaged_integrated_final$CellType <- Idents(killi_youngandaged_integrated_final)
killi_youngandaged_integrated_final$CellType
#####All cells stripchart prep_youngvs. aged diff expr#####
killi_youngandaged_integrated_final <- NormalizeData(killi_youngandaged_integrated_final, verbose=F)
dim(cluster0_DEGs)
cluster0_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "0_AgedK", ident.2 = "0_YoungK", verbose = FALSE)
View(cluster0_DEGs)
cluster0_DEGs$gene <- rownames(cluster0_DEGs)
cluster0_DEGs$clusternumber <- 0
cluster1_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "1_AgedK", ident.2 = "1_YoungK", verbose = FALSE)
head(cluster1_DEGs)
cluster1_DEGs$gene <- rownames(cluster1_DEGs)
cluster1_DEGs$clusternumber <- 1
cluster2_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "2_AgedK", ident.2 = "2_YoungK", verbose = FALSE)
head(cluster2_DEGs)
cluster2_DEGs$gene <- rownames(cluster2_DEGs)
cluster2_DEGs$clusternumber <- 2
cluster3_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "3_AgedK", ident.2 = "3_YoungK", verbose = FALSE)
head(cluster3_DEGs)
cluster3_DEGs$gene <- rownames(cluster3_DEGs)
cluster3_DEGs$clusternumber <- 3
cluster4_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "4_AgedK", ident.2 = "4_YoungK", verbose = FALSE)
head(cluster4_DEGs)
cluster4_DEGs$gene <- rownames(cluster4_DEGs)
cluster4_DEGs$clusternumber <- 4
cluster5_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "5_AgedK", ident.2 = "5_YoungK", verbose = FALSE)
head(cluster5_DEGs)
cluster5_DEGs$gene <- rownames(cluster5_DEGs)
cluster5_DEGs$clusternumber <- 5
cluster6_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "6_AgedK", ident.2 = "6_YoungK", verbose = FALSE)
head(cluster6_DEGs)
cluster6_DEGs$gene <- rownames(cluster6_DEGs)
cluster6_DEGs$clusternumber <- 6
cluster7_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "7_AgedK", ident.2 = "7_YoungK", verbose = FALSE)
head(cluster7_DEGs)
cluster7_DEGs$gene <- rownames(cluster7_DEGs)
cluster7_DEGs$clusternumber <- 7
cluster8_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "8_AgedK", ident.2 = "8_YoungK", verbose = FALSE)
head(cluster8_DEGs)
cluster8_DEGs$gene <- rownames(cluster8_DEGs)
cluster8_DEGs$clusternumber <- 8
cluster9_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "9_AgedK", ident.2 = "9_YoungK", verbose = FALSE)
head(cluster9_DEGs)
cluster9_DEGs$gene <- rownames(cluster9_DEGs)
cluster9_DEGs$clusternumber <- 9
cluster10_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "10_AgedK", ident.2 = "10_YoungK", verbose = FALSE)
head(cluster10_DEGs)
cluster10_DEGs$gene <- rownames(cluster10_DEGs)
cluster10_DEGs$clusternumber <- 10
cluster11_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "11_AgedK", ident.2 = "11_YoungK", verbose = FALSE)
head(cluster11_DEGs)
cluster11_DEGs$gene <- rownames(cluster11_DEGs)
cluster11_DEGs$clusternumber <- 11
cluster12_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "12_AgedK", ident.2 = "12_YoungK", verbose = FALSE)
head(cluster12_DEGs)
cluster12_DEGs$gene <- rownames(cluster12_DEGs)
cluster12_DEGs$clusternumber <- 12
cluster13_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "13_AgedK", ident.2 = "13_YoungK", verbose = FALSE)
head(cluster13_DEGs)
cluster13_DEGs$gene <- rownames(cluster13_DEGs)
cluster13_DEGs$clusternumber <- 13
cluster14_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "14_AgedK", ident.2 = "14_YoungK", verbose = FALSE)
head(cluster14_DEGs)
cluster14_DEGs$gene <- rownames(cluster14_DEGs)
cluster14_DEGs$clusternumber <- 14
cluster15_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "15_AgedK", ident.2 = "15_YoungK", verbose = FALSE)
head(cluster15_DEGs)
cluster15_DEGs$gene <- rownames(cluster15_DEGs)
cluster15_DEGs$clusternumber <- 15
cluster16_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "16_AgedK", ident.2 = "16_YoungK", verbose = FALSE)
head(cluster16_DEGs)
cluster16_DEGs$gene <- rownames(cluster16_DEGs)
cluster16_DEGs$clusternumber <- 16
cluster17_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "17_AgedK", ident.2 = "17_YoungK", verbose = FALSE)
head(cluster17_DEGs)
cluster17_DEGs$gene <- rownames(cluster17_DEGs)
cluster17_DEGs$clusternumber <- 17
cluster18_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "18_AgedK", ident.2 = "18_YoungK", verbose = FALSE)
head(cluster18_DEGs)
cluster18_DEGs$gene <- rownames(cluster18_DEGs)
cluster18_DEGs$clusternumber <- 18
cluster19_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "19_AgedK", ident.2 = "19_YoungK", verbose = FALSE)
head(cluster19_DEGs)
cluster19_DEGs$gene <- rownames(cluster19_DEGs)
cluster19_DEGs$clusternumber <- 19
cluster20_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "20_AgedK", ident.2 = "20_YoungK", verbose = FALSE)
head(cluster20_DEGs)
cluster20_DEGs$gene <- rownames(cluster20_DEGs)
cluster20_DEGs$clusternumber <- 20
cluster21_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "21_AgedK", ident.2 = "21_YoungK", verbose = FALSE)
head(cluster21_DEGs)
cluster21_DEGs$gene <- rownames(cluster21_DEGs)
cluster21_DEGs$clusternumber <- 21
cluster22_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "22_AgedK", ident.2 = "22_YoungK", verbose = FALSE)
head(cluster22_DEGs)
cluster22_DEGs$gene <- rownames(cluster22_DEGs)
cluster22_DEGs$clusternumber <- 22
cluster23_DEGs <- FindMarkers(killi_youngandaged_integrated_final,ident.1 = "23_AgedK", ident.2 = "23_YoungK", verbose = FALSE)
View(cluster23_DEGs)
cluster23_DEGs$gene <- rownames(cluster23_DEGs)
cluster23_DEGs$clusternumber <- 23
cluster23_DEGs
View(merged.table.1)
merged.table.1 <- Reduce(function(...) merge(..., all=T), list(cluster0_DEGs,cluster1_DEGs,cluster2_DEGs,cluster3_DEGs,cluster4_DEGs,cluster5_DEGs,cluster6_DEGs,cluster7_DEGs,cluster8_DEGs,cluster9_DEGs,cluster10_DEGs,cluster11_DEGs,cluster12_DEGs,cluster13_DEGs,cluster14_DEGs,cluster15_DEGs,cluster16_DEGs,cluster17_DEGs,cluster18_DEGs,cluster19_DEGs,cluster20_DEGs,cluster21_DEGs,cluster22_DEGs,cluster23_DEGs))
write_xlsx(cluster0_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_0.xls")
write_xlsx(cluster1_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_1.xls")
write_xlsx(cluster2_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_2.xls")
write_xlsx(cluster3_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_3.xls")
write_xlsx(cluster4_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_4.xls")
write_xlsx(cluster5_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_5.xls")
write_xlsx(cluster6_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_6.xls")
write_xlsx(cluster7_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_7.xls")
write_xlsx(cluster8_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_8.xls")
write_xlsx(cluster9_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_9.xls")
write_xlsx(cluster10_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_10.xls")
write_xlsx(cluster11_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_11.xls")
write_xlsx(cluster12_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_12.xls")
write_xlsx(cluster13_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_13.xls")
write_xlsx(cluster14_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_14.xls")
write_xlsx(cluster15_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_15.xls")
write_xlsx(cluster16_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_16.xls")
write_xlsx(cluster17_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_17.xls")
write_xlsx(cluster18_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_18.xls")
write_xlsx(cluster19_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_19.xls")
write_xlsx(cluster20_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_20.xls")
write_xlsx(cluster21_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_21.xls")
write_xlsx(cluster22_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_22.xls")
write_xlsx(cluster23_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/alldiffexpr_23.xls")
FeaturePlot(killi_youngandaged_integrated_final, features="MKI67","PCNA",cols=c("grey","#FA340D"),pt.size = 1.5, reduction = "tsne")
FeaturePlot(onlyrevelantPCs, features=c("MKI67"), cols=c("grey","#FA340D"),pt.size = 1.5)
?FeaturePlot

######pc makers testing#####
?FindMarkers
IC2_DEGs <- FindMarkers(onlyrevelantPCs, ident.1 = "Intercell.2_AgedK", ident.2 = "Intercell.2_YoungK", verbose = FALSE)
head(IC2_DEGs)
IC2_DEGs$gene <- rownames(IC2_DEGs)
IC3_DEGs <- FindMarkers(onlyrevelantPCs, ident.1 = "Intercell.3_AgedK", ident.2 = "Intercell.3_YoungK", verbose = FALSE)
head(IC3_DEGs)
IC3_DEGs$gene <- rownames(IC3_DEGs)
IC1_DEGs <- FindMarkers(onlyrevelantPCs, ident.1 = "Intercell.1_AgedK", ident.2 = "Intercell.1_YoungK", verbose = FALSE)
head(IC1_DEGs)
IC1_DEGs$gene <- rownames(IC1_DEGs)
IC4_DEGs <- FindMarkers(onlyrevelantPCs, ident.1 = "Intercell.4_AgedK", ident.2 = "Intercell.4_YoungK", verbose = FALSE)
head(IC4_DEGs)
IC4_DEGs$gene <- rownames(IC4_DEGs)
Astro1_DEGs <- FindMarkers(onlyrevelantPCs, ident.1 = "Astro-RG1_AgedK", ident.2 = "Astro-RG1_YoungK", verbose = FALSE)
head(Astro1_DEGs)
Astro1_DEGs$gene <- rownames(Astro1_DEGs)
Astro2_DEGs <- FindMarkers(onlyrevelantPCs, ident.1 = "Astro-RG2_AgedK", ident.2 = "Astro-RG2_YoungK", verbose = FALSE)
head(Astro2_DEGs)
RG4_DEGs$gene <- rownames(RG4_DEGs)
RG4_DEGs <- FindMarkers(onlyrevelantPCs, ident.1 = "EPD-RG4_AgedK", ident.2 = "EPD-RG4_YoungK", verbose = FALSE)
head(RG4_DEGs)
Astro2_DEGs$gene <- rownames(Astro2_DEGs)
NGP2_DEGs <- FindMarkers(onlyrevelantPCs, ident.1 = "NGP.2_AgedK", ident.2 = "NGP.2_YoungK", verbose = FALSE)
head(NGP2_DEGs)
NGP2_DEGs$gene <- rownames(NGP2_DEGs)
NGP2_DEGs <- FindMarkers(onlyrevelantPCs, ident.1 = "NGP.1_AgedK", ident.2 = "NGP.1_YoungK", verbose = FALSE)
head(NGP1_DEGs)
NGP1_DEGs$gene <- rownames(NGP1_DEGs)

IC1vsAst_DEGs <- FindMarkers(onlyrevelantPCs, ident.1 = "Intercell.1", ident.2 = "Astro-RG1", verbose = FALSE)
head(IC1vsAst_DEGs)
IC1vsAst_DEGs$gene <- rownames(IC1vsAst_DEGs)
table(Idents(onlyrevelantPCs))

write_xlsx(NERG3_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/alldiffexpr_NERG3.xls")
write_xlsx(IC2_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/alldiffexpr_IC2.xls")
write_xlsx(IC3_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/alldiffexpr_IC3.xls")
write_xlsx(IC1_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/alldiffexpr_IC1.xls")
write_xlsx(IC4_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/alldiffexpr_IC4.xls")
write_xlsx(Astro1_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/alldiffexpr_Astro1.xls")
write_xlsx(Astro2_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/alldiffexpr_Astro2.xls")
write_xlsx(RG4_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/alldiffexpr_RG4.xls")
write_xlsx(NGP1_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/alldiffexpr_NGP1.xls")
write_xlsx(NGP2_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/alldiffexpr_NGP2.xls")
write_xlsx(IC1vsAst_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/alldiffexpr_IC1vsAst_DEGs.xls")


PCs_filteredgenes <- read_excel("C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/allsubclusters_PC_DEGs.xlsx")

View(PCs_filteredgenes)

killi_youngandaged_integrated_PCs$seurat_clusters <- Idents(killi_youngandaged_integrated_PCs)
killi_youngandaged_integrated_PCs$seurat_clusters
table(Idents(killi_youngandaged_integrated_final))

##check the numbers again killi_youngandaged_integrated_PCs <- subset(killi_youngandaged_integrated_final,idents=c(3,7,12,15))
DefaultAssay(killi_youngandaged_integrated_PCs) <- "integrated"
rm(killi_youngandaged_integrated_PCs)
#PCs_integrated <- SCTransform(PCs_integrated, verbose = TRUE )
killi_youngandaged_integrated_PCs <- RunPCA(killi_youngandaged_integrated_PCs, verbose=FALSE)
ElbowPlot(killi_youngandaged_integrated_PCs, ndims = 40)
killi_youngandaged_integrated_PCs<- FindNeighbors(object = killi_youngandaged_integrated_PCs, 
                               dims = 1:20)

dim(PCs_integrated)
PCs_integrated$seurat_clusters
# Determine the clusters for various resolutions                                
killi_youngandaged_integrated_PCs <- FindClusters(object = killi_youngandaged_integrated_PCs,
                               resolution = 0.4)

#Calc. of TSNE/UMAP
killi_youngandaged_integrated_PCs <- RunTSNE(object = killi_youngandaged_integrated_PCs, reduction= "pca", dims=1:20)
killi_youngandaged_integrated_PCs <- RunUMAP(object = killi_youngandaged_integrated_PCs, reduction= "pca", dims = 1:20)
table(Idents(killi_youngandaged_integrated_PCs))
cDefaultAssay <- "RNA"

Idents()

col_new_ridgeplot_allPCs = c("lightsteelblue","#FAB693", "gray","#F4BB44", "#FFFF8F","#FDDA0D","#E4D00A", "#FA340D", "powderblue","darkcyan")
col_new_ridgeplot_newPCs = c("#F4BB44","#FDDA0D","#FFFF8F","#FAB693","darkcyan","#E4D00A","lightsteelblue", "gray", "powderblue","#FA340D")
relevantPCs_annotations <- c("Intercell.1","Intercell.3","Intercell.2","NGP.1","EPD-RG4","Intercell.4","NE-RG3","NGP.2","Astro-RG2","Astro-RG1")

DefaultAssay(onlyrevelantPCs)
onlyrevelantPCs <- subset(killi_youngandaged_integrated_PCs,idents=c(0,1,2,3,4,5,6,7,8,9))
DimPlot(onlyrevelantPCs, reduction = "tsne", pt.size = 1.8,label=T)
DimPlot(killi_youngandaged_integrated_PCs, reduction = "tsne", pt.size = 1.8,label=T, order=T, split.by="condition") 
p1+p2
onlyrevelantPCs
table(Idents(killi_youngandaged_integrated_PCs))
DefaultAssay(onlyrevelantPCs) <- "RNA"
killi_youngandaged_integrated_PCs <- NormalizeData(killi_youngandaged_integrated_PCs, verbose = FALSE)
DefaultAssay(killi_youngandaged_integrated_PCs)

PCs_YA_integrated_markers <- FindAllMarkers(killi_youngandaged_integrated_PCs, min.pct = 0.25, test.use="MAST", only.pos = F, logfc.threshold = 0.25)
DotPlot(killi_youngandaged_integrated_PCs, features= genes.to.label.clusters, scale=T, cols=c("green","purple")) + RotatedAxis() + theme(axis.text.y = element_text(size=12, color="black")) + theme(axis.text.x = element_text(size=10, color="black")) + theme(axis.text.y = element_text(size=20, color="black")) +  theme(axis.line.x = element_line(size=2, color="black")) + theme(axis.line.y = element_line(size=2, color="black"))

head(PCs_YA_integrated_markers)
youngagedPCs_PCs_cluster_top200 <- PCs_YA_integrated_markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=200, wt = avg_log2FC)
youngPCs_allS_cluster_top
youngPCs_top100genes <- PCs_integrated_clustersmarkers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=100, wt = avg_log2FC)
youngPCs_top100genes
youngPCs_top10genes <- PCs_integrated_clustersmarkers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=10, wt = avg_log2FC)
youngPCs_top10genes
Idents(onlyrevelantPCs) <- "CellType"
table(Idents(onlyrevelantPCs))

#########forsuppl info######
PCs_YA_integrated_markers_suppl <- FindAllMarkers(onlyrevelantPCs, min.pct = 0.25, test.use="MAST", only.pos = F, logfc.threshold = 0.25)
killi_youngandaged_integrated_top500_new222 <- PCs_YA_integrated_markers_suppl %>% group_by(cluster) %>% top_n(n=500, wt = avg_log2FC)
write_xlsx(killi_youngandaged_integrated_top500_new222, col_names = T, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/youngagedPCs_top500_feb2024.xlsx")

###Testing for known inflammaging genes####
overall_aged_upgenes <- c("CCL-C5A","ISG15 (2 OF 2)","BX572630.2","OLA.23920","NFU-G-1-010204","NFU-G-1-017699","NFU-G-1-007454","RNF213","ORLA-UAA","NFU-G-1-014949","B2M (1 OF 2)","LDLRAP1B.1","PVRL2L.1","NFU-G-1-019352","CD68","KRT18","SFTPB","OLA.9601","NFU-G-1-016915","NFU-G-1-002581","RASGEF1BB","ZC3HDC1L","NFU-G-1-021973","OLA.15992","HELZ2.1","PVRL2.2","GABBR1.1","FXYD6","NFU-G-1-010317","HMGB2B.1","CD74","FKBP5.1","EEF1DA","LSR","NFU-G-1-021000.1","BX005256.5","SNORA74","EPD","KRTCAP2","KRT5","FABP7A","NEBL.1","SNCGB","NFU-G-1-003229","TFA","OLA.10208","CASKIN1 (1 OF 4)","KIAA0247","TSC22D3","ANXA5B","PLCB1 (2 OF 2)","APOA2","PABPC1A","PB.6299.1","OCS-05354","CR396586.1","PSAP","BX088712.2.1","NFU-G-1-019948","MIR-944","BAHCC1.2","STMN1A","ADCYAP1R1.1","PPP3R2 (2 OF 3).2","NFU-G-1-019981","6-SEP","RBP4","AHSG","BCL11A (1 OF 3)","LPHN2 (2 OF 2)","TNRC6B (2 OF 2)")
overall_aged_downgenes <- c("HSD17B12B.1","DBI","HSPBP1","SAFB.1","CADM2A","HSPB1","SERBP1 (3 OF 3)","NYAP2","DUSP4","HMGB2B","GGCTA","NFU-G-1-013518","SHC1","BASP1","RRAD","PCDH2G16","DUSP5","ATAD2B.1","HES5 (2 OF 2)","ITM2CB","5-8S-RRNA","FOS","PPP1R15B","BAG3","PLPPR4","NR4A1 (2 OF 2)","SLC25A25A","DNAJA4","HSPA1A.1","HSPA1A","MKNK2","SLC25A6","ZFAND2A","ATF3","HSPA8","OLA.5669","NFU-G-1-023540","09-MAR","NPAS4","DNAJB1B","03-SEP.1","NFU-G-1-026687","FOSL1","PB.4162.1","08-MAR","NFU-G-1-026237","NFU-G-1-006772","LOC107373011","FOSB","DNAJB2","HSP90AA1.2","06-SEP","LOC107396491","LOC107394792","LOC107379568","LOC107391080","LOC107373448","LOC107386688","LOC107374976","LOC107384710","LOC107388042","LOC107387014","LOC107388051","LOC107384722","LOC107381309","LOC107375033","LOC107372965","LOC107381735","LOC107383762","LOC107387769","LOC107394257","LOC107376678","NFIXA","LOC107391378","LOC107393582","LOC107378076","LOC107390365","LOC107393497","LOC107390924","LOC107386604","LOC107380342","LOC107397147","LOC107392070","LOC107374626","LOC107386157","CELF2","LOC107387181","LOC107372823","LOC107379492","LOC107374467","28S-RRNA","LOC107396109","SSU-RRNA-EUKARYA","LOC107377899","LOC107376386","LOC107394272","HSPA8.1","LOC107384901","LOC107379395","LOC107395803","LOC107373838","LOC107381616","LOC107381811","LOC107378566","LOC107383746","LOC107377705","LOC107396261","LOC107394573","LOC107387694","LOC107385042","LOC107394932")


DotPlot(youngisoseq, features= genes.to.label.clusters, scale=T, cols=c("green","purple")) + coord_flip() + RotatedAxis() + theme(axis.text.y = element_text(size=15, color="black")) + theme(axis.text.x = element_text(size=30, color="black")) + theme(axis.text.y = element_text(size=20, color="black")) +  theme(axis.line.x = element_line(size=2, color="black")) + theme(axis.line.y = element_line(size=2, color="black"))

#####Find all markers
write_xlsx(youngagedPCs_PCs_cluster_top200, col_names = T, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/youngagedPCs_clubbed_top200genes_res0.5_fin.xlsx")
write_xlsx(youngPCs_top100genes_renamed, col_names = T, path="C:/Users/u0129074/Documents/rstudio_jan2020/20jan_afterintegration/Young_old_separate/Young_Isoseq_latest/July2020/young_aggr_oct2022/Analysis/youngPCs_top100genes_res0.5_new.xlsx")

DotPlot(killi_youngandaged_integrated_PCs, cluster.idents = F ,split.by="condition",idents=  c("Intercell.1","Intercell.3","Intercell.2","NGP.1","EPD-RG4","Intercell.4","NE-RG3","NGP.2","Astro-RG2","Astro-RG1"),features= c("ZIC2","WNT8B","SOX2","SOX3","HES5 (3 OF 9)","HES5 (2 OF 2)","HMGB2A","STMN1A","ASCL1B","GLUL (2 OF 2)","SLC1A2B","CX43","NDRG4","FABP7A","EPD","APOA1"), cols=c("turquoise","#F8766D")) + RotatedAxis()
youngPCs_annotations <- c("Intercell.1","Intercell.3","Intercell.2","NGP.1","EPD-RG4","Intercell.4","NE-RG3","NGP.2","Astro-RG2","Astro-RG1","some","thing")
names(youngPCs_annotations) <- levels(killi_youngandaged_integrated_PCs)
killi_youngandaged_integrated_PCs <- RenameIdents(killi_youngandaged_integrated_PCs, youngPCs_annotations)

killi_youngandaged_integrated_PCs$CellType <- Idents(killi_youngandaged_integrated_PCs)
killi_youngandaged_integrated_PCs$CellType
killi_youngandaged_integrated_PCs$seurat_clusters <- Idents(killi_youngandaged_integrated_PCs)
killi_youngandaged_integrated_PCs$seurat_clusters
table(Idents(onlyrevelantPCs))

table(onlyrevelantPCs$CellType)

Idents(object = onlyrevelantPCs) <- "RNA_snn_res.0.8"
DotPlot(killi_youngandaged_integrated_final, assay="SCT",split.by = "condition", features=c("SIX3","EMX2","SHH","ZIC2"))
FeaturePlot(onlyrevelantPCs, features=c("CX43"),cols=c("grey","powderblue"), pt.size = 1.5, reduction = "tsne")
FeaturePlot(onlyrevelantPCs, features=c("PCNA"), cols=c("grey","#FDDA0D"),pt.size = 1.5, reduction = "tsne")
FeaturePlot(onlyrevelantPCs, features=c("ZIC2"),cols=c("grey","steelblue"), pt.size = 1.5, reduction = "tsne")
FeaturePlot(onlyrevelantPCs, features=c("SLC1A2B","CX43","PCNA","ZIC2"), label=T, order=T, pt.size = 1, reduction = "tsne")

#####Finish subclustering#####
killi_youngandaged_integrated_final <- AddMetaData(killi_youngandaged_integrated_final, percent.mito_S2, col.name = "percent.mt")
?DimPlot
#normalize and setting RNA assay
DefaultAssay(killi_youngandaged_integrated_final) <- "RNA"
killi_youngandaged_integrated_final <- NormalizeData(killi_youngandaged_integrated_final, verbose = T)
killi_youngandaged_integrated$seurat_clusters
killi_youngandaged_integrated_final$seurat_clusters <- as.factor(killi_youngandaged_integrated_final$seurat_clusters)
head(killi_youngandaged_integrated)

ggplot(killi_youngandaged_integrated@meta.data, split.by="sample", aes(seurat_clusters))+geom_bar(stat="count")

YvsA.Ncbi.res0.4all <- c("ImN1","ImN2","Prog","ImN3","Ex-mN1","MG","Ih-mN1","Vas1","Vas2","Astro-RG2","Ex-mN2","OD","Ih-mN2","NE-RG3","Ih-mN3","Ex-mN3","Vas3")
relevantPCs_annotations <- c("Intercell.1","Intercell.3","Intercell.2","NGP.1","EPD-RG4","Intercell.4","NE-RG3","NGP.2","Astro-RG2","Astro-RG1")
relevantPCs_annotations_withcond <- c("NE-RG3_YoungK", "Intercell.3_YoungK", "EPD-RG4_YoungK", "Intercell.4_YoungK", "Intercell.2_YoungK", "Intercell.1_YoungK", "NGP.2_YoungK", "NGP.1_YoungK","Astro-RG1_YoungK","Astro-RG2_YoungK", "Astro-RG2_AgedK","Intercell.3_AgedK","Intercell.1_AgedK","NGP.2_AgedK","Intercell.2_AgedK","NE-RG3_AgedK", "NGP.1_AgedK","EPD-RG4_AgedK","Astro-RG1_AgedK","Intercell.4_AgedK")
names(relevantPCs_annotations) <- levels(onlyrevelantPCs)
onlyrevelantPCs <- RenameIdents(onlyrevelantPCs, relevantPCs_annotations)

onlyrevelantPCs$CellType <- Idents(onlyrevelantPCs)
onlyrevelantPCs$CellType
onlyrevelantPCs$seurat_clusters <- Idents(onlyrevelantPCs)
onlyrevelantPCs$seurat_clusters
table(Idents(onlyrevelantPCs))
onlyrevelantPCs

onlyrevelantPCs$Celltype.condition <- paste(Idents(onlyrevelantPCs), onlyrevelantPCs$condition, sep = "_")
onlyrevelantPCs$celltype <- Idents(onlyrevelantPCs)
onlyrevelantPCs$celltype
Idents(onlyrevelantPCs) <- "Celltype.condition"
onlyrevelantPCs

killi_youngandaged_integrated_final$Celltype.condition <- paste(Idents(killi_youngandaged_integrated_final), killi_youngandaged_integrated_final$condition, sep = "_")
killi_youngandaged_integrated_final$Celltype <- Idents(killi_youngandaged_integrated_final)
killi_youngandaged_integrated_final$Celltype
Idents(killi_youngandaged_integrated_final) <- "Celltype.condition"
killi_youngandaged_integrated_final
table(Idents(killi_youngandaged_integrated_final))
DefaultAssay(killi_youngandaged_integrated_final)
killi_youngandaged_integrated_final$integrated_snn_res.0.8


####running DEGs for all cell types
MG_DEGs <- FindMarkers(killi_youngandaged_integrated_final, ident.1 = "MG_AgedK", ident.2 = "MG_YoungK", verbose = FALSE)
head(MG_DEGs)
MG_DEGs$gene <- rownames(MG_DEGs)
write_xlsx(MG_DEGs, path="C:/Users/u0129074/Documents/ScRNAseq_paper/withsample2_young/alldiffexpr_MG.xls")
MG_DEGs_significant <- MG_DEGs %>% top_n(n=795, wt = p_val_adj)
View(MG_DEGs_significant)
readxl_example()
MG_filteredgenes <- read_xls(path="C:/Users/u0129074/Documents/ScRNAseq_paper/withsample2_young/alldiffexpr_MG.xls")
head(MG_filteredgenes)

####RUnning FGSEA#####
gene_sets

table(Idents(killi_youngandaged_integrated_final))
#####Running Harmony#####
integrated_harmony    <- merge(killi_young_integrated,agedisoseq)
integrated_harmony <- NormalizeData(integrated_harmony, verbose = F)
#integrated_harmony <- SCTransform(integrated_harmony,verbose = T, method = "glmGamPoi", vars.to.regress = "percent.mt")

integrated_harmony <- FindVariableFeatures(integrated_harmony,  nfeatures = 2000, verbose = F)
integrated_harmony <- ScaleData(integrated_harmony, verbose = F)
integrated_harmony <- RunPCA(integrated_harmony, verbose = F)
integrated_harmony <- RunUMAP(integrated_harmony, reduction = "pca", dims = 1:30, verbose = F)

#expression_table
killi_youngandaged_integrated_final_diffexpr_final <- FindAllMarkers(killi_youngandaged_integrated_final, min.pct = 0.25, only.pos = F, logfc.threshold = 0.25)
head(youngandaged_isoseq_combinedclusters)
killi_youngandaged_integrated_top500 <- killi_youngandaged_integrated_final_diffexpr_final %>% group_by(cluster) %>% top_n(n=500, wt = avg_log2FC)
young_aged_top30_markers <- youngandaged_isoseq_combinedclusters %>% group_by(cluster) %>% top_n(n=30, wt = avg_log2FC)
View(killi_youngandaged_integrated_final_diffexpr)
write_xlsx(killi_youngandaged_integrated_top500, path ="C:/Users/u0129074/Documents/ScRNAseq_paper/Aging paper/Extraanalyses/withsample2_young/y_a_markers_1809_res0.4_top500.xlsx")
write_xlsx(young_aged_top30_markers, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Aug2021/y_a_markers_top30_res1.txt", sep="\t", col.names=TRUE)
young_aged_top30_markers
head(young_aged_isoseq_markers)

dim(SCT_counts(killi_youngandaged_integrated_final))
counts(sce)[1:6, 1:6]

killi_youngandaged_integrated_final@assays
cluster_evaluation_result
filteredgenes_new <- read_xlsx(path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/stripchartdata_all23celltypes.xlsx")
onlytopgenes_allcells <- merged.table.1 %>%
  filter(avg_log2FC > 1.5 |avg_log2FC < -1.5) %>%
  select(p_val, avg_log2FC, gene, p_val_adj, clusternumber)

head(filteredgenes_new)
filteredgenes_PCs
View(onlytopgenes_allcells)
# Using 'label' and 'sample' as our two factors; each column of the output
# corresponds to one unique combination of these two factors.
summed_yanda <- scuttle::aggregateAcrossCells(sce_yanda, 
                                     id=colData(sce_yanda)[,c("CellType", "condition")])
?pseudobulk
metadata_n$cluster_id <- factor(killi_youngandaged_integrated_final@active.ident)
SCT_counts <- killi_youngandaged_integrated_final@assays$sct@counts 
sce_yanda <-  SingleCellExperiment(killi_youngandaged_integrated_final, assays = list(counts = SCT_counts), colData = metadata_n)
sce_yanda
sce_yanda <- sce_yanda[rowSums(counts(sce_yanda)) > 5,]
colData(sce_yanda)
counts(sce_yanda) <- as.matrix(counts(sce_yanda))
reduced_sce <- pseudobulk(sce_yanda, group_by = "condition")

Idents(killi_youngandaged_integrated_final) <- "condition"

metadata$cluster_id <- factor(seurat@active.ident) <- killi_youngandaged_integrated_final@meta.data
rm(SCT_counts)

FeatureScatter(object = killi_youngandaged_integrated_final,feature1 = "nGene")
DimPlot(killi_youngandaged_integrated_final, group.by="nGene",reduction="tsne")
?DimPlot

de_res <- test_de(fit, contrast = `stimstim` + `cellCD4 T cells:stimstim`, 
                  pseudobulk_by = paste0(stim, "-", ind)) 
killi_youngandaged_integrated
my.deg.new <- FindMarkers(killi_youngandaged_integrated, ident.1 = c("ImN","MN","MG","Inter-NGP","Vas1","Vas2","EPD","NA","RG","OD","NGP","Vas3"), ident.2 = NULL, test.use = "MAST",logfc.threshold = )
my.deg.poisson <- FindMarkers(killi_youngandaged_integrated_final, ident.1 = "AgedK", ident.2 = "YoungK", test.use = "poisson")
View(my.deg.poisson)
my.deg.new <- FindMarkers(killi_youngandaged_integrated_final, ident.1 = "AgedK", ident.2 = "YoungK", test.use = "MAST",logfc.threshold = 0.25)
head(my.deg.new)
Idents(killi_youngandaged_integrated_final)
my.deg <- FindMarkers(killi_youngandaged_integrated, ident.1 = c("ImN1","Ex-MN","MG","Inter-NGP","Vas1","Vas2","EPD","NA","RG","OD","NGP","Vas3"), ident.2 = NULL, only.pos = F, test.use = "MAST", min.pct = 0.4)
?FindMarkers
my.deg.next
killi_youngandaged_integrated_final$ident.1
killi_youngandaged_integrated$orig.ident
View(my.deg.new)

###plot scatter plot for different parameters sets
ParameterSetScatterPlot(stable_clusters = stable_clusters,
                        fullsample_idents = fullsample_idents,
                        x_var = "k_param",
                        y_var = "number",
                        facet_rows = "resolution",
                        facet_cols = "pc")

View(my.deg.new)
my.deg.2
my.deg.2.2 <- my.deg.2                                           # Duplicate example data
my.deg.2.2$row_names <- row.names(my.deg.2.2)                     # Apply row.names function
my.deg.2.2  
Idents(killi_youngandaged_integrated) <- killi_youngandaged_integrated$sample
Idents(killi_youngandaged_integrated) <- killi_youngandaged_integrated$CellType

DefaultAssay(my.deg.2) <- "RNA"
my.deg.2@meta.data
keyvals <- ifelse(
  res$log2FoldChange < -2.5, 'royalblue',
  ifelse(res$log2FoldChange > 2.5, 'gold',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'gold'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'

keyvals.colour <- ifelse(
  my.deg.poisson$avg_log2FC < -1.5, 'royalblue',
  ifelse(my.deg.poisson$avg_log2FC > 1.5, 'gold',
         'black'))
keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'gold'] <- 'high'
names(keyvals.colour)[keyvals.colour == 'black'] <- 'mid'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'low'

EnhancedVolcano(my.deg.poisson,rownames(my.deg.poisson),
                x ="avg_log2FC",y ="p_val",
                selectLab = c('CD68','CCL-C5A','ISG15 (2 OF 2)',
                              'CD74','B2M (1 OF 2)','APOEB','KRT18','SFTPB','EPD','HSPA8.1','NFIXA','CELF2','FOSB','FOSL1','ELAVL3','DNAJB2','LOC107381735','HSP90AA1.2','HSPA1A'),
                xlab = bquote(~Log[2]~ 'fold change'),shape = c(6, 4, 2, 11),
                title = 'Custom shape & colour over-ride',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.5,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',col = c('grey30', 'forestgreen', 'royalblue', 'red2'))
head(my.deg.poisson)
EnhancedVolcano(my.deg.new,rownames(my.deg.new),
                x ="avg_log2FC", 
                y ="p_val",
                FCcutoff = 1,
                pCutoff = 0.05,
                pointSize = 2.0,
                labSize = 2.0,
                selectLab = c('CD68','CCL-C5A','ISG15 (2 OF 2)',
                              'CD74','B2M (1 OF 2)','APOEB','KRT18','SFTPB','EPD','HSPA8.1','NFIXA','CELF2','FOSB','FOSL1','ELAVL3','DNAJB2','LOC107381735','HSP90AA1.2','HSPA1A'),
                subtitle = "Young vs Aged",
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1,
                caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
                legendPosition = "right",
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5) + ggplot2::coord_cartesian(xlim=c(-4, 4)) + ggplot2::scale_x_continuous(breaks=seq(-4,4, 1))
my.deg.new$gene <- rownames(my.deg.new)
keyvals.shape <- ifelse(
  my.deg.new$gene %in% i, 'gold',
  'black')
keyvals.shape
keyvals.shape[is.na(keyvals.shape)] <- 'black'
names(keyvals.shape)[keyvals.shape == 'black'] <- 'Others'
names(keyvals.shape)[keyvals.shape == 'gold'] <- 'CD68'
?EnhancedVolcano
View(my.deg.new)
as.numeric(my.deg.near$logFC)
my.deg.near$logFC
str(my.deg.near)
EnhancedVolcano(my.deg.near,my.deg.near$ID,
                x = "logFC",
                y = 'P.Value',
                title = NULL,
                selectLab = c('CD68','CCL-C5A','ISG15 (2 OF 2)','CD74','B2M (1 OF 2)','APOEB','KRT18','SFTPB','EPD','HSPA8.1','NFIXA','CELF2','FOSB','FOSL1','ELAVL3','DNAJB2','LOC107381735','HSP90AA1.2','HSPA1A'),
                subtitle = "Young vs Aged",
                pCutoff = 0.05,boxedLabels = TRUE,
                FCcutoff = 1.2,xlim=c(-6,5.1),labCol = 'black',labSize = 2.0,
                labFace = 'bold',
                pointSize = 4,
                caption = NULL,
                legendPosition = "top",  )

#######usingwilcoxtest
keyvals <- ifelse(my.deg.near$logFC < -1.5, 'orangered2', ifelse(my.deg.near$logFC > 1.5, 'red4', 'grey50'))
keyvals[is.na(keyvals)] <- 'grey50'
names(keyvals)[keyvals == 'orangered2'] <- 'Down-Regulated'
names(keyvals)[keyvals == 'grey50'] <- 'NS'
names(keyvals)[keyvals == 'red4'] <- 'Up-Regulated'
EnhancedVolcano(my.deg.near, lab = my.deg.near$ID, x = 'logFC', y = 'P.Value', selectLab =my.deg.near$ID[which(names(keyvals) %in% c('Down-Regulated', 'Up-Regulated'))], xlab = bquote(~log[2]~ 'fold change'), title = 'G-G+ vs G-G-', pCutoff = 0.05, subtitle = 'Cutoff values (dashed line) at padj=0.05 log2FC=2', FCcutoff = 2, pointSize = 3.0, labSize = 4, shape = c(0, 16, 1, 17), colCustom = keyvals, colAlpha = 1, legendPosition = 'right', legendLabSize = 13, legendIconSize = 5.0, drawConnectors = TRUE, widthConnectors = 0.75, gridlines.major = FALSE, gridlines.minor = FALSE, border = 'partial', borderWidth = 1.5, borderColour = 'black')
EnhancedVolcano(my.deg.poisson,rownames(my.deg.poisson),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = NULL,
                selectLab = c('CD68','CCL-C5A','ISG15 (2 OF 2)','CD74','B2M (1 OF 2)','APOEB','KRT18','SFTPB','EPD','HSPA8.1','NFIXA','CELF2','FOSB','FOSL1','ELAVL3','DNAJB2','LOC107381735','HSP90AA1.2','HSPA1A'),
                subtitle = "Young vs Aged",
                pCutoff = 0.1,boxedLabels = TRUE,
                FCcutoff = 1.0,labCol = 'black',labSize = 2.0,xlim=c(-6,6),
                labFace = 'bold',
                pointSize = 2,
                legendPosition = "top")
                
dev.new(width = , height = y)

dev.new(width = x, height = y)
EnhancedVolcano(my.deg.near,my.deg.near$ID,
                x = "logFC",
                y = 'P.Value',
                selectLab = c('CD68','CCL-C5A','ISG15 (2 OF 2)',
                              'CD74','B2M (1 OF 2)','APOEB','KRT18','SFTPB','EPD','HSPA8.1','NFIXA','CELF2','FOSB','FOSL1','ELAVL3','DNAJB2','LOC107381735','HSP90AA1.2','HSPA1A'),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-14,
                FCcutoff = 0.0,
                pointSize = 3.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,xlim=c(-6,5.1),
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')
?EnhancedVolcano
EnhancedVolcano(my.deg.near, lab = my.deg.near$ID,
                x = "logFC",
                y = 'P.Value', selectLab = c('CD68','CCL-C5A','ISG15 (2 OF 2)',
                                             'CD74','B2M (1 OF 2)','APOEB','KRT18','SFTPB','EPD','HSPA8.1','NFIXA','CELF2','FOSB','FOSL1','ELAVL3','DNAJB2','LOC107381735','HSP90AA1.2','HSPA1A'),
                title = NULL,
                subtitle = "Young vs Aged",
                pCutoff = 0.05,
                FCcutoff = 0,
                labSize = 4,
                colAlpha = 1,
                pointSize = 3,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'))

head(my.deg.2)
str(my.deg.2)
my.deg.2 <- as.list(my.deg.2)

theme_void()
AllyoungIDs <- subset(onlyrevelantPCs, idents=c("Astro-RG1_YoungK","Astro-RG2_YoungK","EPD-RG4_YoungK","NE-RG3_YoungK","NGP.1_YoungK","NGP.2_YoungK","Intercell.1_YoungK","Intercell.2_YoungK","Intercell.3_YoungK","Intercell.4_YoungK"))
AllagedIDs <- subset(onlyrevelantPCs, idents=c("Astro-RG1_AgedK","Astro-RG2_AgedK","EPD-RG4_AgedK","NE-RG3_AgedK","NGP.1_AgedK","NGP.2_AgedK","Intercell.1_AgedK","Intercell.2_AgedK","Intercell.3_AgedK","Intercell.4_AgedK"))

# Make the plot
ggplot(freq_aged, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

freq_aged <- table(Idents(AllagedIDs))
freq_young <- table(Idents(AllyoungIDs))
freq_young
table(Idents(onlyrevelantPCs))
dataPCs <- data.frame(
    category=c("Astro-RG2",	"Astro-RG1",	"Intercell.1",	"Intercell.2",	"Intercell.3",	"Intercell.4",	"NGP.1",	"NGP.2", "EPD-RG4",	"NE-RG3"),
    count=c(136,	9,	348,	501,	112,	5,	116,	76,	173,	134)
  )
table(Idents(AllyoungIDs))
table(Idents(AllagedIDs))
dataPCs$category
dataPCs$count
dataPCs$PCT 
dataPCs$fraction <- dataPCs$count / sum(dataPCs$count)
dataPCs$PCT <- round(dataPCs$fraction * 100, digits =1)

table(Idents(PCs))
dataPCs$ymax = cumsum(dataPCs$fraction)

# Compute the bottom of each rectangle
dataPCs$ymin = c(0, head(dataPCs$ymax, n=-1))

# Compute label position
dataPCs$labelPosition <- (dataPCs$ymax + dataPCs$ymin) / 2
dataPCs$labelPosition
# Compute a good label
dataPCs$label1 <- paste0(c(dataPCs$PCT), "%")
dataPCs$label2 <- paste0(c(dataPCs$category), "%")

pie(dataPCs$PCT, col=col_new_ridgeplot_aged, labels =dataPCs$category )
pie(dataPCs_young$PCT, col=col_new_ridgeplot_aged, labels =dataPCs_young$category )
(pie1,pie2)
?plot_grid
# Make the plot
ggplot(dataPCs, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label2), size=6) +
  scale_fill_manual(values=col_new_ridgeplot_aged) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

col_new_ridgeplot_aged = c("powderblue","#FA340D","#F4BB44","#FFFF8F","#FDDA0D","#E4D00A","#FAB693", "gray","darkcyan","lightsteelblue")
table(Idents (dataPCs))
dataPCs_young <- data.frame(
  category=c("Astro-RG2",	"Astro-RG1",	"Intercell.1",	"Intercell.2",	"Intercell.3",	"Intercell.4",	"NGP.1",	"NGP.2", "EPD-RG4",	"NE-RG3"),
  count=c(38,	158,	270,	80,	260,	217,	139,	117,	74,	74)
)

dataPCs_young$fraction <- dataPCs_young$count / sum(dataPCs_young$count)
dataPCs_young$PCT <- round(dataPCs_young$fraction * 100, digits =1)

dataPCs_young$ymax = cumsum(dataPCs_young$fraction)

# Compute the bottom of each rectangle
dataPCs_young$ymin = c(0, head(dataPCs_young$ymax, n=-1))

# Compute label position
dataPCs_young$labelPosition <- (dataPCs_young$ymax + dataPCs_young$ymin) / 2

# Compute a good label
dataPCs_young$label <- paste0(dataPCs_young$PCT, "%")
dataPCs_young
# Make the plot
ggplot(dataPCs_young, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values=col_new_ridgeplot_aged) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

col_new_ridgeplot_aged = c("powderblue","#FA340D","#F4BB44","#FFFF8F","#FAB693","#E4D00A","#FDDA0D", "gray","darkcyan","steelblue")


pie(table(Idents(killi_youngandaged_integrated_PCs)))
?pie
# Compute percentages
onlyrevelantPCs$fraction <- onlyrevelantPCs$count / sum(onlyrevelantPCs$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, "\n value: ", data$count)

# Make the plot
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

count_table_fresh <- table(killi_youngandaged_integrated_final@meta.data$seurat_clusters, killi_youngandaged_integrated_final@meta.data$sample)
killi_youngandaged_integrated_final@meta.data$seurat_clusters
plot_integrated_clusters(onlyrevelantPCs)
count_table_fresh
View(my.deg.new)
EnhancedVolcano(my.deg.new,lab=rownames(my.deg.new),
                x = 'avg_log2FC',
                y = 'p_val',
                xlab = bquote(~Log[2]~ 'fold change'),
                FCcutoff = 0.1,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                widthConnectors = 0.75) +coord_cartesian(-log10(10E-50), -log10(10E-60))

dataPCs <- data.frame(
  category=c("A", "B", "C"),
  count=c(10, 60, 30)
)

# pseudo-bulk by per donor per cell type
killi.bulk <- AverageExpression(killi_youngandaged_integrated_final, return.seurat = TRUE, group.by = c('condition','seurat_clusters'))
killi.bulk@meta.data
killi.bulk <- NormalizeData(killi.bulk)%>% FindVariableFeatures() %>% ScaleData()%>% RunPCA(npcs = 10)
killi.bulk$seurat_clusters <- rownames(killi.bulk@meta.data)
DimPlot(killi.bulk, group.by = 'seurat_clusters', label = TRUE,  reduction = 'pca')
str(killi.bulk)
View(killi_youngandaged_integrated_final@meta.data)
killi_youngandaged_integrated_final[[1]]@meta.data$condition <- "YoungK"
killi_youngandaged_integrated_final[[2]]@meta.data$condition <- "YoungK"
killi_youngandaged_integrated_final[[3]]@meta.data$treatment <- "AgedK"
View(my.deg.new)
write.xlsx(my.deg.new, file ="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/degexpression_pseudobulk_acrossaged.xlsx")
write.table(my.deg.2, path ="C:/Users/u0129074/Documents/ScRNAseq_paper/Aging paper/Extraanalyses/withsample2_young/degexpression_acrossaged.xls")
my.deg.near <- read.csv(file="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/allcells/agingcombined/allcells_genealtered.csv")
rm(my.deg.2)
?read.csv
my.deg.PCs <- AggregateExpression(killi_youngandaged_integrated_PCs, assay="RNA", group.by="condition", slot='data')

class(my.deg.2)
EnhancedVolcano(my.deg.2, 
                rownames(my.deg.2),
               x ="avg_log2FC", 
              y ="p_val_adj")

str(my.deg.2)
table(Idents(killi_youngandaged_integrated_final))



killi_youngandaged_integrated$nCount_SCT

    condition.diffgenes_MG <- FindMarkers(killi_youngandaged_integrated_final, ident.1 = "MG_AgedK", ident.2= "MG_YoungK", min.pct=0.25, logfc.threshold=0.25)
    write.csv(condition.diffgenes_MG, file="MG.csv")

    df_allfiles <- list.files(path='C:/Users/u0129074/Documents/Rstuff/killifish_sc') %>% 
      lapply(read_csv) %>% 
      bind_rows
str(my.deg.2)
boxplot(avg_log2FC~clusternumber, data=merged.table.1)
str(filteredgenes)
#stripchart(z ~ g, add = TRUE, vertical = TRUE,
           #method = "jitter", col = 3:4, pch = 19)

stripchart(avg_log2FC~clusternumber, add = F,
           data=onlytopgenes_allcells, method = "jitter", jitter=0.4, col="maroon",
           main="Differential expression per cell type",
           xlab="Celltype", vertical=T,
           ylab="Expression",
           pch='a') + RotatedAxis()
?stripchart
View(filteredgenes_new)
group.names=c("IhmN1","IhmN2","IhmN3","ImN1","ImN2","ImN3","ExmN1","ExmN2","ExmN3","OD","Prog","Astro-RG2","NE-RG3","Vas1","Vas2","Vas3")

group.names=c("IhmN1","IhmN2","IhmN3","ImN1","ImN2","ImN3","ExmN1","ExmN2","ExmN3","OD","Prog","Astro-RG2","NE-RG3","Vas1","Vas2","Vas3")
c("ImN1","Ex-mN1","Vas2","EPD","NA","NE-RG3","Astro-RG1","OD","NGP.2","Astro-RG2","Vas3","Ex-mN2","ImN2","Ih-mN1","MG","Intercell","Ex-mN3","NGP.1","Vas1")
data <- data.frame(
  x=PCs_filteredgenes$avg_log2FC,
  y=c("Astro-RG2","Astro-RG1","NE-RG3","NGP1","NGP2","EPD-RG4","Intercell.1","Intercell.2","Intercell.3","Intercell.4"))
filteredgenes$Cluster

# Change baseline
ggplot(data, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=1, yend=y), color="grey") +
  geom_point( color="orange", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) + coord_flip()
  xlab("") +
  ylab("Value of Y")

ifelse(filteredgenes$Topgene=="TRUE", "red", "blue")

youngandaged_isoseq_combinedclusters$Topsgenes <- youngandaged_isoseq_combinedclusters$avg_log2FC > 2.0 || youngandaged_isoseq_combinedclusters$avg_log2FC < -1.5
newgenes <- youngandaged_isoseq_combinedclusters$avg_log2FC > 2.0

filteredgenes$Topgenesonly <- filteredgenes[grepl("TRUE", filteredgenes$Topgene),0 ]
##filteredgenes$Topgene  "TRUE"
View(filteredgenes$Topgenesonly)
View(highlight_df)
  tail(youngandaged_isoseq_combinedclusters$Topgenes)
  #geom_point(data=highlight_df, 
  #           aes(x=avg_log2FC,y=cluster), 
   #          color='red',
    #         size=3)

View(youngandaged_isoseq_combinedclusters$Topgenes)

boxplot(avg_log2FC~Cluster,data=filteredgenes, add=T)
  geom_point(data=filteredgenes, 
             aes(x=avg_log2FC,y=Cluster), 
             color='red',
             size=3)
killi_youngandaged_integrated
library(glmGamPoi)

# filter dataframe to get data to be highligheted
highlight_df <- youngandaged_isoseq_combinedclusters %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)
genes <- highlight_df$gene
?boxplot
text(young_aged_top30_markers$gene, 1.1, labels=young_aged_top30_markers$gene)
youngandaged_isoseq_combinedclusters

data<-replace(youngandaged_isoseq_combinedclusters$avg_log2FC, youngandaged_isoseq_combinedclusters$avg_log2FC>2, "TRUE")
youngandaged_isoseq_combinedclusters$Topgenes[1:455] <- "TRUE"
fit.treat <- treat(youngandaged_isoseq_combinedclusters,lfc=1)
res.treat <- decideTests(fit.treat)
summary(res.treat)

?goana
head(youngandaged_isoseq_combinedclusters)



?stripchart
realcounts <- dataset$nCount_SCT
clusters_new <- dataset$seurat_clusters
dataset <- as.data.frame(killi_youngandaged_integrated_final@meta.data)
FeaturePlot(killi_youngandaged_integrated_final, pt.size=3,features="nGene",reduction="tsne")
head(dataset)
expr_info <- as.data.frame(killi_youngandaged_integrated_final@assays$integrated@counts)
celltype_info <- as.data.frame(killi_youngandaged_integrated_final@meta.data$seurat_clusters)
clusters_new
######Striplot#####
library(cowplot)
theme_set(theme_cowplot())

mg_yandA <- subset(killi_youngandaged_integrated, idents=c("MG"))
Idents(mg_yandA) <- "sample"
avgexpr_mg <- as.data.frame(log1p(AverageExpression(mg_yandA, verbose = FALSE)$RNA))
avgexpr_mg$gene <- rownames(avgexpr_mg)
p1<- ggplot(avgexpr_mg, aes(YoungK, AgedK)) + geom_point() + ggtitle("MG cells")
p1 <- LabelPoints(plot = p1, points = c("APOE","CLU","AIF1"), repel = TRUE)
p1
avgexpr_mg
BiocManager::install("clusterProfiler")
library(clusterProfiler)
clusterProfiler

class(expr_info)
expr_info
?stripchart
killi_youngandaged_integrated$CellType
