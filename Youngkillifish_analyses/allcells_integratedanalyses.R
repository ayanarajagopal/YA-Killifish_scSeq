library(Seurat)
library(SeuratDisk)
library(DelayedArray)
library(ggplot2)
library(dplyr)
library(Matrix)
library(cowplot)
library(scater)
library(glue)
library(MAST)
library(SingleCellExperiment)
library(ggthemes)

raw.data_S1_S2 = list("Sample1"=youngisoseq,"Sample2"=youngisoseq_S2)
youngisoseq$sample <- "Sample1"
youngisoseq_S2$sample <- "Sample2"

dim(youngisoseq)
dim(youngisoseq_S2)
dim(agedisoseq)

combined.features_Y1_Y2 <- SelectIntegrationFeatures(object.list = raw.data_S1_S2, nfeatures = 3000)
raw.data_S1_S2 <- PrepSCTIntegration(object.list = raw.data_S1_S2, anchor.features = combined.features_Y1_Y2)

DefaultAssay(youngisoseq) <- "SCT"
DefaultAssay(youngisoseq_S2) <- "SCT"

int.anchors_Y1_Y2 <- FindIntegrationAnchors(object.list = raw.data_S1_S2, 
                                            anchor.features = combined.features_Y1_Y2, reduction="cca", normalization.method = "SCT")

View(int.anchors_Y1_Y2@anchor.features)
VlnPlot(killi_young_integrated, features = c("percent.mt","percent.ribo", "nGene","nUMI"), group.by="orig.ident",combine=T, stack=T, flip=T, pt.size = 1)
FeatureScatter(object = killi_young_integrated,group.by="orig.ident", feature1 = "nUMI", feature2 = "nGene") +NoLegend()

killi_young_integrated <- IntegrateData(anchorset = int.anchors_Y1_Y2, normalization.method = "SCT", dims = 1:30)

DefaultAssay(killi_young_integrated) <- "integrated"

killi_young_integrated <- ScaleData(killi_young_integrated, verbose = FALSE)
killi_young_integrated <- RunPCA(object = killi_young_integrated, verbose = FALSE)
ElbowPlot(killi_young_integrated, ndims=40) + theme_clean()

# Determine the K-nearest neighbor graph
killi_young_integrated <- FindNeighbors(object = killi_young_integrated, 
                                        dims = 1:30)

# Determine the clusters for various resolutions                                
killi_young_integrated <- FindClusters(object = killi_young_integrated,
                                       resolution = 0.6)

#Calc. of TSNE and perplexity optimization
killi_young_integrated <- RunTSNE(object=killi_young_integrated, perplexity=35, reduction="pca", dims= 1:20)
killi_young_integrated <- RunUMAP(object = killi_young_integrated, reduction= "pca", dims = 1:20)

DimPlot(killi_young_integrated, split.by = "sample", reduction = "tsne", pt.size = 1.8) +NoLegend()
DimPlot(killi_young_integrated, reduction = "tsne", pt.size = 0.6, group.by ="nFeature_RNA",label.size=4)
VlnPlot(killi_young_integrated, stack=T, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.5) + 

killi_young_integrated@meta.data
plot_grid(p1, p2)
dev.off()
dim(killi_young_integrated)

#Normalize and setting RNA assay
DefaultAssay(killi_young_integrated) <- "RNA"
killi_young_integrated <- NormalizeData(killi_young_integrated, verbose = T)

#####Find all markers #####
youngall_isoseq_combinedclusters <- FindAllMarkers(killi_young_integrated, min.pct = 0.25, only.pos = F, logfc.threshold = 0.25)
head(youngall_isoseq_combinedclusters)

youngall_isoseq_combinedclusters_200 <- youngall_isoseq_combinedclusters %>% group_by(cluster) %>% top_n(n=200, wt = avg_log2FC)
youngall_isoseq_combinedclusters_100 <- youngall_isoseq_combinedclusters %>% group_by(cluster) %>% top_n(n=100, wt = avg_log2FC)
youngall_isoseq_combinedclusters_10 <- youngall_isoseq_combinedclusters %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)

young_aged_top30_markers <- youngandaged_isoseq_combinedclusters %>% group_by(cluster) %>% top_n(n=30, wt = avg_log2FC)

head(young_aged_top30_markers)
write.table(young_aged_isoseq_markers, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Aug2021/y_a_markers_res0.6_new.txt", sep="\t", col.names=TRUE)
write_xlsx(youngall_isoseq_combinedclusters_200, col_names = T, path="C:/Users/u0129074/Documents/rstudio_jan2020/20jan_afterintegration/Young_old_separate/Young_Isoseq_latest/July2020/young_aggr_oct2022/Analysis/young_Sample1_2_clusters_res0.6_topmarkers.xlsx")
library(dittoSeq)

genes.to.label.clusters <- c("EPD","APOA1","WNT8B","ZIC2","PDGFRB","KLF2","SOX17","STAB2","PCNA","STMN1A","HMGB2A","CD9B","OLIG1","CX43","FABP7A","GLUL (2 OF 2)","SLC1A2B","APOEB","AIF1","LCP1","SLC17A6B","SV2A","ELAVL3","NEUROD2","MIBP2","TUBB5")

DotPlot(killi_young_integrated, cluster.idents = F, features= genes.to.label.clusters, scale=T, cols=c("green","purple")) + coord_flip() + RotatedAxis() 

dittoScatterPlot(killi_young_integrated, 
                 x.var = "SV2A", y.var="CX43", 
                 assay.x = "RNA", assay.y = "RNA", color.var = "CellType") + scale_color_manual(values = killi_young_integrated$CellType) +
  ggtitle("Scatterplot for CD3/CD20 labelled by celltype")
?dittoScatterPlot
dittoHeatmap(killi_young_integrated, genes=NULL,
             annot.by = c("CellType","sample"), 
             fontsize = 7)
set.seed(9)
DoHeatmap(killi_young_integrated, features = youngall_isoseq_combinedclusters_10$gene,draw.lines = TRUE,
          lines.width = NULL,hjust=0,angle=60, size=5,
          group.bar.height = 0.02) + scale_fill_gradientn(colors = c("black","red"))+ NoLegend()

young_allS.newnames.0.6_new <- c("Intercell.NC", "mN1", "ImN1", "mN2", "mN3", "ImN2", "Inter-cell.PC", "mN4", "MG", "Astro-RG1", "mN5", "mN6", "mN7", "OPC/OD", "NGP", "Vas1", "NA", "Vas2", "mN8", "NE-RG3", "EPD-RG4","Astro-RG2", "Vas3")
names(young_allS.newnames.0.6_new) <- levels(killi_young_integrated)
killi_young_integrated <- RenameIdents(killi_young_integrated, young_allS.newnames.0.6_new)

levels(killi_young_integrated)
names(young_allS.newnames.0.6_new)
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
killi_young_integrated@Idents <- plyr::mapvalues(x = killi_young_integrated@ident, from = current.cluster.ids, to = young_allS.newnames.0.6_new)

killi_young_integrated <- RenameIdents(killi_young_integrated, '0'="Intercell.NC",'1'="Intercell.NC",'2'="ImN1",'3'="mN2",'4'="mN3",'5'="ImN2",'6'="Intercell.PC",'7'="mN4",'8'="MG",'9'="Astro-RG1",'10'= "mN5",'11'= "mN6",'12'= "mN7",'13'= "OPC/OD",'14'= "NGP",'15'= "Vas1",'16'= "NA",'17'= "Vas2",'18'= "mN8",'19'= "NE-RG3",'20'= "EPD-RG4",'21'= "Astro-RG2",'22'= "Vas3")
killi_young_integrated <- RenameIdents(killi_young_integrated, 'Intercell.NC'="NC","mN1"="NC",'ImN1'="NC",'mN2'="NC",'mN3'="NC",'ImN2'="NC",'Intercell.PC'="PC",'mN4'="NC",'MG'="Other",'Astro-RG1'="PC",'mN5'= "NC",'mN6'= "NC",'mN7'= "NC",'OPC/OD'= "Other",'NGP'= "PC",'Vas1'= "Other",'NA'= "Other",'Vas2'= "Other",'mN8'= "NC",'NE-RG3'= "PC",'EPD-RG4'= "PC",'Astro-RG2'= "PC",'Vas3'= "Other")

killi_young_integrated$CellType <- Idents(killi_young_integrated)
killi_young_integrated$CellType
killi_young_integrated$seurat_clusters <- Idents(killi_young_integrated)
killi_young_integrated$seurat_clusters
table(Idents(killi_young_integrated$SCT_snn_res.0.6))

##Plot joint density of PCNA and glul genes; glial/nonglial signature##
plot_density(killi_young_integrated, joint=F,c("SOX2"), reduction="tsne")
plot_density
Plot_Density_Custom(
  killi_young_integrated,
  features=c("PCNA","GLUL (2 OF 2)"),
  joint = FALSE,
  viridis_palette = "magma",
  custom_palette = NULL,
  pt.size = 1,
  aspect_ratio = NULL,
  reduction = NULL,
  combine = TRUE)
###cellcyclescoring
killi_young_integrated <- CellCycleScoring(killi_young_integrated, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
?CellCycleScoring
PCs_integrated[[]]
as_tibble(killi_young_integrated[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()
killi_young_integrated@meta.data %>%
  group_by(CellType,Phase) %>%
  count() %>%
  group_by(CellType) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=CellType,y=percent, split,fill=Phase), cols=c("blue","#ffa10a","#00b300")) +
  geom_col() + RotatedAxis()+
  ggtitle("Percentage of cell cycle phases per cluster")
DimPlot(killi_young_integrated, reduction = "tsne",label=, pt.size = 2.5, group.by="Phase") + theme_test()

saveRDS(killi_young_integrated, file = "C:/Users/u0129074/Documents/ScRNAseq_paper/youngkillifish_allcells.rds")